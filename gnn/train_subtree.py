"""
Subtree-level GNN Training Script

Based on train.py (per-conf GNN), modified to predict subtree-level
log Z difference (free energy correction).

Input:  (fixed_rcs, free_mask) — partial assignment
Output: ΔF = -kT × (log Z_CCD - log Z_emat)  in kcal/mol

At inference:  Z_CCD ∈ Z_emat × [exp(-(pred+q)/kT), exp(-(pred-q)/kT)]
where q is the CP-calibrated prediction error quantile.

Training data is generated from confs.csv by grouping conformations
by partial assignments at varying levels of fixed positions.

Usage:
    python train_subtree.py --data_dir gnn_data/2RL0_all20_4pos/complex --epochs 200
"""

import argparse
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from itertools import combinations


NUM_AA_TYPES = 20
UNKNOWN_AA_IDX = NUM_AA_TYPES  # index 20 for non-standard residues (aa_type_idx == -1)
MAX_CHI = 4
KT = 0.5922  # kcal/mol at 298K

# Standard amino acid 3-letter codes → index (must match Java side)
AA_TYPES = {
    "ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4,
    "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9,
    "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14,
    "SER": 15, "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19,
}


# ============================================================
# Data Loading — Subtree samples from confs.csv
# ============================================================

class SubtreeDataset(Dataset):
    """Each sample: (fixed_rcs, free_mask) → free energy correction ΔF = F_CCD - F_emat."""

    def __init__(self, fixed_rcs, free_masks, labels):
        self.fixed_rcs = torch.tensor(fixed_rcs, dtype=torch.long)
        self.free_masks = torch.tensor(free_masks, dtype=torch.bool)
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return self.fixed_rcs[idx], self.free_masks[idx], self.labels[idx]


def generate_subtree_samples(confs, emat_energies, residuals, num_pos,
                             min_confs=10, max_samples_per_k=50000, seed=42):
    """
    Group existing conformations by partial assignments to create subtree samples.

    For each number of fixed positions k (from 1 to num_pos-1):
      - For each combination of k positions to fix:
        - Group confs by their RC values at those positions
        - For each group with >= min_confs members:
          - label = -kT × (logsumexp(-E_CCD/kT) - logsumexp(-E_emat/kT))
            = free energy correction for the subtree

    Returns: fixed_rcs (N, num_pos), free_masks (N, num_pos), labels (N,)
    """
    rng = np.random.RandomState(seed)
    all_fixed_rcs = []
    all_free_masks = []
    all_labels = []

    # Also include full confs (k=num_pos, no free positions) for multi-task
    print("Generating subtree samples from conformations...")

    for k in range(1, num_pos):
        n_free = num_pos - k
        all_combos = list(combinations(range(num_pos), k))

        # Subsample combos if too many
        if len(all_combos) > 100:
            all_combos = [all_combos[i] for i in rng.choice(len(all_combos), 100, replace=False)]

        k_count = 0
        for fixed_pos in all_combos:
            fixed_pos = list(fixed_pos)
            free_pos = [p for p in range(num_pos) if p not in fixed_pos]

            # Group by fixed positions' RC values
            # Build group keys from fixed columns
            keys = confs[:, fixed_pos]  # (N, k)
            # Use structured array for fast groupby
            key_dtype = np.dtype([(f"p{i}", np.int64) for i in range(k)])
            key_arr = np.array([tuple(row) for row in keys], dtype=key_dtype)
            unique_keys, inverse, counts = np.unique(key_arr, return_inverse=True, return_counts=True)

            for gid in range(len(unique_keys)):
                if counts[gid] < min_confs:
                    continue

                mask = inverse == gid
                g_emat = emat_energies[mask]
                g_res = residuals[mask]

                # Exact log Z difference (free energy correction)
                # label = -kT * (logsumexp(-E_CCD/kT) - logsumexp(-E_emat/kT))
                g_ccd = g_emat + g_res
                def _lse(x, kt=KT):
                    x_div = -x / kt
                    x_max = x_div.max()
                    return -(np.log(np.exp(x_div - x_max).sum()) + x_max) * kt
                lse_emat = _lse(g_emat)
                lse_ccd = _lse(g_ccd)
                label = lse_ccd - lse_emat

                # Build the (fixed_rcs, free_mask) pair
                fixed_rcs_row = np.zeros(num_pos, dtype=np.int64)
                free_mask_row = np.ones(num_pos, dtype=bool)
                for i, p in enumerate(fixed_pos):
                    fixed_rcs_row[p] = int(unique_keys[gid][i])
                    free_mask_row[p] = False

                all_fixed_rcs.append(fixed_rcs_row)
                all_free_masks.append(free_mask_row)
                all_labels.append(label)
                k_count += 1

                if k_count >= max_samples_per_k:
                    break
            if k_count >= max_samples_per_k:
                break

        print(f"  fix {k}/{num_pos} positions ({len(all_combos)} combos): "
              f"{k_count} subtree samples")

    # Also add full confs as subtrees with 0 free positions (multi-task)
    # Use ALL conformations — no cap here
    n_full = len(confs)
    for i in range(n_full):
        all_fixed_rcs.append(confs[i].copy())
        all_free_masks.append(np.zeros(num_pos, dtype=bool))
        all_labels.append(float(residuals[i]))

    print(f"  full confs (k={num_pos}): {n_full} samples")

    fixed_rcs = np.array(all_fixed_rcs)
    free_masks = np.array(all_free_masks)
    labels = np.array(all_labels)
    print(f"Total subtree samples: {len(labels)}")
    return fixed_rcs, free_masks, labels


def load_data(data_dir):
    """Load all exported CSV files and generate subtree training samples."""
    import time as _time
    t0 = _time.time()

    print("Loading confs.csv ...", flush=True)
    confs_df = pd.read_csv(os.path.join(data_dir, "confs.csv"))
    print(f"  confs.csv loaded in {_time.time() - t0:.1f}s", flush=True)

    graph_df = pd.read_csv(os.path.join(data_dir, "graph.csv"))
    meta_df = pd.read_csv(os.path.join(data_dir, "meta.csv"))
    rc_feat_df = pd.read_csv(os.path.join(data_dir, "rc_features.csv"))
    pair_df = pd.read_csv(os.path.join(data_dir, "pairwise_energies.csv"))
    ca_dist_df = pd.read_csv(os.path.join(data_dir, "ca_distances.csv"))

    num_pos = meta_df.shape[0]
    rc_cols = [f"rc_{i}" for i in range(num_pos)]
    confs = confs_df[rc_cols].values.astype(np.int64)

    residuals = confs_df["residual"].values.astype(np.float64)
    emat_energies = confs_df["E_emat"].values.astype(np.float64)
    has_rigid = "E_rigid" in confs_df.columns
    rigid_energies = confs_df["E_rigid"].values.astype(np.float64) if has_rigid else None

    del confs_df

    # Filter (same as train.py)
    n_raw = len(confs)
    mask = np.ones(n_raw, dtype=bool)
    emat_cap = -20.0
    mask &= emat_energies <= emat_cap
    residual_cap = 100.0
    mask &= np.abs(residuals) < residual_cap
    if has_rigid:
        mask &= rigid_energies <= 0.0

    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    print(f"Filtered: {n_raw} -> {len(confs)} conformations")

    # Build graph tables (same as train.py)
    edge_index = torch.tensor(graph_df[["src", "dst"]].values.T, dtype=torch.long)
    num_rcs = meta_df["num_rcs"].values.tolist()
    max_rcs = max(num_rcs)

    aa_table = torch.zeros(num_pos, max_rcs, dtype=torch.long)
    chi_table = torch.zeros(num_pos, max_rcs, MAX_CHI, dtype=torch.float32)

    rc_feat_df = rc_feat_df.sort_values(["pos", "rc"])
    for p in range(num_pos):
        pos_mask = rc_feat_df["pos"].values == p
        pos_data = rc_feat_df[pos_mask]
        n = len(pos_data)
        aa_vals = pos_data["aa_type_idx"].values
        aa_vals = np.where(aa_vals < 0, UNKNOWN_AA_IDX, aa_vals)
        aa_table[p, :n] = torch.tensor(aa_vals, dtype=torch.long)
        chis = pos_data[["chi1", "chi2", "chi3", "chi4"]].fillna(0.0).values
        chi_table[p, :n] = torch.tensor(chis, dtype=torch.float32)

    src_edges, dst_edges = edge_index
    num_edges = src_edges.shape[0]
    pair_table = torch.zeros(num_edges, max_rcs * max_rcs, dtype=torch.float32)

    edge_map = {}
    for e in range(num_edges):
        edge_map[(src_edges[e].item(), dst_edges[e].item())] = e

    p1 = pair_df["pos1"].values
    r1 = pair_df["rc1"].values
    p2 = pair_df["pos2"].values
    r2 = pair_df["rc2"].values
    evals = pair_df["E_pair_min"].values
    for i in range(len(evals)):
        e_idx = edge_map.get((p1[i], p2[i]))
        if e_idx is not None:
            pair_table[e_idx, r1[i] * max_rcs + r2[i]] = evals[i]
        e_idx = edge_map.get((p2[i], p1[i]))
        if e_idx is not None:
            pair_table[e_idx, r2[i] * max_rcs + r1[i]] = evals[i]

    ca_p1 = ca_dist_df["pos1"].values
    ca_p2 = ca_dist_df["pos2"].values
    ca_d = ca_dist_df["distance"].values
    ca_dist_map = {}
    for i in range(len(ca_d)):
        ca_dist_map[(ca_p1[i], ca_p2[i])] = ca_d[i]
        ca_dist_map[(ca_p2[i], ca_p1[i])] = ca_d[i]

    ca_dist_vec = torch.zeros(num_edges, dtype=torch.float32)
    for e in range(num_edges):
        s, d = src_edges[e].item(), dst_edges[e].item()
        ca_dist_vec[e] = ca_dist_map.get((s, d), 0.0)

    elapsed = _time.time() - t0
    print(f"Loaded {len(confs)} conformations, {num_pos} positions, "
          f"{num_edges} directed edges, max_rcs={max_rcs} in {elapsed:.1f}s")

    return {
        "confs": confs,
        "residuals": residuals,
        "emat_energies": emat_energies,
        "edge_index": edge_index,
        "num_rcs": num_rcs,
        "num_pos": num_pos,
        "max_rcs": max_rcs,
        "aa_table": aa_table,
        "chi_table": chi_table,
        "pair_table": pair_table,
        "ca_dist_vec": ca_dist_vec,
    }


# ============================================================
# GNN Model v2
# ============================================================

class SubtreeGNN(nn.Module):
    """
    Subtree-level GNN that predicts the free energy correction ΔF = F_CCD - F_emat
    from a partial assignment (some positions fixed, others free).

    Architecture identical to InteractionGNNv2 except:
      - forward() takes (fixed_rcs, free_mask) instead of (confs)
      - Fixed positions: lookup specific RC embedding (same as before)
      - Free positions: use precomputed mean embedding over all RCs
      - Extra is_fixed feature (1 dim) concatenated to node input
      - Edge features: mean pairwise energy for edges involving free positions
    """

    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec,
                 aa_embed_dim=16, pos_embed_dim=8, node_dim=32,
                 hidden_dim=64, num_layers=3, dropout=0.1):
        super().__init__()
        self.num_pos = num_pos
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.max_rcs = max_rcs
        self.aa_embed_dim = aa_embed_dim

        edge_feat_dim = 2  # [pair_energy, ca_distance]

        # Shared AA type embedding
        self.aa_embedding = nn.Embedding(NUM_AA_TYPES + 1, aa_embed_dim)

        # Position embedding
        self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        chi_input_dim = MAX_CHI * 2

        # +1 for is_fixed feature
        self.node_proj = nn.Sequential(
            nn.Linear(aa_embed_dim + chi_input_dim + pos_embed_dim + 1, node_dim),
            nn.SiLU(),
        )

        # Register lookup tables as buffers
        self.register_buffer("aa_table", aa_table)          # (num_pos, max_rcs)
        self.register_buffer("chi_table", chi_table)        # (num_pos, max_rcs, 4)
        self.register_buffer("pair_table", pair_table)      # (num_edges, max_rcs^2)
        self.register_buffer("ca_dist_vec", ca_dist_vec)    # (num_edges,)
        self.register_buffer("edge_src", edge_index[0])     # (num_edges,)
        self.register_buffer("edge_dst", edge_index[1])     # (num_edges,)

        # Precompute mean embeddings for free positions
        # mean_aa_feat[p] = mean of aa_embedding over all RCs at position p
        # mean_chi_feat[p] = mean of chi sin/cos encoding over all RCs at position p
        # mean_pair[e] = mean pairwise energy over all RC pairs for edge e
        num_rcs_t = torch.tensor(num_rcs, dtype=torch.long)
        self.register_buffer("num_rcs_t", num_rcs_t)

        # These will be computed after aa_embedding is initialized (in _precompute_free_features)
        self.register_buffer("mean_chi_feat", torch.zeros(num_pos, chi_input_dim))

        # Clamp inf in pair_table (steric clashes) to avoid NaN
        PAIR_CLAMP = 100.0  # kcal/mol
        pair_table.clamp_(max=PAIR_CLAMP)

        # Precompute mean chi features (doesn't depend on learned params)
        for p in range(num_pos):
            nr = num_rcs[p]
            chi_raw = chi_table[p, :nr]  # (nr, 4)
            chi_rad = chi_raw * (3.14159265358979 / 180.0)
            chi_sc = torch.cat([torch.sin(chi_rad), torch.cos(chi_rad)], dim=-1)  # (nr, 8)
            self.mean_chi_feat[p] = chi_sc.mean(dim=0)

        # Message passing layers
        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()

        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, node_dim),
            ))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, node_dim),
            ))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        # Readout
        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def _get_mean_aa_feat(self):
        """Compute mean AA embedding per position (depends on learned aa_embedding)."""
        # Called during forward, not cached, so gradients flow through aa_embedding
        feats = []
        for p in range(self.num_pos):
            nr = self.num_rcs_t[p].item()
            aa_indices = self.aa_table[p, :nr]  # (nr,)
            aa_embeds = self.aa_embedding(aa_indices)  # (nr, aa_embed_dim)
            feats.append(aa_embeds.mean(dim=0))  # (aa_embed_dim,)
        return torch.stack(feats, dim=0)  # (num_pos, aa_embed_dim)

    def forward(self, fixed_rcs, free_mask):
        """
        Args:
            fixed_rcs: (batch, num_pos) RC indices; free positions can be any value (ignored)
            free_mask:  (batch, num_pos) bool; True = free position
        Returns:
            delta_F: (batch,) predicted free energy correction ΔF = F_CCD - F_emat (kcal/mol)
        """
        batch = fixed_rcs.shape[0]
        num_pos = self.num_pos
        node_dim = self.node_dim
        device = fixed_rcs.device

        # === Node features ===
        # Fixed positions: lookup specific RC
        aa_idx_fixed = torch.gather(
            self.aa_table.unsqueeze(0).expand(batch, -1, -1),
            2, fixed_rcs.unsqueeze(-1)
        ).squeeze(-1)  # (batch, num_pos)
        aa_feat_fixed = self.aa_embedding(aa_idx_fixed)  # (batch, num_pos, aa_embed_dim)

        chi_raw_fixed = torch.gather(
            self.chi_table.unsqueeze(0).expand(batch, -1, -1, -1),
            2, fixed_rcs.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 1, MAX_CHI)
        ).squeeze(2)  # (batch, num_pos, 4)
        chi_rad_fixed = chi_raw_fixed * (3.14159265358979 / 180.0)
        chi_feat_fixed = torch.cat([torch.sin(chi_rad_fixed), torch.cos(chi_rad_fixed)], dim=-1)

        # Free positions: mean embeddings
        mean_aa_feat = self._get_mean_aa_feat()  # (num_pos, aa_embed_dim)
        aa_feat_free = mean_aa_feat.unsqueeze(0).expand(batch, -1, -1)
        chi_feat_free = self.mean_chi_feat.unsqueeze(0).expand(batch, -1, -1)

        # Mix by free_mask
        mask_3d = free_mask.unsqueeze(-1)  # (batch, num_pos, 1)
        aa_feat = torch.where(mask_3d, aa_feat_free, aa_feat_fixed)
        chi_feat = torch.where(mask_3d, chi_feat_free, chi_feat_fixed)

        # is_fixed feature
        is_fixed = (~free_mask).float().unsqueeze(-1)  # (batch, num_pos, 1)

        # Position embedding
        pos_idx = torch.arange(num_pos, device=device)
        pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)

        # Combine and project
        node_input = torch.cat([aa_feat, chi_feat, pos_feat, is_fixed], dim=-1)
        h = self.node_proj(node_input)  # (batch, num_pos, node_dim)

        # === Edge features ===
        src = self.edge_src
        dst = self.edge_dst
        num_edges = src.shape[0]

        src_exp = src.unsqueeze(0).expand(batch, -1)
        dst_exp = dst.unsqueeze(0).expand(batch, -1)

        # Check if either endpoint is free
        src_free = torch.gather(free_mask, 1, src_exp)  # (batch, num_edges)
        dst_free = torch.gather(free_mask, 1, dst_exp)
        either_free = src_free | dst_free  # (batch, num_edges)

        # Fixed-fixed: lookup specific pair energy
        rc_src = torch.gather(fixed_rcs, 1, src_exp)
        rc_dst = torch.gather(fixed_rcs, 1, dst_exp)
        flat_idx = rc_src * self.max_rcs + rc_dst
        # Advanced indexing avoids materializing [batch, num_edges, max_rcs**2]
        # (that intermediate hits CUDA grid/memory limits on big confspaces — e.g. for
        # max_rcs=189 and num_edges=82 at batch=1149 it is ~13 GB and triggers
        # cudaErrorInvalidConfiguration in the ONNX Expand kernel).
        edge_arange = torch.arange(num_edges, device=device)
        pair_energy_fixed = self.pair_table[edge_arange, flat_idx]  # (batch, num_edges)

        # Free edge: use 0 (unknown pairwise interaction — model infers from node features)
        pair_energy = torch.where(either_free, torch.zeros_like(pair_energy_fixed), pair_energy_fixed)

        ca_dist = self.ca_dist_vec.unsqueeze(0).expand(batch, -1)
        edge_feats = torch.stack([pair_energy, ca_dist], dim=-1)

        # === Message passing (identical to v2) ===
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)
            h_dst = torch.gather(h, 1, dst_idx)

            msg_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)
            messages = self.message_mlps[layer](msg_input)

            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)

            update_input = torch.cat([h, agg], dim=-1)
            h = h + self.update_mlps[layer](update_input)
            h = self.layer_norms[layer](h)

        # === Readout ===
        node_energy = self.node_readout(h).squeeze(-1)
        total_node = node_energy.sum(dim=1)

        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)
        edge_energy = self.edge_readout(edge_input).squeeze(-1)
        total_edge = edge_energy.sum(dim=1) / 2.0

        return total_node + total_edge


# ============================================================
# Training
# ============================================================

def train(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    num_gpus = torch.cuda.device_count()
    print(f"Device: {device}, GPUs available: {num_gpus}")

    data = load_data(args.data_dir)

    # Generate subtree training samples from conformations
    fixed_rcs, free_masks, labels = generate_subtree_samples(
        data["confs"], data["emat_energies"], data["residuals"],
        data["num_pos"], min_confs=args.min_confs,
        max_samples_per_k=args.max_samples_per_k, seed=args.seed,
    )

    # Train/val split
    n = len(labels)
    perm = np.random.RandomState(args.seed).permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    train_ds = SubtreeDataset(fixed_rcs[train_idx], free_masks[train_idx], labels[train_idx])
    val_ds = SubtreeDataset(fixed_rcs[val_idx], free_masks[val_idx], labels[val_idx])
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)

    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")

    # Label statistics
    print(f"Label stats: mean={labels.mean():.4f}, std={labels.std():.4f}, "
          f"min={labels.min():.4f}, max={labels.max():.4f}")

    # Count by number of free positions
    n_free_counts = {}
    for i in range(len(free_masks)):
        nf = int(free_masks[i].sum())
        n_free_counts[nf] = n_free_counts.get(nf, 0) + 1
    print("Samples by #free positions:", dict(sorted(n_free_counts.items())))

    model = SubtreeGNN(
        num_pos=data["num_pos"],
        max_rcs=data["max_rcs"],
        num_rcs=data["num_rcs"],
        aa_table=data["aa_table"],
        chi_table=data["chi_table"],
        pair_table=data["pair_table"],
        edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        aa_embed_dim=args.aa_embed_dim,
        pos_embed_dim=args.pos_embed_dim,
        node_dim=args.node_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
    ).to(device)

    num_params = sum(p.numel() for p in model.parameters())
    num_buffer_elems = sum(b.numel() for b in model.buffers())
    print(f"Model parameters: {num_params:,}")
    print(f"Buffer elements: {num_buffer_elems:,}")

    if num_gpus > 1:
        print(f"Using DataParallel on {num_gpus} GPUs")
        model = nn.DataParallel(model)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_val_rmse = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        train_mse = 0.0
        for fixed_batch, mask_batch, target_batch in train_loader:
            fixed_batch = fixed_batch.to(device)
            mask_batch = mask_batch.to(device)
            target_batch = target_batch.to(device)

            pred = model(fixed_batch, mask_batch)
            if args.loss == "huber":
                loss = nn.functional.huber_loss(pred, target_batch, delta=args.huber_delta)
            else:
                loss = nn.functional.mse_loss(pred, target_batch)
            train_mse += nn.functional.mse_loss(pred, target_batch).item() * len(fixed_batch)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

        train_rmse = (train_mse / len(train_ds)) ** 0.5
        scheduler.step()

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for fixed_batch, mask_batch, target_batch in val_loader:
                fixed_batch = fixed_batch.to(device)
                mask_batch = mask_batch.to(device)
                target_batch = target_batch.to(device)
                pred = model(fixed_batch, mask_batch)
                val_loss += nn.functional.mse_loss(pred, target_batch).item() * len(fixed_batch)
        val_rmse = (val_loss / len(val_ds)) ** 0.5

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            raw_sd = model.module.state_dict() if isinstance(model, nn.DataParallel) else model.state_dict()
            best_state = {k: v.cpu().clone() for k, v in raw_sd.items()}

        if epoch % args.log_every == 0 or epoch == 1:
            print(f"Epoch {epoch:4d}  train_rmse={train_rmse:.4f}  "
                  f"val_rmse={val_rmse:.4f}  best={best_val_rmse:.4f}  "
                  f"lr={scheduler.get_last_lr()[0]:.2e}", flush=True)

        if epoch % 50 == 0:
            ckpt_dir = os.path.join(args.data_dir, "model_subtree")
            os.makedirs(ckpt_dir, exist_ok=True)
            raw_sd = model.module.state_dict() if isinstance(model, nn.DataParallel) else model.state_dict()
            torch.save({
                "model_state": {k: v.cpu().clone() for k, v in raw_sd.items()},
                "optimizer_state": optimizer.state_dict(),
                "scheduler_state": scheduler.state_dict(),
                "epoch": epoch,
                "num_rcs": data["num_rcs"],
                "num_pos": data["num_pos"],
                "max_rcs": data["max_rcs"],
                "args": vars(args),
                "val_rmse": best_val_rmse,
            }, os.path.join(ckpt_dir, "subtree_checkpoint.pt"))

    # Restore best
    base_model = model.module if isinstance(model, nn.DataParallel) else model
    base_model.load_state_dict(best_state)
    base_model.eval()
    model = base_model

    # === Final evaluation by number of free positions ===
    print("\n" + "=" * 70)
    print("Subtree GNN — Validation Results")
    print("=" * 70)

    all_preds, all_targets, all_nfree = [], [], []
    with torch.no_grad():
        for fixed_batch, mask_batch, target_batch in val_loader:
            fixed_batch = fixed_batch.to(device)
            mask_batch = mask_batch.to(device)
            pred = model(fixed_batch, mask_batch)
            all_preds.append(pred.cpu().numpy())
            all_targets.append(target_batch.numpy())
            all_nfree.append(mask_batch.sum(dim=1).cpu().numpy())

    all_preds = np.concatenate(all_preds)
    all_targets = np.concatenate(all_targets)
    all_nfree = np.concatenate(all_nfree).astype(int)

    errors = all_preds - all_targets
    print(f"Overall: N={len(errors)}, MAE={np.mean(np.abs(errors)):.4f}, "
          f"RMSE={np.sqrt(np.mean(errors**2)):.4f}, Bias={np.mean(errors):.4f}")

    for nf in sorted(set(all_nfree)):
        mask = all_nfree == nf
        if mask.sum() < 5:
            continue
        e = errors[mask]
        print(f"  #free={nf}: N={mask.sum()}, MAE={np.mean(np.abs(e)):.4f}, "
              f"RMSE={np.sqrt(np.mean(e**2)):.4f}")

    # CP calibration
    val_errors = np.abs(errors)
    n_cal = len(val_errors)
    cp_stats = {}
    print(f"\nConformal Prediction calibration (N={n_cal}):")
    for alpha in [0.10, 0.05, 0.01, 0.005, 0.001]:
        level = min(np.ceil((1 - alpha) * (n_cal + 1)) / n_cal, 1.0)
        q = float(np.quantile(val_errors, level))
        cp_stats[f"q_alpha_{alpha}"] = q
        print(f"  α={alpha:.3f}  q={q:.4f} kcal/mol")
    cp_stats["n_cal"] = n_cal

    # Save
    out_dir = os.path.join(args.data_dir, "model_subtree")
    os.makedirs(out_dir, exist_ok=True)
    torch.save({
        "model_state": best_state,
        "num_rcs": data["num_rcs"],
        "num_pos": data["num_pos"],
        "max_rcs": data["max_rcs"],
        "args": vars(args),
        "val_rmse": best_val_rmse,
        "cp_stats": cp_stats,
    }, os.path.join(out_dir, "subtree_checkpoint.pt"))

    # Export ONNX
    export_onnx(model, data["num_pos"], out_dir, device)

    print(f"\nModel saved to {out_dir}")


def export_onnx(model, num_pos, out_dir, device):
    """Export subtree model to ONNX with two inputs: fixed_rcs and free_mask."""
    # Move model to CPU for ONNX export to avoid device mismatch
    # between registered buffers and constant-folded tensors
    model = model.cpu()
    model.eval()
    dummy_rcs = torch.zeros(1, num_pos, dtype=torch.long)
    dummy_mask = torch.zeros(1, num_pos, dtype=torch.bool)

    onnx_path = os.path.join(out_dir, "subtree_model.onnx")
    torch.onnx.export(
        model,
        (dummy_rcs, dummy_mask),
        onnx_path,
        input_names=["fixed_rcs", "free_mask"],
        output_names=["residual"],
        dynamic_axes={
            "fixed_rcs": {0: "batch"},
            "free_mask": {0: "batch"},
            "residual": {0: "batch"},
        },
        opset_version=17,
    )
    print(f"ONNX subtree model exported to {onnx_path}")


# ============================================================
# Entry Point
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train subtree-level GNN")
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--loss", type=str, default="huber", choices=["mse", "huber"])
    parser.add_argument("--huber_delta", type=float, default=1.0)
    parser.add_argument("--aa_embed_dim", type=int, default=16)
    parser.add_argument("--pos_embed_dim", type=int, default=8)
    parser.add_argument("--node_dim", type=int, default=32)
    parser.add_argument("--hidden_dim", type=int, default=64)
    parser.add_argument("--num_layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--val_frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--log_every", type=int, default=10)
    parser.add_argument("--min_confs", type=int, default=10,
                        help="Minimum confs per subtree to generate a training sample")
    parser.add_argument("--max_samples_per_k", type=int, default=50000,
                        help="Max subtree samples per level of fixed positions")
    args = parser.parse_args()

    train(args)
