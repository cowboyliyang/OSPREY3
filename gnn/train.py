"""
GNN Energy Surrogate Training Script v2

Improved architecture with:
  - Shared AA type embeddings (instead of position-specific RC embeddings)
  - Chi angle encoding (sin/cos of rotamer dihedral angles)
  - Position embeddings
  - Edge features from emat pairwise energies + Cα distances

Reads data exported by GNNDataExporter:
  confs.csv, graph.csv, meta.csv, rc_features.csv, pairwise_energies.csv, ca_distances.csv

Usage:
    python train.py --data_dir gnn_data/2RL0_flex8 --epochs 200
"""

import argparse
import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader


NUM_AA_TYPES = 20
UNKNOWN_AA_IDX = NUM_AA_TYPES  # index 20 for non-standard residues (aa_type_idx == -1)
MAX_CHI = 4

# Standard amino acid 3-letter codes → index (must match Java side)
AA_TYPES = {
    "ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4,
    "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9,
    "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14,
    "SER": 15, "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19,
}


# ============================================================
# Data Loading
# ============================================================

class ConfDataset(Dataset):
    def __init__(self, confs, targets, emat_energies):
        self.confs = torch.tensor(confs, dtype=torch.long)
        self.targets = torch.tensor(targets, dtype=torch.float32)
        self.emat_energies = torch.tensor(emat_energies, dtype=torch.float32)

    def __len__(self):
        return len(self.confs)

    def __getitem__(self, idx):
        return self.confs[idx], self.targets[idx], self.emat_energies[idx]


def load_data(data_dir):
    """Load all exported CSV files."""
    import time as _time
    t0 = _time.time()

    # Use numpy for the big file, pandas for small metadata files
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
    ccd_energies = confs_df["E_CCD"].values.astype(np.float64)
    has_rigid = "E_rigid" in confs_df.columns
    rigid_energies = confs_df["E_rigid"].values.astype(np.float64) if has_rigid else None

    # Free the big dataframe early
    del confs_df

    # Filter: emat pre-filter (consistent with sampling & inference)
    n_raw = len(confs)
    mask = np.ones(n_raw, dtype=bool)

    # 1. E_emat (minimizing) > -20 → skip
    emat_cap = -20.0
    mask_emat = emat_energies <= emat_cap
    mask &= mask_emat

    # 2. E_rigid > 0 → skip (if column exists)
    if has_rigid:
        rigid_cap = 0.0
        mask_rigid = rigid_energies <= rigid_cap
        mask &= mask_rigid

    # 3. |residual| >= 100 → skip (training quality guard)
    residual_cap = 100.0
    mask_res = np.abs(residuals) < residual_cap
    mask &= mask_res

    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    ccd_energies = ccd_energies[mask]
    if has_rigid:
        rigid_energies = rigid_energies[mask]

    n_kept = len(confs)
    print(f"Filtered: {n_raw} -> {n_kept} conformations (removed {n_raw - n_kept})")
    print(f"  E_emat <= {emat_cap}: removed {int((~mask_emat).sum())}")
    if has_rigid:
        print(f"  E_rigid <= {rigid_cap}: removed {int((~mask_rigid).sum())}")
    print(f"  |residual| < {residual_cap}: removed {int((~mask_res).sum())}")

    edge_index = torch.tensor(graph_df[["src", "dst"]].values.T, dtype=torch.long)
    num_rcs = meta_df["num_rcs"].values.tolist()
    max_rcs = max(num_rcs)

    # RC features: build padded lookup tables (vectorized)
    aa_table = torch.zeros(num_pos, max_rcs, dtype=torch.long)
    chi_table = torch.zeros(num_pos, max_rcs, MAX_CHI, dtype=torch.float32)

    rc_feat_df = rc_feat_df.sort_values(["pos", "rc"])
    for p in range(num_pos):
        pos_mask = rc_feat_df["pos"].values == p
        pos_data = rc_feat_df[pos_mask]
        n = len(pos_data)
        aa_vals = pos_data["aa_type_idx"].values
        # Remap non-standard residues (-1) to UNKNOWN_AA_IDX so embedding lookup is valid
        aa_vals = np.where(aa_vals < 0, UNKNOWN_AA_IDX, aa_vals)
        aa_table[p, :n] = torch.tensor(aa_vals, dtype=torch.long)
        chis = pos_data[["chi1", "chi2", "chi3", "chi4"]].fillna(0.0).values
        chi_table[p, :n] = torch.tensor(chis, dtype=torch.float32)

    # Pairwise energy lookup: vectorized build
    src_edges, dst_edges = edge_index
    num_edges = src_edges.shape[0]
    pair_table = torch.zeros(num_edges, max_rcs * max_rcs, dtype=torch.float32)

    # Build edge-index map for fast lookup
    edge_map = {}
    for e in range(num_edges):
        edge_map[(src_edges[e].item(), dst_edges[e].item())] = e

    # Vectorized: process entire pair_df at once
    p1 = pair_df["pos1"].values
    r1 = pair_df["rc1"].values
    p2 = pair_df["pos2"].values
    r2 = pair_df["rc2"].values
    evals = pair_df["E_pair_min"].values

    for i in range(len(evals)):
        # Try both edge directions
        e_idx = edge_map.get((p1[i], p2[i]))
        if e_idx is not None:
            pair_table[e_idx, r1[i] * max_rcs + r2[i]] = evals[i]
        e_idx = edge_map.get((p2[i], p1[i]))
        if e_idx is not None:
            pair_table[e_idx, r2[i] * max_rcs + r1[i]] = evals[i]

    # Cα distances per edge (vectorized)
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
    print(f"Residual stats: mean={residuals.mean():.4f}, std={residuals.std():.4f}, "
          f"min={residuals.min():.4f}, max={residuals.max():.4f}")

    return {
        "confs": confs,
        "residuals": residuals,
        "emat_energies": emat_energies,
        "ccd_energies": ccd_energies,
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

class InteractionGNNv2(nn.Module):
    """
    Pairwise-energy-aware message-passing GNN.

    Improvements over v1:
      - Shared AA type embedding across all positions
      - Chi angle encoding (sin/cos)
      - Position embedding for local environment
      - Edge features from emat pairwise energies + Cα distance
    """

    def __init__(self, num_pos, max_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec,
                 aa_embed_dim=16, pos_embed_dim=8, node_dim=32,
                 hidden_dim=64, num_layers=3, dropout=0.1):
        super().__init__()
        self.num_pos = num_pos
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.max_rcs = max_rcs

        edge_feat_dim = 2  # [pair_energy, ca_distance]

        # Shared AA type embedding (+1 slot for unknown / non-standard residues)
        self.aa_embedding = nn.Embedding(NUM_AA_TYPES + 1, aa_embed_dim)

        # Position embedding
        self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        # Chi angle: sin/cos for up to 4 chi angles = 8 dims
        chi_input_dim = MAX_CHI * 2

        # Project combined node features to node_dim
        self.node_proj = nn.Sequential(
            nn.Linear(aa_embed_dim + chi_input_dim + pos_embed_dim, node_dim),
            nn.SiLU(),
        )

        # Register lookup tables as buffers (baked into ONNX)
        self.register_buffer("aa_table", aa_table)          # (num_pos, max_rcs)
        self.register_buffer("chi_table", chi_table)        # (num_pos, max_rcs, 4)
        self.register_buffer("pair_table", pair_table)      # (num_edges, max_rcs^2)
        self.register_buffer("ca_dist_vec", ca_dist_vec)    # (num_edges,)
        self.register_buffer("edge_src", edge_index[0])     # (num_edges,)
        self.register_buffer("edge_dst", edge_index[1])     # (num_edges,)

        # Message passing layers with edge features
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

    def forward(self, confs):
        """
        Args:
            confs: (batch, num_pos) RC indices per position
        Returns:
            residual: (batch,) predicted E_CCD - E_emat
        """
        batch = confs.shape[0]
        num_pos = self.num_pos
        node_dim = self.node_dim

        # === Node features ===
        # AA type lookup: aa_table[p, confs[b, p]] for all b, p
        aa_idx = torch.gather(
            self.aa_table.unsqueeze(0).expand(batch, -1, -1),
            2, confs.unsqueeze(-1)
        ).squeeze(-1)  # (batch, num_pos)
        aa_feat = self.aa_embedding(aa_idx)  # (batch, num_pos, aa_embed_dim)

        # Chi angle lookup + sin/cos encoding
        chi_raw = torch.gather(
            self.chi_table.unsqueeze(0).expand(batch, -1, -1, -1),
            2, confs.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 1, MAX_CHI)
        ).squeeze(2)  # (batch, num_pos, 4)
        chi_rad = chi_raw * (3.14159265358979 / 180.0)
        chi_feat = torch.cat([torch.sin(chi_rad), torch.cos(chi_rad)], dim=-1)
        # (batch, num_pos, 8)

        # Position embedding (broadcast)
        pos_idx = torch.arange(num_pos, device=confs.device)
        pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)
        # (batch, num_pos, pos_embed_dim)

        # Combine and project
        node_input = torch.cat([aa_feat, chi_feat, pos_feat], dim=-1)
        h = self.node_proj(node_input)  # (batch, num_pos, node_dim)

        # === Edge features ===
        src = self.edge_src  # (num_edges,)
        dst = self.edge_dst
        num_edges = src.shape[0]

        # Gather RC indices at edge endpoints
        src_exp = src.unsqueeze(0).expand(batch, -1)  # (batch, num_edges)
        dst_exp = dst.unsqueeze(0).expand(batch, -1)
        rc_src = torch.gather(confs, 1, src_exp)  # (batch, num_edges)
        rc_dst = torch.gather(confs, 1, dst_exp)

        # Pairwise energy lookup — advanced indexing avoids the [batch, num_edges,
        # max_rcs**2] intermediate that triggers cudaErrorInvalidConfiguration in
        # the exported ONNX Expand kernel at batch>~700 (see train_subtree.py for
        # the same fix on the subtree model).
        flat_idx = rc_src * self.max_rcs + rc_dst  # (batch, num_edges)
        edge_arange = torch.arange(num_edges, device=confs.device)
        pair_energy = self.pair_table[edge_arange, flat_idx]  # (batch, num_edges)

        # Cα distance (static, same for all samples)
        ca_dist = self.ca_dist_vec.unsqueeze(0).expand(batch, -1)

        edge_feats = torch.stack([pair_energy, ca_dist], dim=-1)
        # (batch, num_edges, 2)

        # === Message passing ===
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)  # (batch, num_edges, node_dim)
            h_dst = torch.gather(h, 1, dst_idx)

            msg_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)
            messages = self.message_mlps[layer](msg_input)

            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)

            update_input = torch.cat([h, agg], dim=-1)
            h = h + self.update_mlps[layer](update_input)
            h = self.layer_norms[layer](h)

        # === Readout ===
        node_energy = self.node_readout(h).squeeze(-1)  # (batch, num_pos)
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
    confs = data["confs"]
    residuals = data["residuals"]
    emat_energies = data["emat_energies"]

    # Train/val split
    n = len(confs)
    perm = np.random.RandomState(args.seed).permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    train_ds = ConfDataset(confs[train_idx], residuals[train_idx], emat_energies[train_idx])
    val_ds = ConfDataset(confs[val_idx], residuals[val_idx], emat_energies[val_idx])
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)

    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")

    model = InteractionGNNv2(
        num_pos=data["num_pos"],
        max_rcs=data["max_rcs"],
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
    start_epoch = 1

    # Resume from checkpoint
    if args.resume:
        print(f"Resuming from {args.resume}")
        ckpt = torch.load(args.resume, map_location=device)
        base = model.module if isinstance(model, nn.DataParallel) else model
        base.load_state_dict(ckpt["model_state"])
        if "optimizer_state" in ckpt:
            optimizer.load_state_dict(ckpt["optimizer_state"])
        if "scheduler_state" in ckpt:
            scheduler.load_state_dict(ckpt["scheduler_state"])
        if "epoch" in ckpt:
            start_epoch = ckpt["epoch"] + 1
        if "val_rmse" in ckpt:
            best_val_rmse = ckpt["val_rmse"]
            best_state = ckpt["model_state"]
        print(f"Resumed at epoch {start_epoch}, best_val_rmse={best_val_rmse:.4f}")

    for epoch in range(start_epoch, args.epochs + 1):
        model.train()
        train_mse = 0.0
        for confs_batch, targets_batch, _ in train_loader:
            confs_batch = confs_batch.to(device)
            targets_batch = targets_batch.to(device)

            pred = model(confs_batch)
            if args.loss == "huber":
                loss = nn.functional.huber_loss(pred, targets_batch, delta=args.huber_delta)
            else:
                loss = nn.functional.mse_loss(pred, targets_batch)
            train_mse += nn.functional.mse_loss(pred, targets_batch).item() * len(confs_batch)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

        train_rmse = (train_mse / len(train_ds)) ** 0.5
        scheduler.step()

        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for confs_batch, targets_batch, _ in val_loader:
                confs_batch = confs_batch.to(device)
                targets_batch = targets_batch.to(device)
                pred = model(confs_batch)
                val_loss += nn.functional.mse_loss(pred, targets_batch).item() * len(confs_batch)
        val_rmse = (val_loss / len(val_ds)) ** 0.5

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            # Strip 'module.' prefix from DataParallel for clean checkpoint
            raw_sd = model.module.state_dict() if isinstance(model, nn.DataParallel) else model.state_dict()
            best_state = {k: v.cpu().clone() for k, v in raw_sd.items()}

        if epoch % args.log_every == 0 or epoch == start_epoch:
            print(f"Epoch {epoch:4d}  train_rmse={train_rmse:.4f}  "
                  f"val_rmse={val_rmse:.4f}  best={best_val_rmse:.4f}  "
                  f"lr={scheduler.get_last_lr()[0]:.2e}", flush=True)

        # Periodic checkpoint every 50 epochs
        if epoch % 50 == 0:
            ckpt_dir = os.path.join(args.data_dir, "model")
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
            }, os.path.join(ckpt_dir, "gnn_checkpoint.pt"))

    # Restore best model (unwrap DataParallel for eval/export)
    base_model = model.module if isinstance(model, nn.DataParallel) else model
    base_model.load_state_dict(best_state)
    base_model.eval()
    model = base_model

    # Final evaluation
    all_preds, all_targets, all_emat = [], [], []
    with torch.no_grad():
        for confs_batch, targets_batch, emat_batch in val_loader:
            confs_batch = confs_batch.to(device)
            pred = model(confs_batch)
            all_preds.append(pred.cpu().numpy())
            all_targets.append(targets_batch.numpy())
            all_emat.append(emat_batch.numpy())

    all_preds = np.concatenate(all_preds)
    all_targets = np.concatenate(all_targets)
    all_emat = np.concatenate(all_emat)

    evaluate_all(all_preds, all_targets, all_emat)

    # ---- Conformal Prediction calibration on val set ----
    # Compute nonconformity scores = |E_GNN - E_CCD| on val set
    val_errors = np.abs(all_preds - all_targets)  # residual-space errors
    n_cal = len(val_errors)
    cp_stats = {}
    for alpha in [0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001]:
        # Finite-sample corrected quantile: ceil((1-α)(n+1)) / n
        level = min(np.ceil((1 - alpha) * (n_cal + 1)) / n_cal, 1.0)
        q = float(np.quantile(val_errors, level))
        cp_stats[f"q_alpha_{alpha}"] = q
        print(f"  CP calibration: α={alpha:.3f}  q={q:.6f} kcal/mol  "
              f"(n_cal={n_cal}, coverage≥{1-alpha:.3f})")

    # Also store raw percentiles for reference
    for p in [50, 90, 95, 99, 99.5, 99.9]:
        cp_stats[f"P{p}"] = float(np.percentile(val_errors, p))
    cp_stats["max"] = float(np.max(val_errors))
    cp_stats["n_cal"] = n_cal

    print(f"\n  Val error percentiles:")
    print(f"    P50={cp_stats['P50']:.6f}  P90={cp_stats['P90']:.6f}  "
          f"P95={cp_stats['P95']:.6f}  P99={cp_stats['P99']:.6f}  "
          f"Max={cp_stats['max']:.6f}")

    # Save
    out_dir = os.path.join(args.data_dir, "model")
    os.makedirs(out_dir, exist_ok=True)

    torch.save({
        "model_state": best_state,
        "num_rcs": data["num_rcs"],
        "num_pos": data["num_pos"],
        "max_rcs": data["max_rcs"],
        "args": vars(args),
        "val_rmse": best_val_rmse,
        "cp_stats": cp_stats,
    }, os.path.join(out_dir, "gnn_checkpoint.pt"))

    try:
        export_onnx(model, data["num_pos"], out_dir, device)
    except Exception as e:
        print(f"ONNX export failed: {e}")
        import traceback
        traceback.print_exc()

    print(f"\nModel saved to {out_dir}")


def evaluate_all(residual_preds, residual_targets, emat_energies):
    n = len(residual_preds)
    gnn_total = emat_energies + residual_preds
    ccd_total = emat_energies + residual_targets
    emat_total = emat_energies

    gnn_delta = gnn_total - ccd_total
    emat_delta = emat_total - ccd_total

    print("\n" + "=" * 70)
    print(f"GNN vs CCD  (validation set, N={n})")
    print("=" * 70)
    print(f"  MAE:      {np.mean(np.abs(gnn_delta)):.6f} kcal/mol")
    print(f"  RMSE:     {np.sqrt(np.mean(gnn_delta**2)):.6f} kcal/mol")
    print(f"  Bias:     {np.mean(gnn_delta):.6f} kcal/mol")
    print(f"  Max|d|:   {np.max(np.abs(gnn_delta)):.6f} kcal/mol")

    try:
        from scipy.stats import spearmanr
        rho, _ = spearmanr(gnn_total, ccd_total)
        print(f"  Spearman: {rho:.4f}")
    except Exception:
        pass

    print(f"\nemat vs CCD  (baseline)")
    print(f"  MAE:      {np.mean(np.abs(emat_delta)):.6f} kcal/mol")
    print(f"  RMSE:     {np.sqrt(np.mean(emat_delta**2)):.6f} kcal/mol")

    gnn_mae = np.mean(np.abs(gnn_delta))
    emat_mae = np.mean(np.abs(emat_delta))
    if emat_mae > 0:
        print(f"\nMAE reduction: {(1 - gnn_mae / emat_mae) * 100:.1f}%")


def export_onnx(model, num_pos, out_dir, device):
    model.eval()
    dummy = torch.zeros(1, num_pos, dtype=torch.long, device=device)

    onnx_path = os.path.join(out_dir, "gnn_model.onnx")
    torch.onnx.export(
        model,
        dummy,
        onnx_path,
        input_names=["confs"],
        output_names=["residual"],
        dynamic_axes={"confs": {0: "batch"}, "residual": {0: "batch"}},
        opset_version=17,
    )
    print(f"ONNX model exported to {onnx_path}")


# ============================================================
# Entry Point
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GNN energy surrogate v2")
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
    parser.add_argument("--resume", type=str, default=None,
                        help="Path to checkpoint to resume from (e.g. gnn_data/.../model/gnn_checkpoint.pt)")
    args = parser.parse_args()

    train(args)
