"""
GNN Energy Surrogate — Ablation Training Script

Ablation variants:
  baseline   — current model as-is
  +E1        — add one-body emat energy as node feature
  +pairsum   — add per-node pairwise sum as node feature
  +emat_total— add E_emat_total to readout
  +bw_loss   — Boltzmann-weighted MSE loss
  +all_feat  — E1 + pairsum + emat_total + bw_loss
  larger     — double hidden_dim, 6 layers
  2x_data    — train on 400K data (needs separate export)
"""

import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from train import NUM_AA_TYPES, MAX_CHI, load_data, evaluate_all

kT = 0.5924  # kcal/mol at 298K


# ============================================================
# Extended Data Loading (adds one-body energies if available)
# ============================================================

def load_data_extended(data_dir):
    """Load standard data + one-body energies if available."""
    data = load_data(data_dir)

    # Try to load one-body energies
    onebody_path = os.path.join(data_dir, "onebody_energies.csv")
    if os.path.exists(onebody_path):
        ob_df = pd.read_csv(onebody_path)
        num_pos = data["num_pos"]
        max_rcs = data["max_rcs"]
        onebody_table = torch.zeros(num_pos, max_rcs, dtype=torch.float32)
        for _, row in ob_df.iterrows():
            p, r = int(row["pos"]), int(row["rc"])
            onebody_table[p, r] = row["E_onebody"]
        data["onebody_table"] = onebody_table
        print(f"  One-body energies loaded: {len(ob_df)} entries")
    else:
        data["onebody_table"] = None
        print(f"  One-body energies not found, skipping")

    return data


# ============================================================
# Dataset with E_CCD for Boltzmann weighting
# ============================================================

class AblationDataset(Dataset):
    def __init__(self, confs, targets, emat_energies, ccd_energies):
        self.confs = torch.tensor(confs, dtype=torch.long)
        self.targets = torch.tensor(targets, dtype=torch.float32)
        self.emat_energies = torch.tensor(emat_energies, dtype=torch.float32)
        self.ccd_energies = torch.tensor(ccd_energies, dtype=torch.float32)

    def __len__(self):
        return len(self.confs)

    def __getitem__(self, idx):
        return self.confs[idx], self.targets[idx], self.emat_energies[idx], self.ccd_energies[idx]


# ============================================================
# Ablation GNN Model
# ============================================================

class AblationGNN(nn.Module):
    def __init__(self, num_pos, max_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, onebody_table=None,
                 aa_embed_dim=16, pos_embed_dim=8, node_dim=64,
                 hidden_dim=128, num_layers=4, dropout=0.1,
                 use_e1=False, use_pairsum=False, use_emat_total=False):
        super().__init__()
        self.num_pos = num_pos
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.max_rcs = max_rcs
        self.use_e1 = use_e1
        self.use_pairsum = use_pairsum
        self.use_emat_total = use_emat_total

        # Shared AA type embedding
        self.aa_embedding = nn.Embedding(NUM_AA_TYPES, aa_embed_dim)
        self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        chi_input_dim = MAX_CHI * 2  # sin/cos

        # Node input dimension
        node_input_dim = aa_embed_dim + chi_input_dim + pos_embed_dim
        if use_e1:
            node_input_dim += 1  # E1 one-body
        if use_pairsum:
            node_input_dim += 1  # pairwise sum

        self.node_proj = nn.Sequential(
            nn.Linear(node_input_dim, node_dim),
            nn.SiLU(),
        )

        edge_feat_dim = 2  # [pair_energy, ca_distance]

        # Register buffers
        self.register_buffer("aa_table", aa_table)
        self.register_buffer("chi_table", chi_table)
        self.register_buffer("pair_table", pair_table)
        self.register_buffer("ca_dist_vec", ca_dist_vec)
        self.register_buffer("edge_src", edge_index[0])
        self.register_buffer("edge_dst", edge_index[1])

        if onebody_table is not None:
            self.register_buffer("onebody_table", onebody_table)
        else:
            self.register_buffer("onebody_table", torch.zeros(num_pos, max_rcs))

        # Message passing
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
        readout_extra = 1 if use_emat_total else 0
        self.node_readout = nn.Sequential(
            nn.Linear(node_dim + readout_extra, hidden_dim),
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

    def forward(self, confs, emat_energies=None):
        batch = confs.shape[0]
        num_pos = self.num_pos
        node_dim = self.node_dim

        # === Node features ===
        aa_idx = torch.gather(
            self.aa_table.unsqueeze(0).expand(batch, -1, -1),
            2, confs.unsqueeze(-1)
        ).squeeze(-1)
        aa_feat = self.aa_embedding(aa_idx)

        chi_raw = torch.gather(
            self.chi_table.unsqueeze(0).expand(batch, -1, -1, -1),
            2, confs.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 1, MAX_CHI)
        ).squeeze(2)
        chi_rad = chi_raw * (3.14159265358979 / 180.0)
        chi_feat = torch.cat([torch.sin(chi_rad), torch.cos(chi_rad)], dim=-1)

        pos_idx = torch.arange(num_pos, device=confs.device)
        pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)

        node_parts = [aa_feat, chi_feat, pos_feat]

        # E1 one-body feature
        if self.use_e1:
            e1 = torch.gather(
                self.onebody_table.unsqueeze(0).expand(batch, -1, -1),
                2, confs.unsqueeze(-1)
            ).squeeze(-1).unsqueeze(-1)  # (batch, num_pos, 1)
            node_parts.append(e1)

        # Pairwise sum per node
        if self.use_pairsum:
            src = self.edge_src
            dst = self.edge_dst
            rc_src = torch.gather(confs, 1, src.unsqueeze(0).expand(batch, -1))
            rc_dst = torch.gather(confs, 1, dst.unsqueeze(0).expand(batch, -1))
            flat_idx = rc_src * self.max_rcs + rc_dst
            pair_e = torch.gather(
                self.pair_table.unsqueeze(0).expand(batch, -1, -1),
                2, flat_idx.unsqueeze(-1)
            ).squeeze(-1)  # (batch, num_edges)
            # scatter_add to destination nodes
            pairsum = torch.zeros(batch, num_pos, device=confs.device)
            dst_exp = dst.unsqueeze(0).expand(batch, -1)
            pairsum.scatter_add_(1, dst_exp, pair_e)
            node_parts.append(pairsum.unsqueeze(-1))  # (batch, num_pos, 1)

        node_input = torch.cat(node_parts, dim=-1)
        h = self.node_proj(node_input)

        # === Edge features ===
        src = self.edge_src
        dst = self.edge_dst
        num_edges = src.shape[0]

        src_exp = src.unsqueeze(0).expand(batch, -1)
        dst_exp = dst.unsqueeze(0).expand(batch, -1)
        rc_src = torch.gather(confs, 1, src_exp)
        rc_dst = torch.gather(confs, 1, dst_exp)

        flat_idx = rc_src * self.max_rcs + rc_dst
        pair_energy = torch.gather(
            self.pair_table.unsqueeze(0).expand(batch, -1, -1),
            2, flat_idx.unsqueeze(-1)
        ).squeeze(-1)

        ca_dist = self.ca_dist_vec.unsqueeze(0).expand(batch, -1)
        edge_feats = torch.stack([pair_energy, ca_dist], dim=-1)

        # === Message passing ===
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
        if self.use_emat_total and emat_energies is not None:
            emat_feat = emat_energies.unsqueeze(-1).unsqueeze(-1).expand(-1, num_pos, 1)
            node_readout_input = torch.cat([h, emat_feat], dim=-1)
        else:
            node_readout_input = h

        node_energy = self.node_readout(node_readout_input).squeeze(-1)
        total_node = node_energy.sum(dim=1)

        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)
        edge_energy = self.edge_readout(edge_input).squeeze(-1)
        total_edge = edge_energy.sum(dim=1) / 2.0

        return total_node + total_edge


# ============================================================
# Per-pair MLP Model (LUTE-style, no message passing)
# ============================================================

class PairwiseMLP(nn.Module):
    """
    Directly predict per-pair and per-node energy corrections.
    residual = Σ_i f_node(rc_i) + Σ_{i<j} f_pair(rc_i, rc_j, E2_ij)

    No message passing — each pair is independent.
    Captures the same physics as LUTE but with neural networks.
    """

    def __init__(self, num_pos, max_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, onebody_table=None,
                 aa_embed_dim=16, pos_embed_dim=8, node_dim=64,
                 hidden_dim=128, num_layers=4, dropout=0.1, **kwargs):
        super().__init__()
        self.num_pos = num_pos
        self.max_rcs = max_rcs

        self.aa_embedding = nn.Embedding(NUM_AA_TYPES, aa_embed_dim)
        self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        chi_input_dim = MAX_CHI * 2
        rc_feat_dim = aa_embed_dim + chi_input_dim + pos_embed_dim  # per-node RC feature

        # Register buffers
        self.register_buffer("aa_table", aa_table)
        self.register_buffer("chi_table", chi_table)
        self.register_buffer("pair_table", pair_table)
        self.register_buffer("ca_dist_vec", ca_dist_vec)
        self.register_buffer("edge_src", edge_index[0])
        self.register_buffer("edge_dst", edge_index[1])
        if onebody_table is not None:
            self.register_buffer("onebody_table", onebody_table)
        else:
            self.register_buffer("onebody_table", torch.zeros(num_pos, max_rcs))

        # Per-node correction: f(rc_feat) → scalar
        self.node_mlp = nn.Sequential(
            nn.Linear(rc_feat_dim + 1, hidden_dim),  # +1 for E1 one-body
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

        # Per-pair correction: f(rc_feat_i, rc_feat_j, E2_ij, ca_dist) → scalar
        self.pair_mlp = nn.Sequential(
            nn.Linear(2 * rc_feat_dim + 2, hidden_dim),  # +2 for E2, ca_dist
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def _rc_features(self, confs):
        """Get per-node RC features: aa_embed + chi_sincos + pos_embed."""
        batch = confs.shape[0]

        aa_idx = torch.gather(
            self.aa_table.unsqueeze(0).expand(batch, -1, -1),
            2, confs.unsqueeze(-1)
        ).squeeze(-1)
        aa_feat = self.aa_embedding(aa_idx)

        chi_raw = torch.gather(
            self.chi_table.unsqueeze(0).expand(batch, -1, -1, -1),
            2, confs.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 1, MAX_CHI)
        ).squeeze(2)
        chi_rad = chi_raw * (3.14159265358979 / 180.0)
        chi_feat = torch.cat([torch.sin(chi_rad), torch.cos(chi_rad)], dim=-1)

        pos_idx = torch.arange(self.num_pos, device=confs.device)
        pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)

        return torch.cat([aa_feat, chi_feat, pos_feat], dim=-1)  # (batch, num_pos, rc_feat_dim)

    def forward(self, confs, emat_energies=None):
        batch = confs.shape[0]
        rc_feat = self._rc_features(confs)  # (batch, num_pos, D)

        # Per-node correction
        e1 = torch.gather(
            self.onebody_table.unsqueeze(0).expand(batch, -1, -1),
            2, confs.unsqueeze(-1)
        ).squeeze(-1).unsqueeze(-1)  # (batch, num_pos, 1)
        node_input = torch.cat([rc_feat, e1], dim=-1)
        node_corr = self.node_mlp(node_input).squeeze(-1).sum(dim=1)  # (batch,)

        # Per-pair correction (use undirected edges only: src < dst)
        src, dst = self.edge_src, self.edge_dst
        mask = src < dst  # undirected: only i<j
        src_u, dst_u = src[mask], dst[mask]

        src_exp = src_u.unsqueeze(0).expand(batch, -1)
        dst_exp = dst_u.unsqueeze(0).expand(batch, -1)

        rc_src = torch.gather(rc_feat, 1,
                              src_exp.unsqueeze(-1).expand(-1, -1, rc_feat.shape[-1]))
        rc_dst = torch.gather(rc_feat, 1,
                              dst_exp.unsqueeze(-1).expand(-1, -1, rc_feat.shape[-1]))

        # E2 pairwise energy for this conformation
        rc_s = torch.gather(confs, 1, src_exp)
        rc_d = torch.gather(confs, 1, dst_exp)

        # Find edge indices for src_u, dst_u
        # pair_table is indexed by directed edge index, so find the right ones
        edge_mask_indices = mask.nonzero(as_tuple=True)[0]
        flat_idx = rc_s * self.max_rcs + rc_d
        pair_e = torch.gather(
            self.pair_table[edge_mask_indices].unsqueeze(0).expand(batch, -1, -1),
            2, flat_idx.unsqueeze(-1)
        ).squeeze(-1).unsqueeze(-1)  # (batch, n_undirected, 1)

        ca_d = self.ca_dist_vec[edge_mask_indices].unsqueeze(0).expand(batch, -1).unsqueeze(-1)

        pair_input = torch.cat([rc_src, rc_dst, pair_e, ca_d], dim=-1)
        pair_corr = self.pair_mlp(pair_input).squeeze(-1).sum(dim=1)  # (batch,)

        return node_corr + pair_corr


# ============================================================
# Ablation configs
# ============================================================

ABLATION_CONFIGS = {
    # === Single-factor ablations ===
    "baseline": dict(
        use_e1=False, use_pairsum=False, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+E1": dict(
        use_e1=True, use_pairsum=False, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+pairsum": dict(
        use_e1=False, use_pairsum=True, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+emat_total": dict(
        use_e1=False, use_pairsum=False, use_emat_total=True, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+bw_loss": dict(
        use_e1=False, use_pairsum=False, use_emat_total=False, bw_loss=True,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+stratified": dict(
        use_e1=False, use_pairsum=False, use_emat_total=False, bw_loss=False,
        stratified=True, node_dim=64, hidden_dim=128, num_layers=4, model_type="gnn",
    ),
    "+target_norm": dict(
        use_e1=False, use_pairsum=False, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4,
        target_norm=True, model_type="gnn",
    ),
    "+huber": dict(
        use_e1=False, use_pairsum=False, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4,
        huber=True, model_type="gnn",
    ),
    # === Combinations ===
    "huber+pairsum": dict(
        use_e1=False, use_pairsum=True, use_emat_total=False, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4,
        huber=True, model_type="gnn",
    ),
    "huber+ps+et+tn": dict(
        use_e1=False, use_pairsum=True, use_emat_total=True, bw_loss=False,
        stratified=False, node_dim=64, hidden_dim=128, num_layers=4,
        huber=True, target_norm=True, model_type="gnn",
    ),
}


# ============================================================
# Training loop
# ============================================================

def train_one(config_name, config, data, args, out_dir):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n{'='*70}")
    print(f" ABLATION: {config_name}")
    print(f"{'='*70}")
    print(f"  Config: {config}")

    confs = data["confs"]
    residuals = data["residuals"]
    emat_energies = data["emat_energies"]
    ccd_energies = data["ccd_energies"]

    # Seed EVERYTHING for reproducibility
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)
    np.random.seed(args.seed)

    # Train/val split
    n = len(confs)
    perm = np.random.RandomState(args.seed).permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    train_ds = AblationDataset(confs[train_idx], residuals[train_idx],
                               emat_energies[train_idx], ccd_energies[train_idx])
    val_ds = AblationDataset(confs[val_idx], residuals[val_idx],
                             emat_energies[val_idx], ccd_energies[val_idx])

    use_stratified = config.get("stratified", False)
    if use_stratified:
        # Weight = 1/(seq_count) so each sequence is sampled equally
        from collections import Counter
        rc_feat_df = pd.read_csv(os.path.join(data["data_dir"], "rc_features.csv"))
        rc_to_aa = {}
        for pos in range(data["num_pos"]):
            pos_rcs = rc_feat_df[rc_feat_df.pos == pos]
            rc_to_aa[pos] = dict(zip(pos_rcs.rc, pos_rcs.aa_type_idx))

        train_confs = confs[train_idx]
        train_seqs = [tuple(rc_to_aa[p][train_confs[i, p]]
                            for p in range(data["num_pos"]))
                      for i in range(len(train_confs))]
        seq_counts = Counter(train_seqs)
        sample_weights = torch.tensor([1.0 / seq_counts[s] for s in train_seqs],
                                      dtype=torch.float64)
        sampler = torch.utils.data.WeightedRandomSampler(
            sample_weights, num_samples=len(train_ds), replacement=True)
        train_loader = DataLoader(train_ds, batch_size=args.batch_size,
                                  sampler=sampler, num_workers=2, pin_memory=True)
        n_seqs = len(seq_counts)
        print(f"  Stratified sampling: {n_seqs} sequences, "
              f"max_count={max(seq_counts.values())}, min_count={min(seq_counts.values())}")
    else:
        train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                                  num_workers=2, pin_memory=True)

    val_loader = DataLoader(val_ds, batch_size=args.batch_size,
                            num_workers=2, pin_memory=True)

    print(f"  Train: {len(train_ds)}, Val: {len(val_ds)}")

    # --- Target normalization ---
    use_target_norm = config.get("target_norm", False)
    if use_target_norm:
        t_mean = float(residuals[train_idx].mean())
        t_std = float(residuals[train_idx].std())
        if t_std < 1e-8:
            t_std = 1.0
        print(f"  Target normalization: mean={t_mean:.4f}, std={t_std:.4f}")
    else:
        t_mean, t_std = 0.0, 1.0

    # --- Model selection ---
    model_type = config.get("model_type", "gnn")
    model_kwargs = dict(
        num_pos=data["num_pos"], max_rcs=data["max_rcs"],
        aa_table=data["aa_table"], chi_table=data["chi_table"],
        pair_table=data["pair_table"], edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"], onebody_table=data["onebody_table"],
        node_dim=config["node_dim"], hidden_dim=config["hidden_dim"],
        num_layers=config["num_layers"], dropout=args.dropout,
    )
    if model_type == "pairwise_mlp":
        model = PairwiseMLP(**model_kwargs).to(device)
    else:
        model = AblationGNN(
            **model_kwargs,
            use_e1=config["use_e1"], use_pairsum=config["use_pairsum"],
            use_emat_total=config["use_emat_total"],
        ).to(device)

    num_params = sum(p.numel() for p in model.parameters())
    print(f"  Parameters: {num_params:,}")
    print(f"  Model type: {model_type}")

    # --- Loss selection ---
    use_bw = config.get("bw_loss", False)
    use_huber = config.get("huber", False)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_val_rmse = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        train_loss = 0.0
        for confs_b, targets_b, emat_b, ccd_b in train_loader:
            confs_b = confs_b.to(device)
            emat_b = emat_b.to(device)
            ccd_b = ccd_b.to(device)

            # Normalize target
            targets_norm = (targets_b - t_mean) / t_std
            targets_norm = targets_norm.to(device)

            pred_norm = model(confs_b, emat_b)

            if use_bw:
                with torch.no_grad():
                    # Log-space to avoid overflow: log_w = -E/kT
                    log_w = -ccd_b / kT
                    log_w = log_w - log_w.max()  # shift for numerical stability
                    w = torch.exp(log_w)
                    w = w / w.mean()  # normalize so mean weight = 1
                    w = w.clamp(max=50.0)
                loss = (w * (pred_norm - targets_norm) ** 2).mean()
            elif use_huber:
                loss = nn.functional.huber_loss(pred_norm, targets_norm, delta=1.0)
            else:
                loss = nn.functional.mse_loss(pred_norm, targets_norm)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            train_loss += loss.item() * len(confs_b)

        train_rmse = (train_loss / len(train_ds)) ** 0.5
        scheduler.step()

        # Validation (always unweighted MSE in ORIGINAL scale for fair comparison)
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for confs_b, targets_b, emat_b, _ in val_loader:
                confs_b = confs_b.to(device)
                targets_b = targets_b.to(device)
                emat_b = emat_b.to(device)
                pred_norm = model(confs_b, emat_b)
                # Denormalize for fair comparison
                pred_orig = pred_norm * t_std + t_mean
                val_loss += nn.functional.mse_loss(pred_orig, targets_b).item() * len(confs_b)

        val_rmse = (val_loss / len(val_ds)) ** 0.5
        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if epoch % args.log_every == 0 or epoch == 1:
            print(f"  Ep {epoch:4d}  train={train_rmse:.4f}  val={val_rmse:.4f}  best={best_val_rmse:.4f}")

    # Restore best and evaluate
    if best_state is None:
        print(f"  WARNING: training failed (no valid checkpoint). Returning NaN.")
        return {
            "config_name": config_name, "num_params": num_params,
            "val_rmse": float("nan"), "mae": float("nan"), "max_err": float("nan"),
        }

    model.load_state_dict(best_state)
    model.eval()

    all_preds, all_targets, all_emat = [], [], []
    with torch.no_grad():
        for confs_b, targets_b, emat_b, _ in val_loader:
            confs_b = confs_b.to(device)
            emat_b = emat_b.to(device)
            pred_norm = model(confs_b, emat_b)
            pred_orig = pred_norm * t_std + t_mean
            all_preds.append(pred_orig.cpu().numpy())
            all_targets.append(targets_b.numpy())
            all_emat.append(emat_b.numpy())

    preds = np.concatenate(all_preds)
    targets = np.concatenate(all_targets)
    emats = np.concatenate(all_emat)

    evaluate_all(preds, targets, emats)

    # Save checkpoint
    model_dir = os.path.join(out_dir, config_name)
    os.makedirs(model_dir, exist_ok=True)
    torch.save({
        "model_state": best_state,
        "config": config,
        "config_name": config_name,
        "num_params": num_params,
        "val_rmse": best_val_rmse,
        "target_norm": {"mean": t_mean, "std": t_std},
    }, os.path.join(model_dir, "checkpoint.pt"))

    return {
        "config_name": config_name,
        "num_params": num_params,
        "val_rmse": best_val_rmse,
        "mae": float(np.mean(np.abs(preds - targets))),
        "max_err": float(np.max(np.abs(preds - targets))),
    }


# ============================================================
# Per-sequence pfunc evaluation
# ============================================================

def eval_per_seq(config_name, model, data, device):
    """Quick per-sequence pfunc evaluation on validation set."""
    confs = data["confs"]
    emat_energies = data["emat_energies"]
    ccd_energies = data["ccd_energies"]

    # Predict all
    conf_tensor = torch.tensor(confs, dtype=torch.long, device=device)
    emat_tensor = torch.tensor(emat_energies, dtype=torch.float32, device=device)

    E_GNN = np.zeros(len(confs))
    bs = 50000
    model.eval()
    with torch.no_grad():
        for s in range(0, len(confs), bs):
            e = min(s + bs, len(confs))
            pred = model(conf_tensor[s:e], emat_tensor[s:e]).cpu().numpy()
            E_GNN[s:e] = emat_energies[s:e] + pred

    # Map rc -> aa
    rc_feat_df = pd.read_csv(os.path.join(data["data_dir"], "rc_features.csv"))
    rc_to_aa = {}
    for pos in range(data["num_pos"]):
        pos_rcs = rc_feat_df[rc_feat_df.pos == pos]
        rc_to_aa[pos] = dict(zip(pos_rcs.rc, pos_rcs.aa_type_idx))

    from collections import defaultdict
    seq_groups = defaultdict(list)
    for i in range(len(confs)):
        sk = tuple(rc_to_aa[p][confs[i, p]] for p in range(data["num_pos"]))
        seq_groups[sk].append(i)

    log_ratios = []
    rel_errs = []
    for seq, idxs in seq_groups.items():
        idxs = np.array(idxs)
        pf_true = np.exp(-ccd_energies[idxs] / kT).sum()
        pf_gnn = np.exp(-E_GNN[idxs] / kT).sum()
        if pf_true > 0 and pf_gnn > 0:
            log_ratios.append(np.log10(pf_gnn / pf_true))
            rel_errs.append(abs(pf_gnn - pf_true) / pf_true)

    log_ratios = np.array(log_ratios)
    rel_errs = np.array(rel_errs)
    return {
        "n_seqs": len(log_ratios),
        "dlog10_mean": float(np.mean(np.abs(log_ratios))),
        "dlog10_median": float(np.median(np.abs(log_ratios))),
        "dlog10_max": float(np.max(np.abs(log_ratios))),
        "dlog10_gt01": int((np.abs(log_ratios) > 0.01).sum()),
        "rel_err_median": float(np.median(rel_errs)),
        "rel_err_p95": float(np.percentile(rel_errs, 95)),
        "rel_err_max": float(np.max(rel_errs)),
    }


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--data_dir_2x", default=None, help="400K data dir for 2x_data ablation")
    parser.add_argument("--out_dir", default=None)
    parser.add_argument("--epochs", type=int, default=500)
    parser.add_argument("--batch_size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--val_frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--log_every", type=int, default=50)
    parser.add_argument("--only", default=None, help="Run only this variant")
    parser.add_argument("--label", default="protein", help="protein or complex")
    args = parser.parse_args()

    data_path = os.path.join(args.data_dir, args.label)
    if args.out_dir is None:
        args.out_dir = os.path.join(args.data_dir, "ablation", args.label)
    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading data from {data_path}")
    data = load_data_extended(data_path)
    data["data_dir"] = data_path

    # Determine which configs to run
    configs = dict(ABLATION_CONFIGS)
    if args.data_dir_2x:
        configs["2x_data"] = dict(ABLATION_CONFIGS["baseline"])  # same arch, more data

    if args.only:
        configs = {args.only: configs[args.only]}

    results = {}
    for name, config in configs.items():
        r = train_one(name, config, data, args, args.out_dir)
        results[name] = r

    # Summary table
    print(f"\n{'='*80}")
    print(f" ABLATION SUMMARY — {args.label}")
    print(f"{'='*80}")
    print(f"{'Variant':<16} {'Params':>10} {'ValRMSE':>10} {'MAE':>10} {'MaxErr':>10}")
    print("-" * 80)
    for name, r in sorted(results.items(), key=lambda x: x[1]["val_rmse"]):
        print(f"{name:<16} {r['num_params']:>10,} {r['val_rmse']:>10.6f} "
              f"{r['mae']:>10.6f} {r['max_err']:>10.4f}")

    # Save results
    results_path = os.path.join(args.out_dir, "ablation_results.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {results_path}")


if __name__ == "__main__":
    main()
