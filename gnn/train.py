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
    confs_df = pd.read_csv(os.path.join(data_dir, "confs.csv"))
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

    # Filter extreme energies
    n_raw = len(confs)
    energy_cap = np.percentile(np.abs(ccd_energies), 95)
    energy_cap = max(energy_cap, 500.0)
    mask = np.abs(ccd_energies) < energy_cap
    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    ccd_energies = ccd_energies[mask]
    print(f"Filtered: {n_raw} -> {len(confs)} conformations "
          f"(removed {n_raw - len(confs)} with |E_CCD| >= {energy_cap:.1f})")

    edge_index = torch.tensor(graph_df[["src", "dst"]].values.T, dtype=torch.long)
    num_rcs = meta_df["num_rcs"].values.tolist()
    max_rcs = max(num_rcs)

    # RC features: build padded lookup tables
    aa_table = torch.zeros(num_pos, max_rcs, dtype=torch.long)
    chi_table = torch.zeros(num_pos, max_rcs, MAX_CHI, dtype=torch.float32)

    for p in range(num_pos):
        pos_df = rc_feat_df[rc_feat_df["pos"] == p].sort_values("rc")
        n = len(pos_df)
        aa_table[p, :n] = torch.tensor(pos_df["aa_type_idx"].values, dtype=torch.long)
        chis = pos_df[["chi1", "chi2", "chi3", "chi4"]].fillna(0.0).values
        chi_table[p, :n] = torch.tensor(chis, dtype=torch.float32)

    # Pairwise energy lookup: padded flat table per directed edge
    src_edges, dst_edges = edge_index
    num_edges = src_edges.shape[0]
    pair_table = torch.zeros(num_edges, max_rcs * max_rcs, dtype=torch.float32)

    # Build dict from CSV
    pair_dict = {}
    for _, row in pair_df.iterrows():
        p1, r1, p2, r2 = int(row["pos1"]), int(row["rc1"]), int(row["pos2"]), int(row["rc2"])
        pair_dict[(p1, r1, p2, r2)] = row["E_pair_min"]

    for e in range(num_edges):
        s, d = src_edges[e].item(), dst_edges[e].item()
        ns, nd = num_rcs[s], num_rcs[d]
        for rs in range(ns):
            for rd in range(nd):
                key = (s, rs, d, rd) if (s, rs, d, rd) in pair_dict else (d, rd, s, rs)
                if key in pair_dict:
                    pair_table[e, rs * max_rcs + rd] = pair_dict[key]

    # Cα distances per edge
    ca_dist_map = {}
    for _, row in ca_dist_df.iterrows():
        p1, p2 = int(row["pos1"]), int(row["pos2"])
        ca_dist_map[(p1, p2)] = row["distance"]
        ca_dist_map[(p2, p1)] = row["distance"]

    ca_dist_vec = torch.zeros(num_edges, dtype=torch.float32)
    for e in range(num_edges):
        s, d = src_edges[e].item(), dst_edges[e].item()
        ca_dist_vec[e] = ca_dist_map.get((s, d), 0.0)

    print(f"Loaded {len(confs)} conformations, {num_pos} positions, "
          f"{num_edges} directed edges, max_rcs={max_rcs}")
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

        # Shared AA type embedding
        self.aa_embedding = nn.Embedding(NUM_AA_TYPES, aa_embed_dim)

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

        # Pairwise energy lookup
        flat_idx = rc_src * self.max_rcs + rc_dst  # (batch, num_edges)
        pair_energy = torch.gather(
            self.pair_table.unsqueeze(0).expand(batch, -1, -1),
            2, flat_idx.unsqueeze(-1)
        ).squeeze(-1)  # (batch, num_edges)

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
    print(f"Device: {device}")

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

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_val_rmse = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        model.train()
        train_loss = 0.0
        for confs_batch, targets_batch, _ in train_loader:
            confs_batch = confs_batch.to(device)
            targets_batch = targets_batch.to(device)

            pred = model(confs_batch)
            loss = nn.functional.mse_loss(pred, targets_batch)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            train_loss += loss.item() * len(confs_batch)

        train_rmse = (train_loss / len(train_ds)) ** 0.5
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
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

        if epoch % args.log_every == 0 or epoch == 1:
            print(f"Epoch {epoch:4d}  train_rmse={train_rmse:.4f}  "
                  f"val_rmse={val_rmse:.4f}  best={best_val_rmse:.4f}  "
                  f"lr={scheduler.get_last_lr()[0]:.2e}")

    # Restore best model
    model.load_state_dict(best_state)
    model.eval()

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
    parser.add_argument("--aa_embed_dim", type=int, default=16)
    parser.add_argument("--pos_embed_dim", type=int, default=8)
    parser.add_argument("--node_dim", type=int, default=32)
    parser.add_argument("--hidden_dim", type=int, default=64)
    parser.add_argument("--num_layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--val_frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--log_every", type=int, default=10)
    args = parser.parse_args()

    train(args)
