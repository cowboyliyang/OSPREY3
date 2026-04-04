"""
GNN Energy Surrogate Training Script

Reads data exported by GNNDataExporter (confs.csv, graph.csv, meta.csv),
trains a message-passing GNN to predict E_CCD - E_emat (residual),
and exports the trained model to ONNX for Java inference.

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


# ============================================================
# Data Loading
# ============================================================

class ConfDataset(Dataset):
    """Dataset of conformations with energies."""

    def __init__(self, confs, targets, emat_energies):
        self.confs = torch.tensor(confs, dtype=torch.long)       # (N, num_pos)
        self.targets = torch.tensor(targets, dtype=torch.float32) # (N,) residual
        self.emat_energies = torch.tensor(emat_energies, dtype=torch.float32)

    def __len__(self):
        return len(self.confs)

    def __getitem__(self, idx):
        return self.confs[idx], self.targets[idx], self.emat_energies[idx]


def load_data(data_dir):
    """Load exported CSV files."""
    confs_df = pd.read_csv(os.path.join(data_dir, "confs.csv"))
    graph_df = pd.read_csv(os.path.join(data_dir, "graph.csv"))
    meta_df = pd.read_csv(os.path.join(data_dir, "meta.csv"))

    # Parse conformations (rc_0, rc_1, ..., rc_{n-1})
    num_pos = meta_df.shape[0]
    rc_cols = [f"rc_{i}" for i in range(num_pos)]
    confs = confs_df[rc_cols].values.astype(np.int64)

    # Targets
    residuals = confs_df["residual"].values.astype(np.float64)
    emat_energies = confs_df["E_emat"].values.astype(np.float64)
    ccd_energies = confs_df["E_CCD"].values.astype(np.float64)

    # Filter out extreme-energy conformations (clashes)
    n_raw = len(confs)
    energy_cap = np.percentile(np.abs(ccd_energies), 95)
    energy_cap = max(energy_cap, 500.0)  # at least 500 kcal/mol
    mask = np.abs(ccd_energies) < energy_cap
    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    ccd_energies = ccd_energies[mask]
    print(f"Filtered: {n_raw} -> {len(confs)} conformations "
          f"(removed {n_raw - len(confs)} with |E_CCD| >= {energy_cap:.1f})")

    # Edge index (already has both directions from exporter)
    edge_index = torch.tensor(graph_df[["src", "dst"]].values.T, dtype=torch.long)

    # Per-position RC counts
    num_rcs = meta_df["num_rcs"].values.tolist()

    print(f"Loaded {len(confs)} conformations, {num_pos} positions, "
          f"{edge_index.shape[1]} directed edges")
    print(f"Residual stats: mean={residuals.mean():.4f}, std={residuals.std():.4f}, "
          f"min={residuals.min():.4f}, max={residuals.max():.4f}")
    print(f"E_CCD range: [{ccd_energies.min():.4f}, {ccd_energies.max():.4f}]")

    return confs, residuals, emat_energies, ccd_energies, edge_index, num_rcs


# ============================================================
# GNN Model
# ============================================================

class InteractionGNN(nn.Module):
    """
    Message-passing GNN on the interaction graph.

    Each node = design position. Node features = learned RC embedding.
    Message passing captures many-body interactions via shared parameters.
    Predicts residual = E_CCD - E_emat per conformation.
    """

    def __init__(self, num_rcs_per_pos, embed_dim=32, hidden_dim=64,
                 num_layers=3, dropout=0.1):
        super().__init__()
        self.num_pos = len(num_rcs_per_pos)
        self.embed_dim = embed_dim
        self.num_layers = num_layers

        # Per-position RC embeddings (each position has its own embedding table)
        self.rc_embeddings = nn.ModuleList([
            nn.Embedding(num_rcs, embed_dim) for num_rcs in num_rcs_per_pos
        ])

        # Message passing layers (shared across all edges)
        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()

        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * embed_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, embed_dim),
            ))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * embed_dim, hidden_dim),
                nn.SiLU(),
                nn.Linear(hidden_dim, embed_dim),
            ))
            self.layer_norms.append(nn.LayerNorm(embed_dim))

        # Readout: per-node contribution to residual energy
        self.node_readout = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

        # Edge readout: pairwise correction
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * embed_dim, hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, confs, edge_index):
        """
        Args:
            confs: (batch, num_pos) RC indices per position
            edge_index: (2, num_edges) interaction graph edges

        Returns:
            residual: (batch,) predicted E_CCD - E_emat
        """
        batch_size = confs.shape[0]

        # Embed each position's RC choice: (batch, num_pos, embed_dim)
        h = torch.stack([
            self.rc_embeddings[p](confs[:, p])
            for p in range(self.num_pos)
        ], dim=1)

        src, dst = edge_index  # (num_edges,) each

        # Message passing
        for layer in range(self.num_layers):
            # Gather source and destination embeddings for all edges
            h_src = h[:, src]  # (batch, num_edges, embed_dim)
            h_dst = h[:, dst]  # (batch, num_edges, embed_dim)

            # Compute messages
            msg_input = torch.cat([h_src, h_dst], dim=-1)  # (batch, num_edges, 2*embed_dim)
            messages = self.message_mlps[layer](msg_input)  # (batch, num_edges, embed_dim)

            # Aggregate messages per node (sum)
            agg = torch.zeros_like(h)  # (batch, num_pos, embed_dim)
            dst_expanded = dst.unsqueeze(0).unsqueeze(-1).expand(
                batch_size, -1, self.embed_dim)
            agg.scatter_add_(1, dst_expanded, messages)

            # Update node embeddings with residual connection
            update_input = torch.cat([h, agg], dim=-1)  # (batch, num_pos, 2*embed_dim)
            h = h + self.update_mlps[layer](update_input)
            h = self.layer_norms[layer](h)

        # Readout: node contributions
        node_energy = self.node_readout(h).squeeze(-1)  # (batch, num_pos)
        total_node = node_energy.sum(dim=1)  # (batch,)

        # Readout: edge corrections
        h_src = h[:, src]
        h_dst = h[:, dst]
        edge_input = torch.cat([h_src, h_dst], dim=-1)
        edge_energy = self.edge_readout(edge_input).squeeze(-1)  # (batch, num_edges)
        # Each undirected edge appears twice; divide by 2
        total_edge = edge_energy.sum(dim=1) / 2.0  # (batch,)

        return total_node + total_edge


# ============================================================
# Training
# ============================================================

def train(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load data
    confs, residuals, emat_energies, ccd_energies, edge_index, num_rcs = \
        load_data(args.data_dir)

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

    # Model
    model = InteractionGNN(
        num_rcs_per_pos=num_rcs,
        embed_dim=args.embed_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
    ).to(device)

    num_params = sum(p.numel() for p in model.parameters())
    print(f"Model parameters: {num_params:,}")

    edge_index = edge_index.to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_val_rmse = float("inf")
    best_state = None

    for epoch in range(1, args.epochs + 1):
        # Train
        model.train()
        train_loss = 0.0
        for confs_batch, targets_batch, _ in train_loader:
            confs_batch = confs_batch.to(device)
            targets_batch = targets_batch.to(device)

            pred = model(confs_batch, edge_index)
            loss = nn.functional.mse_loss(pred, targets_batch)

            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            train_loss += loss.item() * len(confs_batch)

        train_rmse = (train_loss / len(train_ds)) ** 0.5
        scheduler.step()

        # Validate
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for confs_batch, targets_batch, _ in val_loader:
                confs_batch = confs_batch.to(device)
                targets_batch = targets_batch.to(device)
                pred = model(confs_batch, edge_index)
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
    all_preds = []
    all_targets = []
    all_emat = []
    with torch.no_grad():
        for confs_batch, targets_batch, emat_batch in val_loader:
            confs_batch = confs_batch.to(device)
            pred = model(confs_batch, edge_index)
            all_preds.append(pred.cpu().numpy())
            all_targets.append(targets_batch.numpy())
            all_emat.append(emat_batch.numpy())

    all_preds = np.concatenate(all_preds)
    all_targets = np.concatenate(all_targets)
    all_emat = np.concatenate(all_emat)

    # === Full evaluation: GNN vs CCD, emat vs CCD ===
    evaluate_all(all_preds, all_targets, all_emat, args)

    # Save
    out_dir = os.path.join(args.data_dir, "model")
    os.makedirs(out_dir, exist_ok=True)

    # PyTorch checkpoint
    torch.save({
        "model_state": best_state,
        "num_rcs": num_rcs,
        "edge_index": edge_index.cpu(),
        "args": vars(args),
        "val_rmse": best_val_rmse,
    }, os.path.join(out_dir, "gnn_checkpoint.pt"))

    # ONNX export (optional, requires onnx package)
    try:
        export_onnx(model, edge_index, num_rcs, out_dir, device)
    except Exception as e:
        print(f"ONNX export skipped: {e}")

    print(f"\nModel saved to {out_dir}")


def spearman(x, y):
    """Spearman rank correlation."""
    from scipy.stats import spearmanr
    rho, _ = spearmanr(x, y)
    return rho


def top_k_recall(pred, true, k):
    """Fraction of true top-k that appear in predicted top-k."""
    k = min(k, len(pred))
    true_topk = set(np.argsort(true)[:k])
    pred_topk = set(np.argsort(pred)[:k])
    return len(true_topk & pred_topk) / k


def evaluate_all(residual_preds, residual_targets, emat_energies, args):
    """Full comparison: GNN vs CCD and emat vs CCD."""
    n = len(residual_preds)

    # GNN predicted total energy
    gnn_total = emat_energies + residual_preds
    # CCD true total energy
    ccd_total = emat_energies + residual_targets
    # emat baseline (residual = 0)
    emat_total = emat_energies

    # --- GNN vs CCD ---
    gnn_delta = gnn_total - ccd_total
    gnn_mae = np.mean(np.abs(gnn_delta))
    gnn_rmse = np.sqrt(np.mean(gnn_delta ** 2))
    gnn_bias = np.mean(gnn_delta)
    gnn_max = np.max(np.abs(gnn_delta))

    # --- emat vs CCD (baseline) ---
    emat_delta = emat_total - ccd_total
    emat_mae = np.mean(np.abs(emat_delta))
    emat_rmse = np.sqrt(np.mean(emat_delta ** 2))
    emat_bias = np.mean(emat_delta)
    emat_max = np.max(np.abs(emat_delta))

    # Rank correlation and top-K
    try:
        gnn_spearman = spearman(gnn_total, ccd_total)
        emat_spearman = spearman(emat_total, ccd_total)
    except Exception:
        gnn_spearman = emat_spearman = float("nan")

    top_ks = [5, 10, 20]

    print("\n" + "=" * 70)
    print("GNN vs CCD  (validation set, N=%d)" % n)
    print("=" * 70)
    print(f"  MAE:      {gnn_mae:.6f} kcal/mol")
    print(f"  RMSE:     {gnn_rmse:.6f} kcal/mol")
    print(f"  Bias:     {gnn_bias:.6f} kcal/mol")
    print(f"  Max|d|:   {gnn_max:.6f} kcal/mol")
    print(f"  Spearman: {gnn_spearman:.4f}")
    for k in top_ks:
        if k <= n:
            r = top_k_recall(gnn_total, ccd_total, k)
            print(f"  Top-{k} recall: {r * 100:.1f}%")

    print("\n" + "=" * 70)
    print("emat vs CCD  (baseline, N=%d)" % n)
    print("=" * 70)
    print(f"  MAE:      {emat_mae:.6f} kcal/mol")
    print(f"  RMSE:     {emat_rmse:.6f} kcal/mol")
    print(f"  Bias:     {emat_bias:.6f} kcal/mol")
    print(f"  Max|d|:   {emat_max:.6f} kcal/mol")
    print(f"  Spearman: {emat_spearman:.4f}")
    for k in top_ks:
        if k <= n:
            r = top_k_recall(emat_total, ccd_total, k)
            print(f"  Top-{k} recall: {r * 100:.1f}%")

    print("\n" + "=" * 70)
    print("Improvement: GNN over emat baseline")
    print("=" * 70)
    print(f"  MAE  reduction: {(1 - gnn_mae / emat_mae) * 100:.1f}%")
    print(f"  RMSE reduction: {(1 - gnn_rmse / emat_rmse) * 100:.1f}%")

    # Per-conf table (top 30 by CCD energy)
    sort_idx = np.argsort(ccd_total)[:30]
    print("\n" + "=" * 70)
    print("Per-conf comparison (top 30 lowest-energy conformations)")
    print("=" * 70)
    print(f"{'Rank':<5} {'E_CCD':>11} {'E_GNN':>11} {'E_emat':>11} "
          f"{'GNN_err':>9} {'emat_err':>9}")
    print("-" * 70)
    for rank, idx in enumerate(sort_idx):
        print(f"{rank:<5d} {ccd_total[idx]:11.4f} {gnn_total[idx]:11.4f} "
              f"{emat_total[idx]:11.4f} {gnn_delta[idx]:9.4f} "
              f"{emat_delta[idx]:9.4f}")


def export_onnx(model, edge_index, num_rcs, out_dir, device):
    """Export model to ONNX for Java inference."""
    model.eval()
    num_pos = len(num_rcs)

    # Dummy input: single conformation
    dummy_confs = torch.zeros(1, num_pos, dtype=torch.long, device=device)

    # ONNX doesn't handle dynamic edge_index well, so we bake it in
    # by wrapping the model
    class OnnxWrapper(nn.Module):
        def __init__(self, gnn, edge_idx):
            super().__init__()
            self.gnn = gnn
            self.register_buffer("edge_index", edge_idx)

        def forward(self, confs):
            return self.gnn(confs, self.edge_index)

    wrapper = OnnxWrapper(model, edge_index).to(device)
    wrapper.eval()

    onnx_path = os.path.join(out_dir, "gnn_model.onnx")
    torch.onnx.export(
        wrapper,
        dummy_confs,
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
    parser = argparse.ArgumentParser(description="Train GNN energy surrogate")
    parser.add_argument("--data_dir", type=str, required=True,
                        help="Directory with confs.csv, graph.csv, meta.csv")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--embed_dim", type=int, default=32)
    parser.add_argument("--hidden_dim", type=int, default=64)
    parser.add_argument("--num_layers", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--val_frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--log_every", type=int, default=10)
    args = parser.parse_args()

    train(args)
