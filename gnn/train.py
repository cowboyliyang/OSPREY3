"""
GNN Energy Surrogate Training Script v2

Improved architecture with:
  - Shared AA type embeddings (instead of position-specific RC embeddings)
  - Chi angle encoding (sin/cos of rotamer dihedral angles)
  - Position embeddings
  - Edge features from emat pairwise energies + Ca distances

Reads data exported by GNNDataExporter:
  confs.csv, graph.csv, meta.csv, rc_features.csv, pairwise_energies.csv, ca_distances.csv
  (optional) onebody_energies.csv  -- fixed-environment (emat one-body) node feature

Usage:
    python train.py --data_dir gnn_data/2RL0_flex8 --epochs 200

DANCE model-program additions (all OPT-IN; defaults reproduce the original v2):
  --loss zaware         Boltzmann-weighted, anti-optimism asymmetric loss
                        (Z-functional surrogate; weights each conf by its
                         Boltzmann mass, penalizes ENERGY under-prediction more).
  --bracket_floor       softplus the residual so E_GNN >= E_emat (never below the
                        provable lower bound -> structurally cannot fabricate mass).
  --bracket_ceil_lambda soft penalty pushing E_GNN <= E_rigid (within the bracket).
  --uncertainty         heteroscedastic sigma_hat head (2nd ONNX output "log_var");
                        online Java still reads only "residual", so it is backward
                        compatible. sigma_hat is for offline audit routing / CP.
  --no_pos_embed        drop the absolute-position embedding crutch.
  --onebody_feat        add the emat one-body energy as a (transferable) node
                        feature encoding the fixed scaffold environment.
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

# Boltzmann kT at the OSPREY design temperature (298 K), kcal/mol.
# Consistent with log_zhat = -E/kT in the leaf audit logs.
KT_KCAL = 0.5924

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
    """Carries everything the DANCE losses need per conformation."""

    def __init__(self, confs, targets, emat_energies, ccd_energies, rigid_energies):
        self.confs = torch.tensor(confs, dtype=torch.long)
        self.targets = torch.tensor(targets, dtype=torch.float32)
        self.emat_energies = torch.tensor(emat_energies, dtype=torch.float32)
        self.ccd_energies = torch.tensor(ccd_energies, dtype=torch.float32)
        self.rigid_energies = torch.tensor(rigid_energies, dtype=torch.float32)

    def __len__(self):
        return len(self.confs)

    def __getitem__(self, idx):
        return (self.confs[idx], self.targets[idx], self.emat_energies[idx],
                self.ccd_energies[idx], self.rigid_energies[idx])


def load_data(data_dir):
    """Load all exported CSV files."""
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

    onebody_path = os.path.join(data_dir, "onebody_energies.csv")
    onebody_df = pd.read_csv(onebody_path) if os.path.exists(onebody_path) else None

    num_pos = meta_df.shape[0]
    rc_cols = [f"rc_{i}" for i in range(num_pos)]
    confs = confs_df[rc_cols].values.astype(np.int64)

    residuals = confs_df["residual"].values.astype(np.float64)
    emat_energies = confs_df["E_emat"].values.astype(np.float64)
    ccd_energies = confs_df["E_CCD"].values.astype(np.float64)
    has_rigid = "E_rigid" in confs_df.columns
    rigid_energies = (confs_df["E_rigid"].values.astype(np.float64) if has_rigid
                      else np.full(len(confs), np.nan, dtype=np.float64))

    del confs_df

    n_raw = len(confs)
    mask = np.ones(n_raw, dtype=bool)

    emat_cap = -20.0
    mask_emat = emat_energies <= emat_cap
    mask &= mask_emat

    if has_rigid:
        rigid_cap = 0.0
        mask_rigid = rigid_energies <= rigid_cap
        mask &= mask_rigid

    residual_cap = 100.0
    mask_res = np.abs(residuals) < residual_cap
    mask &= mask_res

    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    ccd_energies = ccd_energies[mask]
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

    onebody_table = torch.zeros(num_pos, max_rcs, dtype=torch.float32)
    if onebody_df is not None:
        ob_p = onebody_df["pos"].values
        ob_r = onebody_df["rc"].values
        ob_e = onebody_df["E_onebody"].values
        for i in range(len(ob_e)):
            if 0 <= ob_p[i] < num_pos and 0 <= ob_r[i] < max_rcs:
                onebody_table[ob_p[i], ob_r[i]] = ob_e[i]

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
    print(f"Residual stats: mean={residuals.mean():.4f}, std={residuals.std():.4f}, "
          f"min={residuals.min():.4f}, max={residuals.max():.4f}")
    if onebody_df is not None:
        print(f"One-body env feature: loaded {len(onebody_df)} (pos,rc) entries")
    else:
        print("One-body env feature: onebody_energies.csv not found (zeros)")

    return {
        "confs": confs,
        "residuals": residuals,
        "emat_energies": emat_energies,
        "ccd_energies": ccd_energies,
        "rigid_energies": rigid_energies,
        "has_rigid": has_rigid,
        "edge_index": edge_index,
        "num_rcs": num_rcs,
        "num_pos": num_pos,
        "max_rcs": max_rcs,
        "aa_table": aa_table,
        "chi_table": chi_table,
        "onebody_table": onebody_table,
        "pair_table": pair_table,
        "ca_dist_vec": ca_dist_vec,
    }


# ============================================================
# GNN Model v2
# ============================================================

class InteractionGNNv2(nn.Module):
    """Pairwise-energy-aware message-passing GNN with opt-in DANCE additions."""

    def __init__(self, num_pos, max_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, onebody_table=None,
                 aa_embed_dim=16, pos_embed_dim=8, node_dim=32,
                 hidden_dim=64, num_layers=3, dropout=0.1,
                 use_pos_embed=True, use_onebody=False,
                 bracket_floor=False, uncertainty=False):
        super().__init__()
        self.num_pos = num_pos
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.max_rcs = max_rcs
        self.use_pos_embed = use_pos_embed
        self.use_onebody = use_onebody
        self.bracket_floor = bracket_floor
        self.uncertainty = uncertainty

        edge_feat_dim = 2

        self.aa_embedding = nn.Embedding(NUM_AA_TYPES + 1, aa_embed_dim)

        if use_pos_embed:
            self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        chi_input_dim = MAX_CHI * 2

        node_input_dim = aa_embed_dim + chi_input_dim
        if use_pos_embed:
            node_input_dim += pos_embed_dim
        if use_onebody:
            node_input_dim += 1

        self.node_proj = nn.Sequential(
            nn.Linear(node_input_dim, node_dim),
            nn.SiLU(),
        )

        self.register_buffer("aa_table", aa_table)
        self.register_buffer("chi_table", chi_table)
        self.register_buffer("pair_table", pair_table)
        self.register_buffer("ca_dist_vec", ca_dist_vec)
        self.register_buffer("edge_src", edge_index[0])
        self.register_buffer("edge_dst", edge_index[1])
        if onebody_table is None:
            onebody_table = torch.zeros(num_pos, max_rcs, dtype=torch.float32)
        self.register_buffer("onebody_table", onebody_table)

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

        if uncertainty:
            self.var_readout = nn.Sequential(
                nn.Linear(node_dim, hidden_dim),
                nn.SiLU(),
                nn.Dropout(dropout),
                nn.Linear(hidden_dim, 1),
            )

    def forward(self, confs):
        batch = confs.shape[0]
        num_pos = self.num_pos
        node_dim = self.node_dim

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

        node_parts = [aa_feat, chi_feat]

        if self.use_pos_embed:
            pos_idx = torch.arange(num_pos, device=confs.device)
            pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)
            node_parts.append(pos_feat)

        if self.use_onebody:
            ob = torch.gather(
                self.onebody_table.unsqueeze(0).expand(batch, -1, -1),
                2, confs.unsqueeze(-1)
            )
            node_parts.append(ob)

        node_input = torch.cat(node_parts, dim=-1)
        h = self.node_proj(node_input)

        src = self.edge_src
        dst = self.edge_dst
        num_edges = src.shape[0]

        src_exp = src.unsqueeze(0).expand(batch, -1)
        dst_exp = dst.unsqueeze(0).expand(batch, -1)
        rc_src = torch.gather(confs, 1, src_exp)
        rc_dst = torch.gather(confs, 1, dst_exp)

        flat_idx = rc_src * self.max_rcs + rc_dst
        edge_arange = torch.arange(num_edges, device=confs.device)
        pair_energy = self.pair_table[edge_arange, flat_idx]

        ca_dist = self.ca_dist_vec.unsqueeze(0).expand(batch, -1)

        edge_feats = torch.stack([pair_energy, ca_dist], dim=-1)

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

        node_energy = self.node_readout(h).squeeze(-1)
        total_node = node_energy.sum(dim=1)

        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)
        edge_energy = self.edge_readout(edge_input).squeeze(-1)
        total_edge = edge_energy.sum(dim=1) / 2.0

        residual = total_node + total_edge

        if self.bracket_floor:
            residual = nn.functional.softplus(residual)

        if self.uncertainty:
            graph_embed = h.mean(dim=1)
            log_var = self.var_readout(graph_embed).squeeze(-1)
            return residual, log_var
        return residual


# ============================================================
# DANCE loss
# ============================================================

def dance_loss(pred, target, ccd, emat, rigid, args):
    err = pred - target  # >0: energy too HIGH; <0: too LOW (optimistic)

    delta = args.huber_delta
    abse = err.abs()
    base = torch.where(abse <= delta, 0.5 * err * err, delta * (abse - 0.5 * delta))

    asym = torch.where(err < 0, float(args.asym_under), 1.0)

    if args.boltzmann_lambda > 0.0:
        shifted = -(ccd - ccd.min()) / (KT_KCAL * args.boltzmann_temp)
        bw = torch.softmax(shifted, dim=0)
        nb = pred.shape[0]
        w = (1.0 - args.boltzmann_lambda) / nb + args.boltzmann_lambda * bw
    else:
        w = torch.full_like(pred, 1.0 / pred.shape[0])

    if args.uncertainty_logvar is not None:
        log_var = args.uncertainty_logvar
        nll = 0.5 * torch.exp(-log_var) * (err * err) + 0.5 * log_var
        per_conf = asym * nll
    else:
        per_conf = asym * base

    loss = (w * per_conf).sum()

    if args.bracket_ceil_lambda > 0.0 and torch.isfinite(rigid).all():
        ceil = rigid - emat
        over = torch.relu(pred - ceil)
        loss = loss + args.bracket_ceil_lambda * (over * over).mean()

    return loss


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
    ccd_energies = data["ccd_energies"]
    rigid_energies = data["rigid_energies"]

    n = len(confs)
    perm = np.random.RandomState(args.seed).permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    train_ds = ConfDataset(confs[train_idx], residuals[train_idx], emat_energies[train_idx],
                           ccd_energies[train_idx], rigid_energies[train_idx])
    val_ds = ConfDataset(confs[val_idx], residuals[val_idx], emat_energies[val_idx],
                         ccd_energies[val_idx], rigid_energies[val_idx])
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size)

    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")
    print(f"Config: loss={args.loss} pos_embed={not args.no_pos_embed} "
          f"onebody={args.onebody_feat} bracket_floor={args.bracket_floor} "
          f"uncertainty={args.uncertainty} asym_under={args.asym_under} "
          f"boltzmann_lambda={args.boltzmann_lambda}")

    model = InteractionGNNv2(
        num_pos=data["num_pos"],
        max_rcs=data["max_rcs"],
        aa_table=data["aa_table"],
        chi_table=data["chi_table"],
        pair_table=data["pair_table"],
        edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        onebody_table=data["onebody_table"],
        aa_embed_dim=args.aa_embed_dim,
        pos_embed_dim=args.pos_embed_dim,
        node_dim=args.node_dim,
        hidden_dim=args.hidden_dim,
        num_layers=args.num_layers,
        dropout=args.dropout,
        use_pos_embed=not args.no_pos_embed,
        use_onebody=args.onebody_feat,
        bracket_floor=args.bracket_floor,
        uncertainty=args.uncertainty,
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

    use_dance = (args.loss == "zaware") or args.bracket_ceil_lambda > 0 or args.uncertainty

    for epoch in range(start_epoch, args.epochs + 1):
        model.train()
        train_mse = 0.0
        for confs_batch, targets_batch, emat_batch, ccd_batch, rigid_batch in train_loader:
            confs_batch = confs_batch.to(device)
            targets_batch = targets_batch.to(device)

            out = model(confs_batch)
            if args.uncertainty:
                pred, log_var = out
                args.uncertainty_logvar = log_var
            else:
                pred = out
                args.uncertainty_logvar = None

            if use_dance:
                loss = dance_loss(
                    pred, targets_batch,
                    ccd_batch.to(device), emat_batch.to(device), rigid_batch.to(device),
                    args,
                )
            elif args.loss == "huber":
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
            for confs_batch, targets_batch, emat_batch, ccd_batch, rigid_batch in val_loader:
                confs_batch = confs_batch.to(device)
                targets_batch = targets_batch.to(device)
                out = model(confs_batch)
                pred = out[0] if args.uncertainty else out
                val_loss += nn.functional.mse_loss(pred, targets_batch).item() * len(confs_batch)
        val_rmse = (val_loss / len(val_ds)) ** 0.5

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            raw_sd = model.module.state_dict() if isinstance(model, nn.DataParallel) else model.state_dict()
            best_state = {k: v.cpu().clone() for k, v in raw_sd.items()}

        if epoch % args.log_every == 0 or epoch == start_epoch:
            print(f"Epoch {epoch:4d}  train_rmse={train_rmse:.4f}  "
                  f"val_rmse={val_rmse:.4f}  best={best_val_rmse:.4f}  "
                  f"lr={scheduler.get_last_lr()[0]:.2e}", flush=True)

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

    base_model = model.module if isinstance(model, nn.DataParallel) else model
    base_model.load_state_dict(best_state)
    base_model.eval()
    model = base_model

    all_preds, all_targets, all_emat, all_logvar = [], [], [], []
    with torch.no_grad():
        for confs_batch, targets_batch, emat_batch, ccd_batch, rigid_batch in val_loader:
            confs_batch = confs_batch.to(device)
            out = model(confs_batch)
            if args.uncertainty:
                pred, log_var = out
                all_logvar.append(log_var.cpu().numpy())
            else:
                pred = out
            all_preds.append(pred.cpu().numpy())
            all_targets.append(targets_batch.numpy())
            all_emat.append(emat_batch.numpy())

    all_preds = np.concatenate(all_preds)
    all_targets = np.concatenate(all_targets)
    all_emat = np.concatenate(all_emat)

    evaluate_all(all_preds, all_targets, all_emat)

    val_errors = np.abs(all_preds - all_targets)
    n_cal = len(val_errors)
    cp_stats = {}
    for alpha in [0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001]:
        level = min(np.ceil((1 - alpha) * (n_cal + 1)) / n_cal, 1.0)
        q = float(np.quantile(val_errors, level))
        cp_stats[f"q_alpha_{alpha}"] = q
        print(f"  CP calibration: a={alpha:.3f}  q={q:.6f} kcal/mol  "
              f"(n_cal={n_cal}, coverage>={1-alpha:.3f})")

    for p in [50, 90, 95, 99, 99.5, 99.9]:
        cp_stats[f"P{p}"] = float(np.percentile(val_errors, p))
    cp_stats["max"] = float(np.max(val_errors))
    cp_stats["n_cal"] = n_cal

    print(f"\n  Val error percentiles:")
    print(f"    P50={cp_stats['P50']:.6f}  P90={cp_stats['P90']:.6f}  "
          f"P95={cp_stats['P95']:.6f}  P99={cp_stats['P99']:.6f}  "
          f"Max={cp_stats['max']:.6f}")

    if args.uncertainty and len(all_logvar) > 0:
        sigma = np.exp(0.5 * np.concatenate(all_logvar))
        try:
            from scipy.stats import spearmanr
            rho_u, _ = spearmanr(sigma, val_errors)
            print(f"  sigma_hat vs |error| Spearman: {rho_u:.4f} (audit-routing signal)")
            cp_stats["sigma_err_spearman"] = float(rho_u)
        except Exception:
            pass

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
        export_onnx(model, data["num_pos"], out_dir, device, uncertainty=args.uncertainty)
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


def export_onnx(model, num_pos, out_dir, device, uncertainty=False):
    model.eval()
    dummy = torch.zeros(1, num_pos, dtype=torch.long, device=device)

    if uncertainty:
        output_names = ["residual", "log_var"]
        dynamic_axes = {"confs": {0: "batch"},
                        "residual": {0: "batch"}, "log_var": {0: "batch"}}
    else:
        output_names = ["residual"]
        dynamic_axes = {"confs": {0: "batch"}, "residual": {0: "batch"}}

    onnx_path = os.path.join(out_dir, "gnn_model.onnx")
    torch.onnx.export(
        model,
        dummy,
        onnx_path,
        input_names=["confs"],
        output_names=output_names,
        dynamic_axes=dynamic_axes,
        opset_version=17,
    )
    print(f"ONNX model exported to {onnx_path} (outputs={output_names})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train GNN energy surrogate v2")
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--loss", type=str, default="huber", choices=["mse", "huber", "zaware"])
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
    parser.add_argument("--resume", type=str, default=None)
    parser.add_argument("--boltzmann_lambda", type=float, default=0.0)
    parser.add_argument("--boltzmann_temp", type=float, default=1.0)
    parser.add_argument("--asym_under", type=float, default=1.0)
    parser.add_argument("--bracket_floor", action="store_true")
    parser.add_argument("--bracket_ceil_lambda", type=float, default=0.0)
    parser.add_argument("--uncertainty", action="store_true")
    parser.add_argument("--no_pos_embed", action="store_true")
    parser.add_argument("--onebody_feat", action="store_true")
    args = parser.parse_args()
    args.uncertainty_logvar = None

    train(args)
