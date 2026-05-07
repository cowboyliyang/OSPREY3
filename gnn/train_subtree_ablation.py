"""
Subtree GNN Architecture Ablation

Architectures:
  baseline  — Current SubtreeGNN (MLP message passing, is_fixed flag, free edge=0)
  gat       — Graph Attention Network message passing
  dual_enc  — Separate encoders for fixed/free nodes
  multihead — Multi-head readout by #free positions
  transformer — Graph Transformer (full self-attention)
  transfer  — Transfer learning from per-conf checkpoint

Usage:
    python train_subtree_ablation.py --data_dir ... --arch gat
"""

import argparse
import os
import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
from itertools import combinations


NUM_AA_TYPES = 20
MAX_CHI = 4
KT = 0.5922

AA_TYPES = {
    "ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4,
    "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9,
    "LEU": 10, "LYS": 11, "MET": 12, "PHE": 13, "PRO": 14,
    "SER": 15, "THR": 16, "TRP": 17, "TYR": 18, "VAL": 19,
}


# ============================================================
# Data (shared across all architectures)
# ============================================================

class SubtreeDataset(Dataset):
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
    rng = np.random.RandomState(seed)
    all_fixed_rcs, all_free_masks, all_labels = [], [], []

    print("Generating subtree samples from conformations...")
    for k in range(1, num_pos):
        all_combos = list(combinations(range(num_pos), k))
        if len(all_combos) > 100:
            all_combos = [all_combos[i] for i in rng.choice(len(all_combos), 100, replace=False)]

        k_count = 0
        for fixed_pos in all_combos:
            fixed_pos = list(fixed_pos)
            keys = confs[:, fixed_pos]
            key_dtype = np.dtype([(f"p{i}", np.int64) for i in range(k)])
            key_arr = np.array([tuple(row) for row in keys], dtype=key_dtype)
            unique_keys, inverse, counts = np.unique(key_arr, return_inverse=True, return_counts=True)

            for gid in range(len(unique_keys)):
                if counts[gid] < min_confs:
                    continue
                mask = inverse == gid
                g_emat = emat_energies[mask]
                g_res = residuals[mask]
                neg_e_kt = -g_emat / KT
                neg_e_kt -= neg_e_kt.max()
                w = np.exp(neg_e_kt)
                w /= w.sum()
                label = float((w * g_res).sum())

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
        print(f"  fix {k}/{num_pos} positions ({len(all_combos)} combos): {k_count} subtree samples")

    # Full confs — no cap
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

    n_raw = len(confs)
    mask = np.ones(n_raw, dtype=bool)
    mask &= emat_energies <= -20.0
    mask &= np.abs(residuals) < 100.0
    if has_rigid:
        mask &= rigid_energies <= 0.0
    confs = confs[mask]
    residuals = residuals[mask]
    emat_energies = emat_energies[mask]
    print(f"Filtered: {n_raw} -> {len(confs)} conformations")

    edge_index = torch.tensor(graph_df[["src", "dst"]].values.T, dtype=torch.long)
    num_rcs = meta_df["num_rcs"].values.tolist()
    max_rcs = max(num_rcs)

    aa_table = torch.zeros(num_pos, max_rcs, dtype=torch.long)
    chi_table = torch.zeros(num_pos, max_rcs, MAX_CHI, dtype=torch.float32)
    rc_feat_df = rc_feat_df.sort_values(["pos", "rc"])
    for p in range(num_pos):
        pos_data = rc_feat_df[rc_feat_df["pos"].values == p]
        n = len(pos_data)
        aa_table[p, :n] = torch.tensor(pos_data["aa_type_idx"].values, dtype=torch.long)
        chis = pos_data[["chi1", "chi2", "chi3", "chi4"]].fillna(0.0).values
        chi_table[p, :n] = torch.tensor(chis, dtype=torch.float32)

    src_edges, dst_edges = edge_index
    num_edges = src_edges.shape[0]
    pair_table = torch.zeros(num_edges, max_rcs * max_rcs, dtype=torch.float32)
    edge_map = {}
    for e in range(num_edges):
        edge_map[(src_edges[e].item(), dst_edges[e].item())] = e

    evals = pair_df["E_pair_min"].values
    p1, r1 = pair_df["pos1"].values, pair_df["rc1"].values
    p2, r2 = pair_df["pos2"].values, pair_df["rc2"].values
    for i in range(len(evals)):
        e_idx = edge_map.get((p1[i], p2[i]))
        if e_idx is not None:
            pair_table[e_idx, r1[i] * max_rcs + r2[i]] = evals[i]
        e_idx = edge_map.get((p2[i], p1[i]))
        if e_idx is not None:
            pair_table[e_idx, r2[i] * max_rcs + r1[i]] = evals[i]

    # Clamp inf
    pair_table.clamp_(max=100.0)

    ca_dist_vec = torch.zeros(num_edges, dtype=torch.float32)
    ca_d = ca_dist_df["distance"].values
    ca_p1, ca_p2 = ca_dist_df["pos1"].values, ca_dist_df["pos2"].values
    ca_dist_map = {}
    for i in range(len(ca_d)):
        ca_dist_map[(ca_p1[i], ca_p2[i])] = ca_d[i]
        ca_dist_map[(ca_p2[i], ca_p1[i])] = ca_d[i]
    for e in range(num_edges):
        ca_dist_vec[e] = ca_dist_map.get((src_edges[e].item(), dst_edges[e].item()), 0.0)

    elapsed = _time.time() - t0
    print(f"Loaded {len(confs)} confs, {num_pos} pos, {num_edges} edges, max_rcs={max_rcs} in {elapsed:.1f}s")

    return {
        "confs": confs, "residuals": residuals, "emat_energies": emat_energies,
        "edge_index": edge_index, "num_rcs": num_rcs, "num_pos": num_pos,
        "max_rcs": max_rcs, "aa_table": aa_table, "chi_table": chi_table,
        "pair_table": pair_table, "ca_dist_vec": ca_dist_vec,
    }


# ============================================================
# Shared: node/edge feature extraction mixin
# ============================================================

class SubtreeFeatureMixin:
    """Shared logic for extracting node & edge features from (fixed_rcs, free_mask)."""

    def _init_features(self, num_pos, max_rcs, num_rcs, aa_table, chi_table,
                       pair_table, edge_index, ca_dist_vec,
                       aa_embed_dim, pos_embed_dim):
        self.num_pos = num_pos
        self.max_rcs = max_rcs
        self.aa_embed_dim = aa_embed_dim
        chi_input_dim = MAX_CHI * 2

        self.aa_embedding = nn.Embedding(NUM_AA_TYPES, aa_embed_dim)
        self.pos_embedding = nn.Embedding(num_pos, pos_embed_dim)

        self.register_buffer("aa_table", aa_table)
        self.register_buffer("chi_table", chi_table)
        self.register_buffer("pair_table", pair_table)
        self.register_buffer("ca_dist_vec", ca_dist_vec)
        self.register_buffer("edge_src", edge_index[0])
        self.register_buffer("edge_dst", edge_index[1])
        self.register_buffer("num_rcs_t", torch.tensor(num_rcs, dtype=torch.long))

        # Precompute mean chi features
        mean_chi = torch.zeros(num_pos, chi_input_dim)
        for p in range(num_pos):
            nr = num_rcs[p]
            chi_raw = chi_table[p, :nr]
            chi_rad = chi_raw * (math.pi / 180.0)
            chi_sc = torch.cat([torch.sin(chi_rad), torch.cos(chi_rad)], dim=-1)
            mean_chi[p] = chi_sc.mean(dim=0)
        self.register_buffer("mean_chi_feat", mean_chi)

        self._feat_input_dim = aa_embed_dim + chi_input_dim + pos_embed_dim + 1  # +1 for is_fixed

    def _get_mean_aa_feat(self):
        feats = []
        for p in range(self.num_pos):
            nr = self.num_rcs_t[p].item()
            aa_embeds = self.aa_embedding(self.aa_table[p, :nr])
            feats.append(aa_embeds.mean(dim=0))
        return torch.stack(feats, dim=0)

    def _extract_features(self, fixed_rcs, free_mask):
        """Returns: node_input (B, P, feat_dim), edge_feats (B, E, 2), h_is_fixed (B, P, 1)"""
        batch = fixed_rcs.shape[0]
        device = fixed_rcs.device
        num_pos = self.num_pos

        # Node features
        aa_idx_fixed = torch.gather(
            self.aa_table.unsqueeze(0).expand(batch, -1, -1), 2, fixed_rcs.unsqueeze(-1)
        ).squeeze(-1)
        aa_feat_fixed = self.aa_embedding(aa_idx_fixed)

        chi_raw_fixed = torch.gather(
            self.chi_table.unsqueeze(0).expand(batch, -1, -1, -1),
            2, fixed_rcs.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 1, MAX_CHI)
        ).squeeze(2)
        chi_rad_fixed = chi_raw_fixed * (math.pi / 180.0)
        chi_feat_fixed = torch.cat([torch.sin(chi_rad_fixed), torch.cos(chi_rad_fixed)], dim=-1)

        mean_aa_feat = self._get_mean_aa_feat()
        mask_3d = free_mask.unsqueeze(-1)
        aa_feat = torch.where(mask_3d, mean_aa_feat.unsqueeze(0).expand(batch, -1, -1), aa_feat_fixed)
        chi_feat = torch.where(mask_3d, self.mean_chi_feat.unsqueeze(0).expand(batch, -1, -1), chi_feat_fixed)

        is_fixed = (~free_mask).float().unsqueeze(-1)
        pos_idx = torch.arange(num_pos, device=device)
        pos_feat = self.pos_embedding(pos_idx).unsqueeze(0).expand(batch, -1, -1)

        node_input = torch.cat([aa_feat, chi_feat, pos_feat, is_fixed], dim=-1)

        # Edge features
        src, dst = self.edge_src, self.edge_dst
        src_exp = src.unsqueeze(0).expand(batch, -1)
        dst_exp = dst.unsqueeze(0).expand(batch, -1)

        src_free = torch.gather(free_mask, 1, src_exp)
        dst_free = torch.gather(free_mask, 1, dst_exp)
        either_free = src_free | dst_free

        rc_src = torch.gather(fixed_rcs, 1, src_exp)
        rc_dst = torch.gather(fixed_rcs, 1, dst_exp)
        flat_idx = rc_src * self.max_rcs + rc_dst
        pair_energy_fixed = torch.gather(
            self.pair_table.unsqueeze(0).expand(batch, -1, -1), 2, flat_idx.unsqueeze(-1)
        ).squeeze(-1)
        pair_energy = torch.where(either_free, torch.zeros_like(pair_energy_fixed), pair_energy_fixed)
        ca_dist = self.ca_dist_vec.unsqueeze(0).expand(batch, -1)
        edge_feats = torch.stack([pair_energy, ca_dist], dim=-1)

        return node_input, edge_feats, is_fixed


# ============================================================
# 1. Baseline — same as train_subtree.py
# ============================================================

class SubtreeBaseline(SubtreeFeatureMixin, nn.Module):
    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        edge_feat_dim = 2

        self.node_proj = nn.Sequential(nn.Linear(self._feat_input_dim, node_dim), nn.SiLU())

        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1))

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, _ = self._extract_features(fixed_rcs, free_mask)
        h = self.node_proj(node_input)

        src, dst = self.edge_src, self.edge_dst
        node_dim = self.node_dim
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)
            h_dst = torch.gather(h, 1, dst_idx)
            messages = self.message_mlps[layer](torch.cat([h_src, h_dst, edge_feats], dim=-1))
            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)
            h = h + self.update_mlps[layer](torch.cat([h, agg], dim=-1))
            h = self.layer_norms[layer](h)

        node_energy = self.node_readout(h).squeeze(-1).sum(dim=1)
        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_energy = self.edge_readout(torch.cat([h_src, h_dst, edge_feats], dim=-1)).squeeze(-1).sum(dim=1) / 2.0
        return node_energy + edge_energy


# ============================================================
# 2. GAT — Graph Attention message passing
# ============================================================

class SubtreeGAT(SubtreeFeatureMixin, nn.Module):
    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1,
                 num_heads=4, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.num_heads = num_heads
        assert node_dim % num_heads == 0
        head_dim = node_dim // num_heads
        edge_feat_dim = 2

        self.node_proj = nn.Sequential(nn.Linear(self._feat_input_dim, node_dim), nn.SiLU())

        # GAT layers
        self.W_q = nn.ModuleList()
        self.W_k = nn.ModuleList()
        self.W_v = nn.ModuleList()
        self.W_e = nn.ModuleList()  # edge feature projection
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()

        for _ in range(num_layers):
            self.W_q.append(nn.Linear(node_dim, node_dim, bias=False))
            self.W_k.append(nn.Linear(node_dim, node_dim, bias=False))
            self.W_v.append(nn.Linear(node_dim, node_dim, bias=False))
            self.W_e.append(nn.Linear(edge_feat_dim, num_heads, bias=False))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1))

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, _ = self._extract_features(fixed_rcs, free_mask)
        h = self.node_proj(node_input)

        src, dst = self.edge_src, self.edge_dst
        num_edges = src.shape[0]
        node_dim = self.node_dim
        num_heads = self.num_heads
        head_dim = node_dim // num_heads
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            q = self.W_q[layer](h).view(batch, self.num_pos, num_heads, head_dim)
            k = self.W_k[layer](h).view(batch, self.num_pos, num_heads, head_dim)
            v = self.W_v[layer](h).view(batch, self.num_pos, num_heads, head_dim)

            # Gather for edges
            src_idx_4d = src.unsqueeze(0).unsqueeze(-1).unsqueeze(-1).expand(batch, -1, num_heads, head_dim)
            dst_idx_4d = dst.unsqueeze(0).unsqueeze(-1).unsqueeze(-1).expand(batch, -1, num_heads, head_dim)

            q_dst = torch.gather(q, 1, dst_idx_4d)  # (B, E, H, D)
            k_src = torch.gather(k, 1, src_idx_4d)
            v_src = torch.gather(v, 1, src_idx_4d)

            # Attention scores
            attn = (q_dst * k_src).sum(dim=-1) / math.sqrt(head_dim)  # (B, E, H)
            # Add edge bias
            edge_bias = self.W_e[layer](edge_feats)  # (B, E, H)
            attn = attn + edge_bias

            # Softmax per destination node per head
            attn = attn - 1e9 * torch.zeros_like(attn)  # placeholder for masking if needed
            # Sparse softmax: for each dst node, softmax over its incoming edges
            dst_idx_h = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, num_heads)
            attn_max = torch.zeros(batch, self.num_pos, num_heads, device=h.device).fill_(-1e9)
            attn_max.scatter_reduce_(1, dst_idx_h, attn, reduce='amax', include_self=False)
            attn = attn - torch.gather(attn_max, 1, dst_idx_h)
            attn_exp = torch.exp(attn)
            attn_sum = torch.zeros(batch, self.num_pos, num_heads, device=h.device)
            attn_sum.scatter_add_(1, dst_idx_h, attn_exp)
            attn_norm = attn_exp / (torch.gather(attn_sum, 1, dst_idx_h) + 1e-8)

            # Aggregate
            weighted_v = v_src * attn_norm.unsqueeze(-1)  # (B, E, H, D)
            agg = torch.zeros(batch, self.num_pos, num_heads, head_dim, device=h.device)
            agg.scatter_add_(1, dst_idx_4d, weighted_v)
            agg = agg.view(batch, self.num_pos, node_dim)

            h = h + self.update_mlps[layer](torch.cat([h, agg], dim=-1))
            h = self.layer_norms[layer](h)

        node_energy = self.node_readout(h).squeeze(-1).sum(dim=1)
        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_energy = self.edge_readout(torch.cat([h_src, h_dst, edge_feats], dim=-1)).squeeze(-1).sum(dim=1) / 2.0
        return node_energy + edge_energy


# ============================================================
# 3. DualEncoder — separate fixed/free encoders
# ============================================================

class SubtreeDualEnc(SubtreeFeatureMixin, nn.Module):
    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        edge_feat_dim = 2

        # Two separate encoders
        feat_dim_no_flag = self._feat_input_dim - 1  # without is_fixed
        self.fixed_proj = nn.Sequential(
            nn.Linear(feat_dim_no_flag, node_dim), nn.SiLU(), nn.Linear(node_dim, node_dim))
        self.free_proj = nn.Sequential(
            nn.Linear(feat_dim_no_flag, node_dim), nn.SiLU(), nn.Linear(node_dim, node_dim))

        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1))

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, _ = self._extract_features(fixed_rcs, free_mask)
        # node_input has is_fixed as last dim — strip it
        node_feat = node_input[..., :-1]

        h_fixed = self.fixed_proj(node_feat)
        h_free = self.free_proj(node_feat)
        mask_3d = free_mask.unsqueeze(-1)
        h = torch.where(mask_3d, h_free, h_fixed)

        src, dst = self.edge_src, self.edge_dst
        node_dim = self.node_dim
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)
            h_dst = torch.gather(h, 1, dst_idx)
            messages = self.message_mlps[layer](torch.cat([h_src, h_dst, edge_feats], dim=-1))
            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)
            h = h + self.update_mlps[layer](torch.cat([h, agg], dim=-1))
            h = self.layer_norms[layer](h)

        node_energy = self.node_readout(h).squeeze(-1).sum(dim=1)
        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_energy = self.edge_readout(torch.cat([h_src, h_dst, edge_feats], dim=-1)).squeeze(-1).sum(dim=1) / 2.0
        return node_energy + edge_energy


# ============================================================
# 4. MultiHead — separate readout heads by #free positions
# ============================================================

class SubtreeMultiHead(SubtreeFeatureMixin, nn.Module):
    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        edge_feat_dim = 2

        self.node_proj = nn.Sequential(nn.Linear(self._feat_input_dim, node_dim), nn.SiLU())

        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        # One readout head per possible #free (0 to num_pos-1)
        self.node_readouts = nn.ModuleList()
        self.edge_readouts = nn.ModuleList()
        for _ in range(num_pos):  # nfree = 0, 1, ..., num_pos-1
            self.node_readouts.append(nn.Sequential(
                nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1)))
            self.edge_readouts.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
                nn.Dropout(dropout), nn.Linear(hidden_dim, 1)))

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, _ = self._extract_features(fixed_rcs, free_mask)
        h = self.node_proj(node_input)

        src, dst = self.edge_src, self.edge_dst
        node_dim = self.node_dim
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)
            h_dst = torch.gather(h, 1, dst_idx)
            messages = self.message_mlps[layer](torch.cat([h_src, h_dst, edge_feats], dim=-1))
            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)
            h = h + self.update_mlps[layer](torch.cat([h, agg], dim=-1))
            h = self.layer_norms[layer](h)

        # Route each sample to its head based on n_free
        n_free = free_mask.sum(dim=1).long()  # (batch,)
        output = torch.zeros(batch, device=h.device)

        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_input = torch.cat([h_src, h_dst, edge_feats], dim=-1)

        for nf in range(self.num_pos):
            mask = n_free == nf
            if not mask.any():
                continue
            idx = mask.nonzero(as_tuple=True)[0]
            h_sub = h[idx]
            ne = self.node_readouts[nf](h_sub).squeeze(-1).sum(dim=1)
            ee = self.edge_readouts[nf](edge_input[idx]).squeeze(-1).sum(dim=1) / 2.0
            output[idx] = ne + ee

        return output


# ============================================================
# 5. Graph Transformer — full self-attention
# ============================================================

class SubtreeTransformer(SubtreeFeatureMixin, nn.Module):
    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1,
                 num_heads=4, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.num_heads = num_heads
        edge_feat_dim = 2

        self.node_proj = nn.Sequential(nn.Linear(self._feat_input_dim, node_dim), nn.SiLU())

        # Pairwise bias: project edge features to attention bias
        # Build full NxN pair feature matrix from sparse edges
        self.register_buffer("_edge_pair_idx",
            torch.stack([edge_index[0], edge_index[1]], dim=0))  # (2, E)

        self.pair_bias_proj = nn.ModuleList()
        self.transformer_layers = nn.ModuleList()
        for _ in range(num_layers):
            self.pair_bias_proj.append(nn.Linear(edge_feat_dim, num_heads))
            self.transformer_layers.append(nn.TransformerEncoderLayer(
                d_model=node_dim, nhead=num_heads, dim_feedforward=hidden_dim,
                dropout=dropout, activation='gelu', batch_first=True, norm_first=True))

        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        # Edge readout uses graph edges
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1))

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, _ = self._extract_features(fixed_rcs, free_mask)
        h = self.node_proj(node_input)  # (B, P, D)

        num_pos = self.num_pos
        src, dst = self.edge_src, self.edge_dst
        node_dim = self.node_dim

        for layer in range(self.num_layers):
            # Build attention bias from edge features
            bias = self.pair_bias_proj[layer](edge_feats)  # (B, E, H)
            # Scatter into (B, H, P, P)
            attn_bias = torch.zeros(batch, self.num_heads, num_pos, num_pos, device=h.device)
            src_exp = src.unsqueeze(0).expand(batch, -1)  # (B, E)
            dst_exp = dst.unsqueeze(0).expand(batch, -1)
            for head in range(self.num_heads):
                b_h = bias[:, :, head]  # (B, E)
                idx = src_exp * num_pos + dst_exp  # (B, E)
                attn_bias_flat = attn_bias[:, head].reshape(batch, -1).clone()  # (B, P*P)
                attn_bias_flat.scatter_add_(1, idx, b_h)
                attn_bias[:, head] = attn_bias_flat.view(batch, num_pos, num_pos)

            # Manually apply transformer with bias
            ln = self.transformer_layers[layer].norm1
            h_normed = ln(h)
            # Self attention with bias
            nhead = self.num_heads
            head_dim = node_dim // nhead
            sa = self.transformer_layers[layer].self_attn
            q = sa.in_proj_weight[:node_dim] @ h_normed.transpose(-1, -2)
            k = sa.in_proj_weight[node_dim:2*node_dim] @ h_normed.transpose(-1, -2)
            v = sa.in_proj_weight[2*node_dim:] @ h_normed.transpose(-1, -2)
            q = q.transpose(-1, -2).view(batch, num_pos, nhead, head_dim).transpose(1, 2)
            k = k.transpose(-1, -2).view(batch, num_pos, nhead, head_dim).transpose(1, 2)
            v = v.transpose(-1, -2).view(batch, num_pos, nhead, head_dim).transpose(1, 2)
            if sa.in_proj_bias is not None:
                bq, bk, bv = sa.in_proj_bias.chunk(3)
                q = q + bq.view(1, nhead, 1, head_dim)
                k = k + bk.view(1, nhead, 1, head_dim)
                v = v + bv.view(1, nhead, 1, head_dim)

            attn_w = (q @ k.transpose(-1, -2)) / math.sqrt(head_dim)
            attn_w = attn_w + attn_bias
            attn_w = F.softmax(attn_w, dim=-1)
            attn_out = (attn_w @ v).transpose(1, 2).contiguous().view(batch, num_pos, node_dim)
            attn_out = sa.out_proj(attn_out)
            h = h + attn_out

            # FFN part
            h = h + self.transformer_layers[layer].linear2(
                F.gelu(self.transformer_layers[layer].linear1(
                    self.transformer_layers[layer].norm2(h))))

        node_energy = self.node_readout(h).squeeze(-1).sum(dim=1)
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_energy = self.edge_readout(torch.cat([h_src, h_dst, edge_feats], dim=-1)).squeeze(-1).sum(dim=1) / 2.0
        return node_energy + edge_energy


# ============================================================
# 6. Transfer — load per-conf checkpoint, add free-position modules
# ============================================================

class SubtreeTransfer(SubtreeFeatureMixin, nn.Module):
    """Load pretrained per-conf GNN weights, add free-position handling, fine-tune."""

    def __init__(self, num_pos, max_rcs, num_rcs, aa_table, chi_table, pair_table,
                 edge_index, ca_dist_vec, aa_embed_dim=16, pos_embed_dim=8,
                 node_dim=32, hidden_dim=64, num_layers=3, dropout=0.1,
                 pretrained_path=None, freeze_backbone=False, **kwargs):
        nn.Module.__init__(self)
        self._init_features(num_pos, max_rcs, num_rcs, aa_table, chi_table,
                           pair_table, edge_index, ca_dist_vec, aa_embed_dim, pos_embed_dim)
        self.node_dim = node_dim
        self.num_layers = num_layers
        self.freeze_backbone = freeze_backbone
        edge_feat_dim = 2

        # Per-conf backbone (same as InteractionGNNv2, without is_fixed)
        feat_dim_no_flag = self._feat_input_dim - 1
        self.node_proj_backbone = nn.Sequential(nn.Linear(feat_dim_no_flag, node_dim), nn.SiLU())

        # Subtree adapter: projects is_fixed signal into node_dim
        self.free_adapter = nn.Sequential(
            nn.Linear(1, node_dim), nn.SiLU())

        self.message_mlps = nn.ModuleList()
        self.update_mlps = nn.ModuleList()
        self.layer_norms = nn.ModuleList()
        for _ in range(num_layers):
            self.message_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.update_mlps.append(nn.Sequential(
                nn.Linear(2 * node_dim, hidden_dim), nn.SiLU(),
                nn.Linear(hidden_dim, node_dim)))
            self.layer_norms.append(nn.LayerNorm(node_dim))

        self.node_readout = nn.Sequential(
            nn.Linear(node_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.edge_readout = nn.Sequential(
            nn.Linear(2 * node_dim + edge_feat_dim, hidden_dim), nn.SiLU(),
            nn.Dropout(dropout), nn.Linear(hidden_dim, 1))

        # Load pretrained weights
        if pretrained_path and os.path.exists(pretrained_path):
            print(f"Loading pretrained weights from {pretrained_path}")
            ckpt = torch.load(pretrained_path, map_location='cpu')
            pretrained_sd = ckpt['model_state']
            self._load_pretrained(pretrained_sd)

    def _load_pretrained(self, pretrained_sd):
        """Map InteractionGNNv2 weights to SubtreeTransfer."""
        mapping = {
            'aa_embedding.weight': 'aa_embedding.weight',
            'pos_embedding.weight': 'pos_embedding.weight',
            'node_proj.0.weight': 'node_proj_backbone.0.weight',
            'node_proj.0.bias': 'node_proj_backbone.0.bias',
        }
        # Message/update/readout layers have same names
        for layer in range(self.num_layers):
            for part in ['message_mlps', 'update_mlps']:
                for sub in ['0.weight', '0.bias', '2.weight', '2.bias']:
                    key = f'{part}.{layer}.{sub}'
                    mapping[key] = key
            mapping[f'layer_norms.{layer}.weight'] = f'layer_norms.{layer}.weight'
            mapping[f'layer_norms.{layer}.bias'] = f'layer_norms.{layer}.bias'

        for sub in ['0.weight', '0.bias', '3.weight', '3.bias']:
            mapping[f'node_readout.{sub}'] = f'node_readout.{sub}'
            mapping[f'edge_readout.{sub}'] = f'edge_readout.{sub}'

        my_sd = self.state_dict()
        loaded = 0
        for src_key, dst_key in mapping.items():
            if src_key in pretrained_sd and dst_key in my_sd:
                if pretrained_sd[src_key].shape == my_sd[dst_key].shape:
                    my_sd[dst_key] = pretrained_sd[src_key]
                    loaded += 1
                else:
                    print(f"  Shape mismatch: {src_key} {pretrained_sd[src_key].shape} vs {dst_key} {my_sd[dst_key].shape}")
            elif src_key in pretrained_sd:
                print(f"  Missing in model: {dst_key}")

        self.load_state_dict(my_sd, strict=False)
        print(f"  Loaded {loaded} pretrained parameter tensors")

        if self.freeze_backbone:
            # Freeze everything except free_adapter
            for name, param in self.named_parameters():
                if 'free_adapter' not in name:
                    param.requires_grad = False
            trainable = sum(p.numel() for p in self.parameters() if p.requires_grad)
            print(f"  Frozen backbone. Trainable params: {trainable}")

    def forward(self, fixed_rcs, free_mask):
        batch = fixed_rcs.shape[0]
        node_input, edge_feats, is_fixed = self._extract_features(fixed_rcs, free_mask)

        # Split: backbone features (without is_fixed) + adapter
        node_feat = node_input[..., :-1]
        h = self.node_proj_backbone(node_feat) + self.free_adapter(is_fixed)

        src, dst = self.edge_src, self.edge_dst
        node_dim = self.node_dim
        src_idx = src.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)
        dst_idx = dst.unsqueeze(0).unsqueeze(-1).expand(batch, -1, node_dim)

        for layer in range(self.num_layers):
            h_src = torch.gather(h, 1, src_idx)
            h_dst = torch.gather(h, 1, dst_idx)
            messages = self.message_mlps[layer](torch.cat([h_src, h_dst, edge_feats], dim=-1))
            agg = torch.zeros_like(h)
            agg.scatter_add_(1, dst_idx, messages)
            h = h + self.update_mlps[layer](torch.cat([h, agg], dim=-1))
            h = self.layer_norms[layer](h)

        node_energy = self.node_readout(h).squeeze(-1).sum(dim=1)
        h_src = torch.gather(h, 1, src_idx)
        h_dst = torch.gather(h, 1, dst_idx)
        edge_energy = self.edge_readout(torch.cat([h_src, h_dst, edge_feats], dim=-1)).squeeze(-1).sum(dim=1) / 2.0
        return node_energy + edge_energy


# ============================================================
# Architecture registry
# ============================================================

ARCH_MAP = {
    'baseline': SubtreeBaseline,
    'gat': SubtreeGAT,
    'dual_enc': SubtreeDualEnc,
    'multihead': SubtreeMultiHead,
    'transformer': SubtreeTransformer,
    'transfer': SubtreeTransfer,
}


# ============================================================
# Training
# ============================================================

def train(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    num_gpus = torch.cuda.device_count()
    print(f"Arch: {args.arch} | Device: {device}, GPUs: {num_gpus}")

    data = load_data(args.data_dir)

    fixed_rcs, free_masks, labels = generate_subtree_samples(
        data["confs"], data["emat_energies"], data["residuals"],
        data["num_pos"], min_confs=args.min_confs,
        max_samples_per_k=args.max_samples_per_k, seed=args.seed)

    n = len(labels)
    perm = np.random.RandomState(args.seed).permutation(n)
    n_val = max(1, int(n * args.val_frac))
    val_idx, train_idx = perm[:n_val], perm[n_val:]

    train_ds = SubtreeDataset(fixed_rcs[train_idx], free_masks[train_idx], labels[train_idx])
    val_ds = SubtreeDataset(fixed_rcs[val_idx], free_masks[val_idx], labels[val_idx])
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=4, pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size,
                            num_workers=4, pin_memory=True)

    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")
    print(f"Label stats: mean={labels.mean():.4f}, std={labels.std():.4f}, "
          f"min={labels.min():.4f}, max={labels.max():.4f}")

    n_free_counts = {}
    for i in range(len(free_masks)):
        nf = int(free_masks[i].sum())
        n_free_counts[nf] = n_free_counts.get(nf, 0) + 1
    print("Samples by #free:", dict(sorted(n_free_counts.items())))

    # Build model
    model_cls = ARCH_MAP[args.arch]
    model_kwargs = dict(
        num_pos=data["num_pos"], max_rcs=data["max_rcs"], num_rcs=data["num_rcs"],
        aa_table=data["aa_table"], chi_table=data["chi_table"],
        pair_table=data["pair_table"], edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        aa_embed_dim=args.aa_embed_dim, pos_embed_dim=args.pos_embed_dim,
        node_dim=args.node_dim, hidden_dim=args.hidden_dim,
        num_layers=args.num_layers, dropout=args.dropout,
    )
    if args.arch == 'transfer':
        model_kwargs['pretrained_path'] = args.pretrained_path
        model_kwargs['freeze_backbone'] = args.freeze_backbone
    if args.arch in ('gat', 'transformer'):
        model_kwargs['num_heads'] = args.num_heads

    model = model_cls(**model_kwargs).to(device)

    num_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Parameters: {num_params:,} (trainable: {trainable_params:,})")

    if num_gpus > 1:
        print(f"Using DataParallel on {num_gpus} GPUs")
        model = nn.DataParallel(model)

    optimizer = torch.optim.AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=args.lr, weight_decay=args.wd)
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
                loss = F.huber_loss(pred, target_batch, delta=args.huber_delta)
            else:
                loss = F.mse_loss(pred, target_batch)
            train_mse += F.mse_loss(pred, target_batch).item() * len(fixed_batch)

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
                val_loss += F.mse_loss(pred, target_batch).item() * len(fixed_batch)
        val_rmse = (val_loss / len(val_ds)) ** 0.5

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            raw_sd = model.module.state_dict() if isinstance(model, nn.DataParallel) else model.state_dict()
            best_state = {k: v.cpu().clone() for k, v in raw_sd.items()}

        if epoch % args.log_every == 0 or epoch == 1:
            print(f"[{args.arch}] Epoch {epoch:4d}  train={train_rmse:.4f}  "
                  f"val={val_rmse:.4f}  best={best_val_rmse:.4f}  "
                  f"lr={scheduler.get_last_lr()[0]:.2e}", flush=True)

    # Restore best
    base_model = model.module if isinstance(model, nn.DataParallel) else model
    base_model.load_state_dict(best_state)
    base_model.eval()
    model = base_model

    # Final eval by #free
    print(f"\n{'='*70}")
    print(f"[{args.arch}] Validation Results")
    print(f"{'='*70}")

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
        print(f"  #free={nf}: N={mask.sum()}, MAE={np.mean(np.abs(e)):.4f}, RMSE={np.sqrt(np.mean(e**2)):.4f}")

    # Save
    out_dir = os.path.join(args.data_dir, f"model_subtree_{args.arch}")
    os.makedirs(out_dir, exist_ok=True)
    torch.save({
        "model_state": best_state,
        "arch": args.arch,
        "num_rcs": data["num_rcs"],
        "num_pos": data["num_pos"],
        "max_rcs": data["max_rcs"],
        "args": vars(args),
        "val_rmse": best_val_rmse,
    }, os.path.join(out_dir, "subtree_checkpoint.pt"))

    # Export ONNX
    try:
        model.eval()
        dummy_rcs = torch.zeros(1, data["num_pos"], dtype=torch.long, device=device)
        dummy_mask = torch.zeros(1, data["num_pos"], dtype=torch.bool, device=device)
        onnx_path = os.path.join(out_dir, "subtree_model.onnx")
        torch.onnx.export(
            model, (dummy_rcs, dummy_mask), onnx_path,
            input_names=["fixed_rcs", "free_mask"],
            output_names=["residual"],
            dynamic_axes={"fixed_rcs": {0: "batch"}, "free_mask": {0: "batch"}, "residual": {0: "batch"}},
            opset_version=17)
        print(f"ONNX exported to {onnx_path}")
    except Exception as e:
        print(f"ONNX export failed: {e}")

    print(f"\nModel saved to {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--arch", type=str, required=True, choices=list(ARCH_MAP.keys()))
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--batch_size", type=int, default=4096)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--loss", type=str, default="huber", choices=["mse", "huber"])
    parser.add_argument("--huber_delta", type=float, default=1.0)
    parser.add_argument("--aa_embed_dim", type=int, default=16)
    parser.add_argument("--pos_embed_dim", type=int, default=8)
    parser.add_argument("--node_dim", type=int, default=32)
    parser.add_argument("--hidden_dim", type=int, default=64)
    parser.add_argument("--num_layers", type=int, default=3)
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--val_frac", type=float, default=0.15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--log_every", type=int, default=5)
    parser.add_argument("--min_confs", type=int, default=10)
    parser.add_argument("--max_samples_per_k", type=int, default=50000)
    # Transfer-specific
    parser.add_argument("--pretrained_path", type=str, default=None)
    parser.add_argument("--freeze_backbone", action="store_true")
    args = parser.parse_args()
    train(args)
