"""Reload checkpoint + val set, compute CP quantiles, save back."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ""))
sys.path.insert(0, "/home/users/lz280/IdeaProjects/OSPREY3/gnn")

import numpy as np
import torch
from train import load_data, InteractionGNNv2

def calibrate(base_dir):
    print(f"\n{'='*60}")
    print(f"  {base_dir}")
    print(f"{'='*60}")

    data = load_data(base_dir)
    ckpt_path = os.path.join(base_dir, "model", "gnn_checkpoint.pt")
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    args_dict = ckpt["args"]

    model = InteractionGNNv2(
        num_pos=ckpt["num_pos"], max_rcs=ckpt["max_rcs"],
        aa_table=data["aa_table"], chi_table=data["chi_table"],
        pair_table=data["pair_table"], edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        aa_embed_dim=args_dict.get("aa_embed_dim", 16),
        pos_embed_dim=args_dict.get("pos_embed_dim", 8),
        node_dim=args_dict.get("node_dim", 64),
        hidden_dim=args_dict["hidden_dim"],
        num_layers=args_dict["num_layers"],
        dropout=args_dict.get("dropout", 0.1),
    )
    model.load_state_dict(ckpt["model_state"])
    model.eval()

    confs = data["confs"]
    residuals = data["residuals"]
    emat = data["emat_energies"]
    N = len(confs)

    seed = args_dict.get("seed", 42)
    val_frac = args_dict.get("val_frac", 0.15)
    perm = np.random.RandomState(seed).permutation(N)
    n_val = max(1, int(N * val_frac))
    val_idx = perm[:n_val]

    print(f"  N={N}, n_val={n_val}, seed={seed}")

    # Inference on val set
    confs_t = torch.tensor(confs[val_idx], dtype=torch.long)
    with torch.no_grad():
        preds = []
        for i in range(0, len(confs_t), 4096):
            preds.append(model(confs_t[i:i+4096]).numpy())
        preds = np.concatenate(preds)

    targets = residuals[val_idx]
    val_errors = np.abs(preds - targets)
    n_cal = len(val_errors)

    print(f"\n  Val set error distribution (|pred_residual - true_residual|):")
    for p in [50, 90, 95, 99, 99.5, 99.9]:
        print(f"    P{p}: {np.percentile(val_errors, p):.6f}")
    print(f"    Max:  {np.max(val_errors):.6f}")
    print(f"    MAE:  {np.mean(val_errors):.6f}")
    print(f"    RMSE: {np.sqrt(np.mean(val_errors**2)):.6f}")

    cp_stats = {}
    print(f"\n  CP quantiles (finite-sample corrected):")
    for alpha in [0.20, 0.10, 0.05, 0.02, 0.01, 0.005, 0.001]:
        level = min(np.ceil((1 - alpha) * (n_cal + 1)) / n_cal, 1.0)
        q = float(np.quantile(val_errors, level))
        cp_stats[f"q_alpha_{alpha}"] = q
        print(f"    α={alpha:.3f}  q={q:.6f} kcal/mol  coverage≥{1-alpha:.1%}")

    for p in [50, 90, 95, 99, 99.5, 99.9]:
        cp_stats[f"P{p}"] = float(np.percentile(val_errors, p))
    cp_stats["max"] = float(np.max(val_errors))
    cp_stats["n_cal"] = n_cal

    # Save back
    ckpt["cp_stats"] = cp_stats
    torch.save(ckpt, ckpt_path)
    print(f"\n  Saved cp_stats to {ckpt_path}")

os.chdir("/home/users/lz280/IdeaProjects/OSPREY3")
for d in [
    "gnn_data/2RL0_all20_4pos_merged/protein",
    "gnn_data/2RL0_all20_4pos_merged/complex",
]:
    if os.path.exists(os.path.join(d, "model", "gnn_checkpoint.pt")):
        calibrate(d)
