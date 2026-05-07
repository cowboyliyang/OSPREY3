#!/usr/bin/env python3
"""Export all ablation subtree checkpoints to ONNX."""
import os
import sys
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from train_subtree_ablation import ARCH_MAP
from train_subtree import load_data


def export_one(ckpt_path, data_dir, device="cpu"):
    ckpt = torch.load(ckpt_path, map_location=device)
    arch = ckpt["arch"]
    saved_args = ckpt["args"]
    num_pos = ckpt["num_pos"]

    data = load_data(data_dir)

    model_cls = ARCH_MAP[arch]
    model_kwargs = dict(
        num_pos=data["num_pos"], max_rcs=data["max_rcs"], num_rcs=data["num_rcs"],
        aa_table=data["aa_table"], chi_table=data["chi_table"],
        pair_table=data["pair_table"], edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        aa_embed_dim=saved_args.get("aa_embed_dim", 16),
        pos_embed_dim=saved_args.get("pos_embed_dim", 8),
        node_dim=saved_args.get("node_dim", 32),
        hidden_dim=saved_args.get("hidden_dim", 64),
        num_layers=saved_args.get("num_layers", 3),
        dropout=saved_args.get("dropout", 0.1),
    )
    if arch == 'transfer':
        model_kwargs['pretrained_path'] = None
        model_kwargs['freeze_backbone'] = False
    if arch in ('gat', 'transformer'):
        model_kwargs['num_heads'] = saved_args.get("num_heads", 4)

    model = model_cls(**model_kwargs).to(device)
    model.load_state_dict(ckpt["model_state"])
    model.eval()

    dummy_rcs = torch.zeros(1, num_pos, dtype=torch.long, device=device)
    dummy_mask = torch.zeros(1, num_pos, dtype=torch.bool, device=device)

    onnx_path = os.path.join(os.path.dirname(ckpt_path), "subtree_model.onnx")
    torch.onnx.export(
        model, (dummy_rcs, dummy_mask), onnx_path,
        input_names=["fixed_rcs", "free_mask"],
        output_names=["residual"],
        dynamic_axes={
            "fixed_rcs": {0: "batch"},
            "free_mask": {0: "batch"},
            "residual": {0: "batch"},
        },
        opset_version=17,
    )
    print(f"  Exported: {onnx_path}")


def main():
    base = "gnn_data/2RL0_all20_4pos_merged"
    # GAT uses scatter_reduce which is not supported by ONNX export; skip it
    archs = ["baseline", "transformer", "dual_enc", "multihead", "transfer"]

    for cs in ["complex", "protein"]:
        data_dir = os.path.join(base, cs)

        for arch in archs:
            model_dir = os.path.join(data_dir, f"model_subtree_{arch}")
            ckpt_file = os.path.join(model_dir, "subtree_checkpoint.pt")
            onnx_file = os.path.join(model_dir, "subtree_model.onnx")
            if not os.path.exists(ckpt_file):
                print(f"[{cs}/{arch}] No checkpoint, skipping")
                continue
            if os.path.exists(onnx_file):
                print(f"[{cs}/{arch}] ONNX exists, skipping")
                continue
            print(f"[{cs}/{arch}] Exporting...")
            export_one(ckpt_file, data_dir)

    print("\nDone!")


if __name__ == "__main__":
    main()
