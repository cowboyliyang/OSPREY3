#!/usr/bin/env python3
"""Re-export subtree ONNX models from saved checkpoints (no retraining needed)."""
import sys, os, argparse, torch

sys.path.insert(0, os.path.dirname(__file__))
from train_subtree import SubtreeGNN, export_onnx

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", help="Path to subtree_checkpoint.pt")
    parser.add_argument("--out_dir", help="Output directory (default: same as checkpoint)")
    args = parser.parse_args()

    ckpt = torch.load(args.checkpoint, map_location="cpu")
    hparams = ckpt["args"]
    num_pos = ckpt["num_pos"]
    max_rcs = ckpt["max_rcs"]
    num_rcs = ckpt["num_rcs"]
    state = ckpt["model_state"]

    # Extract buffer tensors from state_dict to construct model
    aa_table = state["aa_table"]
    chi_table = state["chi_table"]
    pair_table = state["pair_table"]
    ca_dist_vec = state["ca_dist_vec"]
    edge_index = torch.stack([state["edge_src"], state["edge_dst"]], dim=0)

    model = SubtreeGNN(
        num_pos=num_pos, max_rcs=max_rcs, num_rcs=num_rcs,
        aa_table=aa_table, chi_table=chi_table, pair_table=pair_table,
        edge_index=edge_index, ca_dist_vec=ca_dist_vec,
        aa_embed_dim=hparams["aa_embed_dim"], pos_embed_dim=hparams["pos_embed_dim"],
        node_dim=hparams["node_dim"], hidden_dim=hparams["hidden_dim"],
        num_layers=hparams["num_layers"], dropout=hparams["dropout"],
    )
    model.load_state_dict(state)

    out_dir = args.out_dir or os.path.dirname(args.checkpoint)
    export_onnx(model, num_pos, out_dir, device=torch.device("cpu"))
    print(f"Done. val_rmse was {ckpt['val_rmse']:.6f}")

if __name__ == "__main__":
    main()
