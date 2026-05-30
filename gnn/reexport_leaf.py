#!/usr/bin/env python3
"""Re-export leaf GNN ONNX models from saved checkpoints (no retraining)."""
import sys, os, argparse, torch

sys.path.insert(0, os.path.dirname(__file__))
from train import InteractionGNNv2, export_onnx

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", help="Path to gnn_checkpoint.pt")
    parser.add_argument("--out_dir", help="Output directory (default: same as checkpoint)")
    args = parser.parse_args()

    ckpt = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    h = ckpt["args"]
    state = ckpt["model_state"]

    aa_table = state["aa_table"]
    chi_table = state["chi_table"]
    pair_table = state["pair_table"]
    ca_dist_vec = state["ca_dist_vec"]
    edge_index = torch.stack([state["edge_src"], state["edge_dst"]], dim=0)

    model = InteractionGNNv2(
        num_pos=ckpt["num_pos"], max_rcs=ckpt["max_rcs"],
        aa_table=aa_table, chi_table=chi_table, pair_table=pair_table,
        edge_index=edge_index, ca_dist_vec=ca_dist_vec,
        aa_embed_dim=h["aa_embed_dim"], pos_embed_dim=h["pos_embed_dim"],
        node_dim=h["node_dim"], hidden_dim=h["hidden_dim"],
        num_layers=h["num_layers"], dropout=h["dropout"],
    )
    # Older checkpoints have aa_embedding of size [20, D]; current model expects
    # [21, D] (extra row for unknown AAs). Pad missing rows with zeros so the
    # weights match the current architecture.
    aa_w = state.get("aa_embedding.weight")
    expected = model.aa_embedding.weight.shape[0]
    if aa_w is not None and aa_w.shape[0] < expected:
        pad = torch.zeros(expected - aa_w.shape[0], aa_w.shape[1], dtype=aa_w.dtype)
        state["aa_embedding.weight"] = torch.cat([aa_w, pad], dim=0)
    model.load_state_dict(state)

    out_dir = args.out_dir or os.path.dirname(args.checkpoint)
    export_onnx(model, ckpt["num_pos"], out_dir, device=torch.device("cpu"))
    print(f"Done. val_rmse was {ckpt['val_rmse']:.6f}")

if __name__ == "__main__":
    main()
