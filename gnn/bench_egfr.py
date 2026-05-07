"""
Benchmark EGFR GNN inference speed (simulating resistor calls).
Tests leaf (per-conf) and subtree models for protein and complex.
"""
import numpy as np
import time, os, sys, gc
import torch

sys.path.insert(0, os.path.dirname(__file__))
from train import load_data, InteractionGNNv2

BASE = "gnn_data/egfr_5pos_erlotinib"
N_BENCH = 5000

def load_leaf_model(data_dir, data):
    ckpt = torch.load(os.path.join(data_dir, "model", "gnn_checkpoint.pt"),
                      map_location="cpu", weights_only=False)
    args = ckpt["args"]
    model = InteractionGNNv2(
        num_pos=ckpt["num_pos"], max_rcs=ckpt["max_rcs"],
        aa_table=data["aa_table"], chi_table=data["chi_table"],
        pair_table=data["pair_table"], edge_index=data["edge_index"],
        ca_dist_vec=data["ca_dist_vec"],
        aa_embed_dim=args.get("aa_embed_dim", 32),
        pos_embed_dim=args.get("pos_embed_dim", 16),
        node_dim=args.get("node_dim", 64),
        hidden_dim=args["hidden_dim"],
        num_layers=args["num_layers"],
        dropout=args.get("dropout", 0.1),
    )
    model.load_state_dict(ckpt["model_state"])
    model.eval()
    return model

def bench_leaf(model, confs_np, label, device):
    N = min(N_BENCH, len(confs_np))
    ct = torch.tensor(confs_np[:N], dtype=torch.long, device=device)

    # warmup
    with torch.no_grad():
        for _ in range(20):
            model(ct[:500])

    # batch=1 latency
    n_single = min(500, N)
    with torch.no_grad():
        if device.type == "cuda": torch.cuda.synchronize()
        t0 = time.perf_counter()
        for i in range(n_single):
            model(ct[i:i+1])
        if device.type == "cuda": torch.cuda.synchronize()
        t_s = time.perf_counter() - t0
    ms1 = t_s / n_single * 1000

    # batch=4096 throughput
    with torch.no_grad():
        if device.type == "cuda": torch.cuda.synchronize()
        t0 = time.perf_counter()
        for start in range(0, N, 4096):
            model(ct[start:start+4096])
        if device.type == "cuda": torch.cuda.synchronize()
        t_b = time.perf_counter() - t0
    ms_b = t_b / N * 1000

    # full batch
    with torch.no_grad():
        if device.type == "cuda": torch.cuda.synchronize()
        t0 = time.perf_counter()
        model(ct[:N])
        if device.type == "cuda": torch.cuda.synchronize()
        t_full = time.perf_counter() - t0
    ms_full = t_full / N * 1000

    print(f"  [{label}] {N} confs, {confs_np.shape[1]} pos, {ckpt['max_rcs']} max_rcs")
    print(f"    batch=1:     {ms1:.4f} ms/conf  ({n_single} total)")
    print(f"    batch=4096:  {ms_b:.4f} ms/conf  ({t_b:.3f}s total)")
    print(f"    batch={N}:      {ms_full:.4f} ms/conf  ({t_full:.4f}s total)")
    print(f"    speedup:     {t_s/t_full:.1f}x (batch={N} vs single)")
    return {"label": label, "ms_single": ms1, "ms_batch4096": ms_b, "ms_full": ms_full}

def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device} ({torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})")
    print("=" * 60)

    results = []
    for dname in ["protein", "complex"]:
        data_dir = f"{BASE}/{dname}"
        print(f"\nLoading {dname} data...", flush=True)
        data = load_data(data_dir)
        confs = data["confs"]
        ckpt_file = os.path.join(data_dir, "model", "gnn_checkpoint.pt")

        global ckpt
        ckpt = torch.load(ckpt_file, map_location="cpu", weights_only=False)

        model = load_leaf_model(data_dir, data)
        model.to(device)
        model.edge_index = model.edge_index.to(device)
        model.ca_dist_vec = model.ca_dist_vec.to(device)
        model.pair_table = model.pair_table.to(device)
        model.aa_table = model.aa_table.to(device)
        model.chi_table = model.chi_table.to(device)

        r = bench_leaf(model, confs, f"leaf-{dname}", device)
        results.append(r)

        del model, data, confs
        gc.collect()
        if device.type == "cuda":
            torch.cuda.empty_cache()

    print("\n" + "=" * 60)
    print("Summary — Resistor inference per-conformation")
    print("-" * 60)
    print(f"{'Model':<18s} {'1x (ms)':>10s} {'4096x (ms)':>12s} {'5k (ms)':>10s}")
    print("-" * 60)
    for r in results:
        print(f"{r['label']:<18s} {r['ms_single']:10.4f} {r['ms_batch4096']:12.4f} {r['ms_full']:10.4f}")
    print("-" * 60)
    print("In resistor: ~1-10 confs per step → use batch=1 latency")
    print(f"Typical resistor step with 1 conf: {results[0]['ms_single']:.2f} ms (protein)")
    print(f"                                      {results[1]['ms_single']:.2f} ms (complex)")

if __name__ == "__main__":
    main()
