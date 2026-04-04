"""
Benchmark: batch vs single inference overhead for GNN model (PyTorch).
Measures the fundamental question: is N×single == 1×batch(N)?
"""
import numpy as np
import time
import torch
import sys, os

sys.path.insert(0, os.path.dirname(__file__))
from train import load_data, InteractionGNNv2

BASE_DIR = "gnn_data/2RL0_markstar_13pos/protein"

def load_model(base_dir, data):
    ckpt = torch.load(
        os.path.join(base_dir, "model", "gnn_checkpoint.pt"),
        map_location="cpu", weights_only=False,
    )
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
    return model

def run_bench(model, confs_all, device_tag, N=5000):
    """Run all benchmarks on given model/device."""
    confs_np = confs_all[:N]
    device = next(model.parameters()).device
    confs_t = torch.tensor(confs_np, dtype=torch.long, device=device)

    print(f"\n{'='*65}", flush=True)
    print(f"  [{device_tag}] {N} confs, {confs_np.shape[1]} positions, device={device}", flush=True)
    print(f"{'='*65}", flush=True)

    # Warmup
    with torch.no_grad():
        for _ in range(30):
            model(confs_t[:100])
        if device.type == "cuda":
            torch.cuda.synchronize()

    # ---- A: batch=1 ----
    N_single = min(2000, N)
    with torch.no_grad():
        if device.type == "cuda": torch.cuda.synchronize()
        t0 = time.perf_counter()
        for i in range(N_single):
            model(confs_t[i:i+1])
        if device.type == "cuda": torch.cuda.synchronize()
        t_single = time.perf_counter() - t0
    ms1 = t_single / N_single * 1000
    projected = t_single / N_single * N
    print(f"\n[A] batch=1 × {N_single}: {t_single:.3f}s, {ms1:.3f} ms/conf")

    # ---- B: various batch sizes ----
    for bs in [4, 16, 64, 256, 1024, 4096, N]:
        if bs > N: continue
        n_calls = (N + bs - 1) // bs
        with torch.no_grad():
            if device.type == "cuda": torch.cuda.synchronize()
            t0 = time.perf_counter()
            for start in range(0, N, bs):
                model(confs_t[start:min(start+bs, N)])
            if device.type == "cuda": torch.cuda.synchronize()
            t = time.perf_counter() - t0
        ms = t / N * 1000
        print(f"[B] batch={bs:>5d} × {n_calls:>4d}: {t:.4f}s, {ms:.4f} ms/conf, {projected/t:.1f}x vs single")

    # ---- C: single-batch scaling ----
    print(f"\n--- Single-batch scaling ---", flush=True)
    for n in [1000, 2000, 5000, 10000]:
        if n > len(confs_all): break
        ct = torch.tensor(confs_all[:n], dtype=torch.long, device=device)
        with torch.no_grad():
            if device.type == "cuda": torch.cuda.synchronize()
            t0 = time.perf_counter()
            model(ct)
            if device.type == "cuda": torch.cuda.synchronize()
            t = time.perf_counter() - t0
        print(f"[C] N={n:>7d}: {t:.4f}s, {t/n*1000:.4f} ms/conf")

    # ---- D: scattered vs continuous ----
    print(f"\n--- Scattered (10 bursts × 500, 50ms gap) vs Continuous ---")
    bs = 64
    with torch.no_grad():
        if device.type == "cuda": torch.cuda.synchronize()
        t0 = time.perf_counter()
        for start in range(0, N, bs):
            model(confs_t[start:min(start+bs, N)])
        if device.type == "cuda": torch.cuda.synchronize()
        t_cont = time.perf_counter() - t0

    t_scat = 0
    with torch.no_grad():
        for i in range(10):
            chunk = confs_t[i*500:(i+1)*500]
            if device.type == "cuda": torch.cuda.synchronize()
            t0 = time.perf_counter()
            for start in range(0, 500, bs):
                model(chunk[start:start+bs])
            if device.type == "cuda": torch.cuda.synchronize()
            t_scat += time.perf_counter() - t0
            if i < 9: time.sleep(0.05)

    print(f"  Continuous:  {t_cont:.4f}s")
    print(f"  Scattered:   {t_scat:.4f}s (compute only)")
    print(f"  Ratio: {t_scat/t_cont:.3f}x")

    # ---- E: burst batch=1 vs batch=500 ----
    print(f"\n--- Burst: batch=1 vs batch=500 ---")
    t1 = 0
    with torch.no_grad():
        for i in range(10):
            chunk = confs_t[i*500:(i+1)*500]
            if device.type == "cuda": torch.cuda.synchronize()
            t0 = time.perf_counter()
            for j in range(500):
                model(chunk[j:j+1])
            if device.type == "cuda": torch.cuda.synchronize()
            t1 += time.perf_counter() - t0
    t500 = 0
    with torch.no_grad():
        for i in range(10):
            chunk = confs_t[i*500:(i+1)*500]
            if device.type == "cuda": torch.cuda.synchronize()
            t0 = time.perf_counter()
            model(chunk)
            if device.type == "cuda": torch.cuda.synchronize()
            t500 += time.perf_counter() - t0
    print(f"  batch=1:   {t1:.4f}s ({t1/N*1000:.3f} ms/conf)")
    print(f"  batch=500: {t500:.4f}s ({t500/N*1000:.4f} ms/conf)")
    print(f"  Speedup: {t1/t500:.1f}x")


def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    print("Loading data...", flush=True)
    data = load_data(BASE_DIR)
    # Only keep first 10k confs to save memory
    confs_all = data["confs"][:10000].copy()
    # Free the big data dict
    del data["confs"], data["emat_energies"], data["ccd_energies"]
    import gc; gc.collect()
    print(f"Using {len(confs_all)} conformations for benchmark", flush=True)

    # Build model once
    model = load_model(BASE_DIR, data)

    # ==================== CPU (1 thread) ====================
    torch.set_num_threads(1)
    model.cpu()
    run_bench(model, confs_all, "CPU-1thread")

    # ==================== CPU (multi thread) ====================
    n_cpu = min(os.cpu_count() or 4, 8)
    torch.set_num_threads(n_cpu)
    print(f"\n(set torch threads to {torch.get_num_threads()})", flush=True)
    run_bench(model, confs_all, f"CPU-{n_cpu}thread")

    # ==================== GPU (if available) ====================
    if torch.cuda.is_available():
        dev = torch.device("cuda:0")
        print(f"\nGPU: {torch.cuda.get_device_name(0)}", flush=True)
        model.to(dev)
        model.edge_index = model.edge_index.to(dev)
        model.ca_dist_vec = model.ca_dist_vec.to(dev)
        model.pair_table = model.pair_table.to(dev)
        model.aa_table = model.aa_table.to(dev)
        model.chi_table = model.chi_table.to(dev)
        run_bench(model, confs_all, "GPU-A5000")
    else:
        print("\nNo GPU available, skipping GPU benchmark.")

if __name__ == "__main__":
    main()
