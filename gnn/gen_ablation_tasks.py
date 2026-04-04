"""Generate task list for distributed ablation: one line per (variant, seed, label)."""
VARIANTS = [
    "baseline", "+E1", "+pairsum", "+emat_total", "+bw_loss",
    "+stratified", "+target_norm", "+huber",
    "huber+pairsum", "huber+ps+et+tn",
]
SEEDS = [42, 123, 456, 789, 1024]
LABELS = ["protein", "complex"]

tasks = []
for label in LABELS:
    for variant in VARIANTS:
        for seed in SEEDS:
            tasks.append(f"{variant}\t{seed}\t{label}")

with open("gnn/ablation_tasks.txt", "w") as f:
    for t in tasks:
        f.write(t + "\n")

print(f"Generated {len(tasks)} tasks")
