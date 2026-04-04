"""Per-sequence partition function error analysis."""
import pandas as pd
import numpy as np
import torch, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from train import InteractionGNNv2, load_data

kT = 0.5924
base_dir = sys.argv[1] if len(sys.argv) > 1 else 'gnn_data/2RL0_all20_4pos'

for label in ['protein', 'complex']:
    base = f'{base_dir}/{label}'
    data = load_data(base)
    ckpt = torch.load(f'{base}/model/gnn_checkpoint.pt', map_location='cpu')
    a = ckpt['args']
    
    model = InteractionGNNv2(
        data['num_pos'], data['max_rcs'],
        data['aa_table'], data['chi_table'], data['pair_table'],
        data['edge_index'], data['ca_dist_vec'],
        aa_embed_dim=a['aa_embed_dim'], pos_embed_dim=a['pos_embed_dim'],
        node_dim=a['node_dim'], hidden_dim=a['hidden_dim'],
        num_layers=a['num_layers'], dropout=a['dropout']
    )
    model.load_state_dict(ckpt['model_state'])
    model.eval()
    
    confs = data['confs']
    E_emat = data['emat_energies']
    E_CCD = data['ccd_energies']
    
    # Batch inference to save memory
    E_GNN = np.zeros(len(confs))
    bs = 50000
    for s in range(0, len(confs), bs):
        e = min(s + bs, len(confs))
        with torch.no_grad():
            pred = model(torch.tensor(confs[s:e], dtype=torch.long)).squeeze().numpy()
        E_GNN[s:e] = E_emat[s:e] + pred
    
    # Map rc -> aa
    rc_feat_df = pd.read_csv(f'{base}/rc_features.csv')
    rc_to_aa = {}
    for pos in range(data['num_pos']):
        pos_rcs = rc_feat_df[rc_feat_df.pos == pos]
        rc_to_aa[pos] = dict(zip(pos_rcs.rc, pos_rcs.aa_type_idx))
    
    # Group by sequence
    from collections import defaultdict
    seq_groups = defaultdict(list)
    for i in range(len(confs)):
        sk = tuple(rc_to_aa[p][confs[i,p]] for p in range(data['num_pos']))
        seq_groups[sk].append(i)
    
    results = []
    for seq, idxs in seq_groups.items():
        idxs = np.array(idxs)
        pf_true = np.exp(-E_CCD[idxs] / kT).sum()
        pf_gnn = np.exp(-E_GNN[idxs] / kT).sum()
        if pf_true > 0 and pf_gnn > 0:
            results.append({
                'seq': seq, 'n': len(idxs),
                'log10_pf': np.log10(pf_true),
                'log_ratio': np.log10(pf_gnn / pf_true),
                'rel_err': abs(pf_gnn - pf_true) / pf_true,
                'mae': abs(E_GNN[idxs] - E_CCD[idxs]).mean(),
                'max_err': abs(E_GNN[idxs] - E_CCD[idxs]).max()
            })
    
    R = pd.DataFrame(results)
    print(f'\n{"="*60}')
    print(f' {label.upper()} — Per-seq pfunc ({len(R)} seqs, all n)')
    print(f' (from {len(seq_groups)} total sequences)')
    print(f'{"="*60}')
    
    # Overall stats
    print(f'Confs per seq: min={R.n.min()} max={R.n.max()} median={R.n.median():.0f} mean={R.n.mean():.1f}')
    
    print(f'\nlog10(pf_GNN / pf_true):')
    print(f'  mean={R.log_ratio.mean():.4f}  std={R.log_ratio.std():.4f}  median={R.log_ratio.median():.4f}')
    print(f'  min={R.log_ratio.min():.4f}  max={R.log_ratio.max():.4f}')
    for t in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0]:
        n = (R.log_ratio.abs() > t).sum()
        print(f'  |Dlog10| > {t}: {n}/{len(R)} ({n/len(R)*100:.1f}%)')
    
    print(f'\nRelative error |pf_GNN-pf_true|/pf_true:')
    print(f'  mean={R.rel_err.mean()*100:.2f}%  median={R.rel_err.median()*100:.2f}%')
    for p in [0.9, 0.95, 0.99, 1.0]:
        v = R.rel_err.quantile(p) if p < 1.0 else R.rel_err.max()
        lbl = f'P{int(p*100)}' if p < 1.0 else 'max'
        print(f'  {lbl}={v*100:.2f}%')
    
    print(f'\nWorst 10:')
    for _, r in R.nlargest(10, 'rel_err').iterrows():
        print(f'  aa={r.seq} n={r.n:3d} log10pf={r.log10_pf:6.1f} Dlog10={r.log_ratio:+.4f} rel={r.rel_err*100:6.2f}% mae={r.mae:.4f} max_e={r.max_err:.3f}')
    
    # By n_confs buckets
    print(f'\nError by #confs per seq:')
    for lo, hi in [(1,1),(2,5),(6,20),(21,100),(101,10000)]:
        sub = R[(R.n >= lo) & (R.n <= hi)]
        if len(sub) > 0:
            print(f'  n={lo}-{hi}: {len(sub)} seqs, med_rel={sub.rel_err.median()*100:.2f}% max_rel={sub.rel_err.max()*100:.1f}% med_|Dlog10|={sub.log_ratio.abs().median():.4f}')
    
    # By pfunc quartile
    print(f'\nError by pfunc quartile:')
    R_s = R.sort_values('log10_pf')
    for i, q in enumerate(np.array_split(R_s, 4)):
        print(f'  Q{i+1} (log10pf {q.log10_pf.min():.1f}~{q.log10_pf.max():.1f}): med_rel={q.rel_err.median()*100:.2f}% max={q.rel_err.max()*100:.1f}% med_|Dlog10|={q.log_ratio.abs().median():.4f}')
    
    del model, data
    import gc; gc.collect()
