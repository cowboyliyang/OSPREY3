# BranchMARKStar residual-budget cutoff

## Goal

BranchMARKStar should not cut interactions only because they are far away or individually below a fixed threshold. A cutoff is easier to justify if it directly budgets the full energy missing from the sparse branch graph.

For a state `s` and conformation `c`, define

```text
E_full_s(c) = E_sparse_s(c) + R_s(c)
```

If the cut edges for state `s` have worst-case absolute pair energies `r_e`, then

```text
|R_s(c)| <= rho_s = sum_{e cut in s} r_e
```

for every conformation. This is the residual budget.

## Partition-function perturbation

Because Boltzmann weight is exponential,

```text
exp(-rho_s/RT) Z_sparse_s <= Z_full_s <= exp(rho_s/RT) Z_sparse_s
```

If BranchMARKStar returns sparse pfunc bounds

```text
L_s <= Z_sparse_s <= U_s
```

then a conservative full-Hamiltonian certificate is

```text
exp(-rho_s/RT) L_s <= Z_full_s <= exp(rho_s/RT) U_s
```

This bound is valid but can be loose. At 298.15 K with OSPREY's classic gas constant, `RT = 0.593050165 kcal/mol`, so each `1 kcal/mol` residual is a `10^0.732` multiplicative Z factor.

## Kstar perturbation

For

```text
K = Z_complex / (Z_protein Z_ligand)
```

the conservative sparse-to-full ratio is controlled by the sum of state residuals:

```text
K_full in [
  exp(-(rho_C + rho_P + rho_L)/RT) K_sparse_lower,
  exp( (rho_C + rho_P + rho_L)/RT) K_sparse_upper
]
```

Equivalently, to certify Kstar within a multiplicative factor `F`, the total residual budget must satisfy

```text
rho_C + rho_P + rho_L <= RT ln F
```

Useful reference points:

```text
factor 10     -> total rho <= 1.366 kcal/mol
factor 100    -> total rho <= 2.731 kcal/mol
factor 1e3    -> total rho <= 4.097 kcal/mol
factor 1e4    -> total rho <= 5.462 kcal/mol
factor 1e6    -> total rho <= 8.193 kcal/mol
factor 1e8    -> total rho <= 10.924 kcal/mol
```

## What the 11644178 numbers imply

From the observed sparse/full Z ratios:

```text
protein:  full/sparse ~= 2.59e3 to 8.09e3 -> rho_P ~= 4.66 to 5.34 kcal/mol
complex:  full/sparse ~= 1.38e3 to 4.40e3 -> rho_C ~= 4.29 to 4.97 kcal/mol
```

Ignoring ligand residual, preserving that level of sparsification corresponds to

```text
rho_C + rho_P ~= 8.95 to 10.31 kcal/mol
```

which is a Kstar uncertainty factor of roughly

```text
10^6.6 to 10^7.6
```

So the old speedup appears to live in a residual regime that is useful as a sparse objective or heuristic, but too loose to be a tight full-Hamiltonian Kstar certificate by itself.

## Better cutoff rule

The implemented default is now residual-budget cutoff:

1. Start with the complete position-pair graph.
2. Compute each edge risk as the max absolute pairwise energy over rigid and minimizing emats.
3. Sort candidate cuts from lowest risk to highest risk, breaking ties by larger template distance.
4. Cut edges while the sum of cut risks stays below `branchmarkstar.cutoff.residualBudget`.
5. Keep the graph connected when `branchmarkstar.cutoff.keepConnected=true`.

This makes the theoretical knob explicit: increasing `residualBudget` buys branch-decomposition speed, while decreasing it buys a tighter sparse-to-full Kstar certificate.

## Current caveat

The fast BranchMARKStar path still mixes sparse branch decomposition with full-emat scoring/minimized leaves. The residual theorem is the clean certificate for a consistent sparse Hamiltonian and for any pfunc bound computed against that sparse Hamiltonian. The code exposes optional full-bound inflation with `branchmarkstar.certifyFullBounds=true`, but the main practical use right now is to report the residual factor and choose a cutoff that is defensible.
