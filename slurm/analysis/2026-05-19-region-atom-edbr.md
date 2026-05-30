# Region-Atom EDBR: deterministic two-sided refinement for BranchMARK*

## 一句话

BranchMARK* 现在的问题不是 partition-function 求和不准，而是每个 conformation mass 的
energy sandwich 太松：

```text
E_minemat(c) <= E_true(c) <= E_rigid(c)
```

SPARSE DP table 已经把

```text
sum_c exp(-E_rigid(c)/RT)
sum_c exp(-E_minemat(c)/RT)
```

基本精确地加出来了。epsilon 仍然停在 0.97--0.99，说明继续优化 branch
decomposition / summation order 只能带来有限收益。真正需要的是：一次昂贵 oracle
call 产生一个可复用的 certificate，批量收紧一大片 conformation slice 的上下界。

新的 refinement unit 应该从 “fix 一组 residues” 升级为：

```text
separator-conditioned region atom = 带 boundary 条件的局部 partition-function bound table
```

它把 leaf-by-leaf minimization 变成 local-slice-by-local-slice 的 bound refinement。

## EDBR 框架是否仍然正确

正确。EDBR 最抽象的不变量是：

```text
C = disjoint union of blocks B in Q
L(B) <= Z(B) <= U(B)
Z^- = sum_B L(B)
Z^+ = sum_B U(B)
```

只要每一步操作满足下面任意一种：

1. partition refinement：把一个 block 拆成 disjoint children，并给每个 child
   valid bounds；
2. bound refinement：不拆 block，但把某些 block 的 L 抬高或 U 降低，同时保持
   admissibility；

那么全局始终有：

```text
Z^- <= Z <= Z^+
```

停止条件

```text
1 - Z^- / Z^+ <= epsilon
```

就是 deterministic certificate。

Branch decomposition 是一种 partition refinement。Region atom 是一种 bound
refinement。二者可以并存。

## 为什么现有 internal correction 是错的

当前 BranchMARK* internal correction 的危险逻辑是：

```java
double confCorrection = correctionMatrix.confE(node.partialConf);
node.confLowerBound = Math.min(confCorrection, node.confUpperBound);
node.recomputeZBounds(bc);
```

问题是 `correctionMatrix.confE(node.partialConf)` 只是已经固定 residues 的 corrected
g-score，不是这个 internal node 所代表的所有 completions 的 energy lower bound。
internal node 的合法 bound 应该包含未固定 residues 的 residual lower bound / subtree
contribution。

如果要对 assigned-only correction 做安全更新，形式应该是：

```text
E^-_new(node) = E^-_old(node) + delta_g
```

或者直接在 Z upper mass 上乘 factor，而不是把 full node lower energy 替换成 partial
g-score。更好的做法是 correction-aware DP / region table replacement。

## Region atom 是什么

一个 region atom 是：

```text
a = (R, B)
```

其中：

- R 是内部 region，比如 3--6 个相互作用强的 residues；
- B 是 boundary / separator，包含所有把 R 和外界连接起来的 residues；
- 对每个 boundary assignment b，atom 存一项 local bound table：

```text
b -> [L_R(b), U_R(b)]
```

这里

```text
Z_R(b) = sum_{x_R} exp(-E_R(x_R; b)/RT)
```

atom 的证书条件是：

```text
L_R(b) <= Z_R(b) <= U_R(b)    for all boundary states b
```

所以 region atom 不是一个 full conformation，也不是一个 scalar correction。它是一张
按 boundary state 索引的局部条件 partition-function 上下界表。

## Boundary 怎么处理

Boundary 不是被 region 消掉的变量。Boundary 是 table index。

pairwise graph 下最安全定义是：

```text
B = N(R) \ R
```

如果有 higher-order factor，则规则是：

```text
任何 energy factor 只要同时碰到 R 和外界，外界变量都必须放进 B。
```

这样 B 才真正分离 R 和 outside：

```text
B separates R from O = V \ (R union B)
```

region energy 只拥有：

```text
E_R(x_R; b)
  = one-body terms in R
  + pair/higher terms fully inside R
  + interaction terms crossing R and B
```

不能包含：

- B 自己的 one-body；
- B-B terms；
- outside terms。

这些留给全局或 outside factor，避免 double counting。

## Region atom 合法性的证明

固定 boundary state b。由于 B 分离 R 和 outside O，全局能量可以写成：

```text
E(x_R, b, x_O) = E_R(x_R; b) + E_O(x_O; b) + E_B(b)
```

其中每个 energy term 被唯一归属，没有 double counting。

于是 boundary slice 的 partition function 是：

```text
Z(b)
 = sum_{x_O} exp(-(E_O(x_O;b)+E_B(b))/RT)
   * sum_{x_R} exp(-E_R(x_R;b)/RT)
 = Z_O(b) * exp(-E_B(b)/RT) * Z_R(b)
```

如果 atom 给出：

```text
L_R(b) <= Z_R(b) <= U_R(b)
```

因为 outside factor 非负，所以：

```text
Z_O(b) exp(-E_B(b)/RT) L_R(b)
 <= Z(b) <=
Z_O(b) exp(-E_B(b)/RT) U_R(b)
```

这证明 region atom 可以合法嵌回全局 Z bound。

更底层地，如果 atom 给每个 local assignment 一个 energy sandwich：

```text
E^-_R(x_R;b) <= E_R(x_R;b) <= E^+_R(x_R;b)
```

则自动有：

```text
sum_{x_R} exp(-E^+_R(x_R;b)/RT)
 <= Z_R(b) <=
sum_{x_R} exp(-E^-_R(x_R;b)/RT)
```

## two-sided refinement 的正确抽象

不要强求一个 atom 同时给 additive gamma 和 eta。更安全的抽象是维护两族 certificate：

```text
E_lower(c) = max_a E_lower_a(c) <= E_true(c)
E_upper(c) = min_b E_upper_b(c) >= E_true(c)
```

于是：

```text
Z^- = sum_c exp(-E_upper(c)/RT)
Z^+ = sum_c exp(-E_lower(c)/RT)
```

这样有两个好处：

1. lower-energy certificate 用 max 合并，天然避免 additive correction double counting；
2. feasible upper-energy witness 用 min 合并，天然保持 admissibility。

只有在证明 disjoint ownership / LP packing / no double counting 后，才应该把 corrections
相加。

## Region atom 如何产生新信息

初始 energy sandwich 已经有：

```text
E_minemat(c) <= E_true(c) <= E_rigid(c)
```

但它太松。新增信息必须把 sandwich 变窄：

```text
E_lower_new(c) >= E_lower_old(c)
E_upper_new(c) <= E_upper_old(c)
```

Region atom 的新增信息来自局部 oracle：

```text
For each boundary state b:
    compute tighter [L_R(b), U_R(b)]
```

可选 oracle：

1. local exact enumeration with improved tuple/region minimization bounds；
2. local branch DP on the induced region；
3. local WMB / mini-bucket upper bound plus explicit or mean-field lower bound；
4. local feasible witness table for E_upper；
5. local lower-energy convex/relaxed certificate for E_lower。

核心不是某个 leaf 的真实能量，而是一个 boundary slice 的 reusable local certificate。

## 为什么这可能比 MARK* 好

MARK* leaf minimization 是：

```text
[l_c, u_c] -> exact w_c
```

关闭的 gap 是：

```text
Delta_leaf = u_c - l_c
```

但它只影响一个 conformation。

Region atom 覆盖一个 slice S。如果它能在该 slice 上把 lower energy 抬高 gamma，
把 upper energy 降低 eta，则 gap 收缩大致是：

```text
Delta_atom
  ~= U_S * (1 - exp(-gamma/RT))
    + L_S * (exp(eta/RT) - 1)
```

其中：

```text
U_S = sum_{c in S} exp(-E_lower(c)/RT)
L_S = sum_{c in S} exp(-E_upper(c)/RT)
```

只要：

```text
Delta_atom / cost_atom > Delta_leaf / cost_leaf
```

它就比 MARK* 的 leaf-by-leaf minimization 更划算。

当前实验里几千次 minimization 只把 epsilon 从约 0.996 降到 0.988，说明 single-leaf
refinement 的 amortization 很差。如果能找到覆盖高-mass slice 的 region atom，即使只
给出 0.2--0.5 kcal/mol 的 deterministic tightening，也可能等价于很多 leaf
minimization。

## 和 branch decomposition 的关系

Branch decomposition 不是被废掉，而是角色改变：

原来：

```text
branch decomposition = refinement unit
```

新的视角：

```text
branch/tree decomposition = boundary system + fast Z/mass query engine
region atom = refinement unit
```

Branch decomposition 可以继续用来：

1. 给 region 选择 separator boundary；
2. 查询一个 boundary slice 的 current Z mass；
3. 把 local region table 接回全局 DP；
4. 避免 region/outside double counting；
5. 控制 table size。

## Algorithm sketch

### State

维护：

```text
Q or DP factor graph
current region tables
global Z^- and Z^+
candidate region atoms
```

### Candidate generation

从当前 certificate gap 最大的地方生成 region candidates：

```text
score(R,B) = current local gap mass / estimated oracle cost
```

候选 region 可以来自：

1. high gap branch edge / separator state；
2. high gap tuple/triple/quad；
3. interaction graph 上高 beta edge bundle；
4. top queue node 的 pending subtree；
5. leaf minimization logs 中反复出现的 high residual local motif。

### Atom construction

对候选 `(R,B)`，对每个 boundary state `b`：

```text
compute or approximate local [L_R_new(b), U_R_new(b)]
verify:
    L_R_old(b) <= L_R_new(b) <= U_R_new(b) <= U_R_old(b)
or at least:
    L_R_new(b) <= true Z_R(b) <= U_R_new(b)
```

如果 table 太大，则拒绝该 candidate。

### Bound update

把旧 local table 替换为更紧 table：

```text
L_R(b) <- max(L_R_old(b), L_R_new(b))
U_R(b) <- min(U_R_old(b), U_R_new(b))
```

然后用 DP/WMB 重新传播到 root，得到 global：

```text
Z^-_new, Z^+_new
```

### Selection rule

选择单位成本 gap 收缩最大的 atom：

```text
rho(atom) =
    [(Z^+_old - Z^+_new) + (Z^-_new - Z^-_old)] / cost(atom)
```

更接近 epsilon 的目标函数可以用：

```text
rho_eps(atom) =
    [epsilon_old - epsilon_new] / cost(atom)
```

## 最小可行实验

不要一上来实现完整 two-sided system。先做一个 diagnostic / proof-of-concept：

1. 从当前 sparse DP table 中找 top-k high-gap edge/subtree/boundary states。
2. 选很小的 region：

```text
R = 3 or 4 residues
B = graph neighbors or branch separator
```

要求：

```text
|B| small enough that boundary table fits
```

3. 先只做 table-level what-if：

```text
假设某 region slice 的 U 降低 0.1/0.3/0.5 kcal/mol
假设某 region slice 的 L 提高 0.1/0.3/0.5 kcal/mol
```

看全局 epsilon 能下降多少，换算成多少 leaf minimizations。

4. 如果 what-if 显示 region atom amortization 足够大，再做真实 local oracle。

## 主要风险

### Risk 1: upper-side feasible witness 难复用

降低 E_upper、提高 Z^- 需要可复用 feasible witness。如果 witness 依赖外部
continuous geometry，而 boundary 只存 RC assignment，则不一定合法。

解决：

1. boundary 必须包含所有影响 local feasible energy 的状态；
2. 或者 local upper witness 对 boundary continuous DOFs 取 worst-case envelope；
3. 或者先只做 lower-side atom，再配合 selected leaf minimization 修 lower Z。

### Risk 2: boundary table explosion

atom table size 是：

```text
prod_{i in B} q_i
```

如果 boundary 太大，这条路不可行。

解决：

1. 只选小 separator；
2. 用 branch/tree decomposition 的 M-set 当 boundary；
3. 对大 boundary 用 WMB/mini-bucket，而不是 exact table；
4. 拒绝 cost/gain 不划算的 atom。

### Risk 3: double counting

region factor replacement 必须有 energy-term ownership。每个 one-body、pairwise、
higher-order term 只能归属一次。

解决：

1. 明确 owner rule；
2. region atom 只替换自己 owned terms 的 local contribution；
3. boundary terms 不进 region table；
4. correction 默认 max/min 合并，只有证明后才 additive。

### Risk 4: atom 太保守

如果 local bound 没比 emat table 紧多少，atom 不如 leaf minimization。

解决：

先做 what-if 和 offline profiling，只有 high expected gain 的 region 才计算。

## Theory statement draft

### Definition: region atom

Let `G=(V,E)` be the interaction graph. A region atom is a pair `(R,B)` such that
`B` separates `R` from `O=V\(R union B)`. For each boundary assignment `b`, the
atom stores `[L_R(b), U_R(b)]`.

### Assumption: local soundness

For all boundary states `b`,

```text
L_R(b) <= Z_R(b) <= U_R(b)
```

where `Z_R(b)` is the local conditional partition function over region variables
using only the energy terms owned by `R` and interactions between `R` and `B`.

### Theorem: sound replacement

Replacing the previous local region table by a region atom preserving local
soundness preserves the global EDBR invariant:

```text
Z^- <= Z <= Z^+
```

### Proof sketch

Boundary separation implies factorization of each boundary slice:

```text
Z(b) = outside_factor(b) * Z_R(b)
```

Since outside factor is nonnegative, multiplying the local inequality by the
outside factor preserves the inequality. Summing over all boundary states yields
the global sandwich.

## Overall assessment

This idea is high-risk, high-upside.

It is better aligned with the actual bottleneck than replacing branch
decomposition by another graph decomposition, because it attacks the energy-bound
gap rather than the summation schedule.

It is provable if implemented as separator-conditioned factor replacement with
strict energy ownership and local soundness.

It can beat MARK* only if region atoms have good amortization:

```text
gap closed per oracle second > leaf minimization gap closed per oracle second
```

The first concrete milestone should not be a full implementation. It should be a
diagnostic that estimates, from current DP mass tables, how much epsilon would
drop if selected high-mass regions received 0.1/0.3/0.5 kcal/mol two-sided
tightening. If that number is large, the region-atom route is worth building.

