# PACK* H200 实验执行说明

## 固定版本与实验顺序

先固定本次提交的 Git commit，再编译和提交任务。所有编译、预检、实验和结果分析均通过 Slurm；源码可以放 home，依赖缓存、构建快照、PDB、EMAT、临时文件和结果必须放集群 scratch。以下命令中的账号、partition、GPU 类型名称和 scratch 路径需要按 H200 集群填写。

生产配置：FP64 DP、CPU CCD、CPU sampling、joint moment learning、`budget-forward` triple selection、最多 3 个 triple / 3 条 fill edge / 1,000,000 个 triple assignments。epsilon=0.683，delta=0.05，seed=42。完整参数在 `slurm/h200/production.properties`；不要用旧脚本的参数列表覆盖它。GPU CCD 已在 A5000 上试过，完整流程更慢，本轮保持 CPU CCD。

新增 triple fill edges 后，候选成本预览和实际 proposal 都重新计算 weighted Hicks 分解，并使用与初始阶段相同的 `rootSplit=gpubytes` 选根规则，带入本次运行的 host/GPU 内存预算。无需新增边时保留初始根。日志中的 `rootSplit` 是分解树的切分边编号，`rootSelectionMs` 记录选根耗时；候选预览不分配 DP 表或枚举数组。

优先顺序：

1. 编译和 GPU 正确性测试；46 个 design 的位置/序列预检；2xgy +0 完整试跑。
2. 12 个系统、46 个 flex design 的完整 P/L/PL 工作负载，对接 Duke 已提交的 MARK* 基线。
3. 原始 38 系统完整基准，避免只展示筛选出的高加速比系统。
4. 固定 8 个系统上的学习消融，以及固定系统上的 1/2 H200 扩展测试。先保留 seed=42 的主结果，再用 43/44 检查重复性。

## 12 系统、46 个 frontier design

当前正式集合为 7 个保留系统的 20 个 design，加上 5 个新增系统的 26 个 design，共 46 项。旧集合的 4z80、2rl0、2q2a 已取消，不纳入本次 H200 正式对照。`slurm/h200/frontier.tsv` 与下表对应，包含完整 mutable/flexible 位点、PDB SHA256、期望序列数、MARK* job 和逐项资源配置；H200 使用连续 index 0–45。不要沿用旧版 32 项的 array 范围，也不要重新选择加点顺序。

同系统所有 design 固定 mutable 位点及其顺序，允许标准 20 种氨基酸，枚举 WT 加至多一个同时突变，序列数为 1 + 19 × mutable 位点数；两算法都关闭 stability filter，逐序列计算 Protein、Ligand、Complex 三个状态。flexible 位点保留 WT 氨基酸身份；负 delta 从原列表末端按蛋白侧顺序移除 WT flex，正 delta 使用已经冻结的加点顺序。以下完整列表为最终实际 specs。4znc 的原列表重复出现 F731，本次只保留其第一次出现，因此原始 WT flex 为 7 个。

2xgy、4u3s 历史 MARK* 原始设计超过一周，1a0r 历史原始设计在两周超时，因此三者只用 -2、-1、0，不再向上扩展。其余系统依据历史耗时和加点扫描选择并包含每个中间层级。1a0r 的历史日志出现过 bounds 警告，不能把历史超时归因于该警告；新实验仍须核查数值状态。该 12 系统集合经过性能筛选，不能替代原始 38 系统的整体结果。

### 系统与任务索引

| 系统 | H200 array index | flex delta | design 数 | 每 design 序列数 |
|---|---|---|---:|---:|
| 2xgy | 0–2 | -2, -1, +0 | 3 | 39 |
| 4u3s | 3–5 | -2, -1, +0 | 3 | 39 |
| 3ma2 | 6–7 | +0, +1 | 2 | 39 |
| 4wyu | 8–9 | +0, +1 | 2 | 39 |
| 5a6y | 10–12 | +0, +1, +2 | 3 | 39 |
| 3k3q | 13–15 | +0, +1, +2 | 3 | 39 |
| 4wyq | 16–19 | +0, +1, +2, +3 | 4 | 20 |
| 1a0r | 20–22 | -2, -1, +0 | 3 | 20 |
| 4wwi | 23–26 | +0, +1, +2, +3 | 4 | 20 |
| 5dc4 | 27–31 | +0, +1, +2, +3, +4 | 5 | 39 |
| 4znc | 32–37 | +0, +1, +2, +3, +4, +5 | 6 | 20 |
| 5dc0 | 38–45 | +0, +1, +2, +3, +4, +5, +6, +7 | 8 | 39 |

### 全部 design 的精确 specs

下面的 mutable、flexible 均为完整列表，分号分隔；顺序必须保留。`flex` 只计 WT flexible 位点，`总点位` 包含 mutable。`m2/m1/p0/p1` 分别对应相对原始 flex 的 -2/-1/0/+1。H200 index 已重新连续编号，MARK* job 保留原 array/task 编号；应按 design_id 对齐结果，不能直接对齐两个集群的 array index。

| H200 index | design_id | mutable | flexible（完整） | flex / 总点位 | 序列数 | MARK* job |
|---:|---|---|---|---:|---:|---|
| 0 | 2xgy_flex_m2 | B189;B184 | B248;B242;B277;B192;B190;B186;B183;A83;A76 | 9 / 11 | 39 | 12631847_0 |
| 1 | 2xgy_flex_m1 | B189;B184 | B248;B242;B277;B192;B190;B186;B183;B250;A83;A76 | 10 / 12 | 39 | 12631847_1 |
| 2 | 2xgy_flex_p0 | B189;B184 | B248;B242;B277;B192;B190;B186;B183;B250;B251;A83;A76 | 11 / 13 | 39 | 12631847_2 |
| 3 | 4u3s_flex_m2 | B280;B277 | B309;B276;A32;A76;A145;A97;A94;A93;A95 | 9 / 11 | 39 | 12631847_3 |
| 4 | 4u3s_flex_m1 | B280;B277 | B309;B276;B273;A32;A76;A145;A97;A94;A93;A95 | 10 / 12 | 39 | 12631847_4 |
| 5 | 4u3s_flex_p0 | B280;B277 | B309;B276;B273;B274;A32;A76;A145;A97;A94;A93;A95 | 11 / 13 | 39 | 12631847_5 |
| 6 | 3ma2_flex_p0 | A316;A254 | B341;B343;A276;A290;A291;A271;A277;A287;A318;A274;A255;A315;A286 | 13 / 15 | 39 | 12631847_6 |
| 7 | 3ma2_flex_p1 | A316;A254 | B341;B343;A276;A290;A291;A271;A277;A287;A318;A274;A255;A315;A286;A317 | 14 / 16 | 39 | 12631847_7 |
| 8 | 4wyu_flex_p0 | A94;A70 | B212;B213;A12;A90;A96;A93;A92;A72;A69;A88 | 10 / 12 | 39 | 12631847_8 |
| 9 | 4wyu_flex_p1 | A94;A70 | B212;B213;A12;A90;A96;A93;A92;A72;A69;A88;A95 | 11 / 13 | 39 | 12631847_9 |
| 10 | 5a6y_flex_p0 | D455;D452 | D356;D359;D431;D355;D429;D450;D453;D396;C274;D363;C247;C275;C334;C331 | 14 / 16 | 39 | 12631847_10 |
| 11 | 5a6y_flex_p1 | D455;D452 | D356;D359;D431;D355;D429;D450;D453;D396;C274;D363;C247;C275;C334;C331;D454 | 15 / 17 | 39 | 12631847_11 |
| 12 | 5a6y_flex_p2 | D455;D452 | D356;D359;D431;D355;D429;D450;D453;D396;C274;D363;C247;C275;C334;C331;D454;D451 | 16 / 18 | 39 | 12631847_12 |
| 13 | 3k3q_flex_p0 | C444;C440 | B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439 | 15 / 17 | 39 | 12631847_13 |
| 14 | 3k3q_flex_p1 | C444;C440 | B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439;C453 | 16 / 18 | 39 | 12631847_14 |
| 15 | 3k3q_flex_p2 | C444;C440 | B254;B316;B253;B263;B358;B284;C416;C456;C454;C446;C436;C420;C413;C437;C439;C453;C455 | 17 / 19 | 39 | 12631847_15 |
| 16 | 4wyq_flex_p0 | E394 | D288;D289;D285;E367;E371;E373;E365;E391;E390 | 9 / 10 | 20 | 12631847_28 |
| 17 | 4wyq_flex_p1 | E394 | D288;D289;D285;E367;E371;E373;E365;E391;E390;E372 | 10 / 11 | 20 | 12631847_29 |
| 18 | 4wyq_flex_p2 | E394 | D288;D289;D285;E367;E371;E373;E365;E391;E390;E372;E370 | 11 / 12 | 20 | 12631847_30 |
| 19 | 4wyq_flex_p3 | E394 | D288;D289;D285;E367;E371;E373;E365;E391;E390;E372;E370;D284 | 12 / 13 | 20 | 12631847_31 |
| 20 | 1a0r_flex_m2 | G390 | G387;G392;B280;B281;B283;B279;B284 | 7 / 8 | 20 | 12632191_0 |
| 21 | 1a0r_flex_m1 | G390 | G387;G392;G386;B280;B281;B283;B279;B284 | 8 / 9 | 20 | 12632191_1 |
| 22 | 1a0r_flex_p0 | G390 | G387;G392;G386;G393;B280;B281;B283;B279;B284 | 9 / 10 | 20 | 12632191_2 |
| 23 | 4wwi_flex_p0 | B60 | E382;E384;E540;E378;E538;B62;B64;B65;B58 | 9 / 10 | 20 | 12632191_3 |
| 24 | 4wwi_flex_p1 | B60 | E382;E384;E540;E378;E538;B62;B64;B65;B58;E541 | 10 / 11 | 20 | 12632191_4 |
| 25 | 4wwi_flex_p2 | B60 | E382;E384;E540;E378;E538;B62;B64;B65;B58;E541;B59 | 11 / 12 | 20 | 12632191_5 |
| 26 | 4wwi_flex_p3 | B60 | E382;E384;E540;E378;E538;B62;B64;B65;B58;E541;B59;B66 | 12 / 13 | 20 | 12632191_6 |
| 27 | 5dc4_flex_p0 | A95;A28 | B176;B172;B173;B174;A25;A26;A97;A21;A46 | 9 / 11 | 39 | 12632191_7 |
| 28 | 5dc4_flex_p1 | A95;A28 | B176;B172;B173;B174;A25;A26;A97;A21;A46;A24 | 10 / 12 | 39 | 12632191_8 |
| 29 | 5dc4_flex_p2 | A95;A28 | B176;B172;B173;B174;A25;A26;A97;A21;A46;A24;A29 | 11 / 13 | 39 | 12632191_9 |
| 30 | 5dc4_flex_p3 | A95;A28 | B176;B172;B173;B174;A25;A26;A97;A21;A46;A24;A29;A47 | 12 / 14 | 39 | 12632191_10 |
| 31 | 5dc4_flex_p4 | A95;A28 | B176;B172;B173;B174;A25;A26;A97;A21;A46;A24;A29;A47;B171 | 13 / 15 | 39 | 12632191_11 |
| 32 | 4znc_flex_p0 | E534 | F716;F733;F731;F718;E499;E528;E497 | 7 / 8 | 20 | 12632191_12 |
| 33 | 4znc_flex_p1 | E534 | F716;F733;F731;F718;E499;E528;E497;E533 | 8 / 9 | 20 | 12632191_13 |
| 34 | 4znc_flex_p2 | E534 | F716;F733;F731;F718;E499;E528;E497;E533;E500 | 9 / 10 | 20 | 12632191_14 |
| 35 | 4znc_flex_p3 | E534 | F716;F733;F731;F718;E499;E528;E497;E533;E500;F732 | 10 / 11 | 20 | 12632191_15 |
| 36 | 4znc_flex_p4 | E534 | F716;F733;F731;F718;E499;E528;E497;E533;E500;F732;E535 | 11 / 12 | 20 | 12632191_16 |
| 37 | 4znc_flex_p5 | E534 | F716;F733;F731;F718;E499;E528;E497;E533;E500;F732;E535;E496 | 12 / 13 | 20 | 12632191_17 |
| 38 | 5dc0_flex_p0 | B116;B138 | B181;B137;A51;A52 | 4 / 6 | 39 | 12632191_18 |
| 39 | 5dc0_flex_p1 | B116;B138 | B181;B137;A51;A52;B136 | 5 / 7 | 39 | 12632191_19 |
| 40 | 5dc0_flex_p2 | B116;B138 | B181;B137;A51;A52;B136;A53 | 6 / 8 | 39 | 12632191_20 |
| 41 | 5dc0_flex_p3 | B116;B138 | B181;B137;A51;A52;B136;A53;B117 | 7 / 9 | 39 | 12632191_21 |
| 42 | 5dc0_flex_p4 | B116;B138 | B181;B137;A51;A52;B136;A53;B117;B115 | 8 / 10 | 39 | 12632191_22 |
| 43 | 5dc0_flex_p5 | B116;B138 | B181;B137;A51;A52;B136;A53;B117;B115;B139 | 9 / 11 | 39 | 12632191_23 |
| 44 | 5dc0_flex_p6 | B116;B138 | B181;B137;A51;A52;B136;A53;B117;B115;B139;B182 | 10 / 12 | 39 | 12632283_24 |
| 45 | 5dc0_flex_p7 | B116;B138 | B181;B137;A51;A52;B136;A53;B117;B115;B139;B182;A50 | 11 / 13 | 39 | 12632283_25 |

### PDB 文件与 SHA256

每个系统固定使用 `structures/<system>.min.reduce.renum.pdb`。残基编号来自该预处理后的 PDB，不能换成原始 PDB 的编号或重新预处理后直接沿用 specs。

| 系统 | PDB SHA256 |
|---|---|
| 2xgy | `15b6cc16fc3a32d759168c7d8a1c7f583684980f641d87ecc08687070e5de2db` |
| 4u3s | `96a204e521a8e3b9f1ee81384d85e85d1250ca75e7c00027dc6cd4bf47b7e913` |
| 3ma2 | `35574de4dd7f77096cb7dce027f6edd29efb1441a4b4155163e66e0f11756bad` |
| 4wyu | `db330b3e3968865b8225c3f7735df44d39828c3bd032bc2f2ce3d3def7c77d53` |
| 5a6y | `53066d5b532226fcdeddb6a19e6b183dedcedcdc51578da23a8d53dd19e0090a` |
| 3k3q | `c304d4f14e5f9ea391f67d824dac985117354bcf840418913063edf0c0473bf4` |
| 4wyq | `adc46edd5b0e67e235edd709ef040a657aa4b1141cf7f4573ca6c79caadbfbc0` |
| 1a0r | `1dafb36daca93ee3e6db23e7bb648d96f7f67e14fb5bce94f2efa036a3e85b2d` |
| 4wwi | `e3fe494b90ae378dd254eb049d59e71262172051f3e125553950c305ced7bd74` |
| 5dc4 | `f1dc1758bd76057ce195cabfc51efcafce37d661b9fc5cc8c24a8375649c07d1` |
| 4znc | `19d4423b65055cbb79be8cde77db8b6bf0da4e5d2ced3c0c1c809dfb710525db` |
| 5dc0 | `9c5456f069c4bef6af20afcb920c8044e420e17fe493a785ffb1be09cd0233db` |

### MARK* 基线与资源记录

全部 MARK* design 使用 FP64、CPU CCD、16 CPU、无 GPU、epsilon=0.683、关闭 stability filter、独立新建 EMAT、14 天时限。通常申请 192 GiB 内存、160 GiB Java heap；这是保守资源配置，并非经测量证明的最低需求。`5dc0_flex_p6` 和 `5dc0_flex_p7` 为缩短排队改用 128 GiB 内存、96 GiB heap，比较时单独标明这两项资源差异。

| 批次 | 正式任务 | 原始冻结输入包 | 结果目录 |
|---|---|---|---|
| 保留 20 项 | `12631847_0–15`、`12631847_28–31` | `/usr/xtmp/lz280/packstar_flex_frontier10_20260917/prep_12631844/package` | `/usr/xtmp/lz280/packstar_flex_frontier10_20260917/A12631847/T<原task>` |
| 新增中的 24 项 | `12632191_0–23` | `/usr/xtmp/lz280/packstar_flex_frontier12_20260918/prep_12632172/package` | `/usr/xtmp/lz280/packstar_flex_frontier12_20260918/A12632191/T<原task>` |
| 降内存的 2 项 | `12632283_24–25` | 同上 | `/usr/xtmp/lz280/packstar_flex_frontier12_20260918/A12632283/T<原task>` |

旧 `12631847_16–27`（4z80、2rl0、2q2a）以及旧 `12632191_24–25` 均已取消。不要把这些被替代或退出集合的结果混入正式 46 项。每个正式任务可通过 `design.json`、`run_manifest.json`、`command.sh`、`protocol.json`、`wall.time` 和 CSV 核验配置与状态。

### 传输到 H200 的输入包

已从以上两个冻结包提取正式 46 项，保留完整 PDB、已通过的位点/序列预检及校验记录，合并包为：

```
/usr/xtmp/lz280/packstar_flex_frontier12_20260918/handoff_12632294/package
```

包内有 12 个 PDB，共 6,894,938 bytes；46 个 design 各有序列清单及预检日志，另附少量配置和校验文件，总量约 7 MiB。46 项合计 1,471 个序列工作负载，每项均需计算 P/L/PL。合并审计 Slurm job `12632294` 已通过原输入 SHA256、完整位点、序列数及同系统序列顺序检查；原有冻结包保持原样。

将该目录的内容原样复制到 H200 scratch 的 `INPUT_ROOT`，即 `INPUT_ROOT/designs.tsv`、`INPUT_ROOT/structures/`、`INPUT_ROOT/preflight/` 和 `INPUT_ROOT/SHA256SUMS`。仅复制这个合并包；两份原始包中的 task index 不等于新 H200 index。传输和 `sha256sum -c SHA256SUMS` 校验通过集群允许的 Slurm 数据搬运任务执行。无需迁移历史 EMAT、Duke classpath 或整个项目运行目录。数据包不提交 Git；所有运行配置和此文档随 main 提交。

38 系统配置在 `slurm/h200/baseline38.csv`。另行转移 Duke 的 `/usr/xtmp/lz280/dance_bench/pdbs_prepped/<pdb>/<pdb>.min.reduce.renum.pdb` 中这 38 个文件，保持目录层级；设置 `INPUT_ROOT` 为其上级 `pdbs_prepped`。先在 Slurm 中统计这 38 个文件的总字节数和 SHA256，再传输并核验；不要复制整个原始数据目录。CSV 中历史 `num_seq` 不是本轮实际序列数，本轮是 WT 加至多一个同时突变：1 + 19 × mutable 位点数。

## H200 上构建

要求 JDK 17、Python 3、CUDA 12.x 的 nvcc、兼容 NVIDIA driver、Git、tar、GNU time；Maven/Gradle 依赖需要网络或预置缓存。加载集群对应 module 后设置：

```bash
export REPO=/path/to/OSPREY3
export SCRATCH=/path/to/scratch/packstar
export JAVA_HOME=/path/to/jdk-17
export NVCC=/path/to/cuda/bin/nvcc
export OSPREY_CUDA_ARCHS=90
export GRADLE_USER_HOME="$SCRATCH/gradle-home"
export BUILD_ROOT="$SCRATCH/build-release-1"
export INPUT_ROOT="$SCRATCH/frontier12-package"
export RESULT_ROOT="$SCRATCH/results"
export ACCOUNT=your_account
export GPU_PARTITION=your_h200_partition
# 日志目录必须在 sbatch 前存在；这里只创建空目录。
mkdir -p "$SCRATCH/logs"
cd "$REPO"
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:1 \
  --output="$SCRATCH/logs/build_%j.out" --error="$SCRATCH/logs/build_%j.err" \
  slurm/h200/build.slurm
```

构建脚本在 scratch 建立新的源码快照、重新编译全部 CUDA kernels 和 Java/Kotlin，执行 PACK* 测试及四条实际 GPU DP 路径的数值对照。GPU 测试被跳过也算失败。成功才写 `BUILD_ROOT/READY`。不要使用 Duke 的绝对 classpath，不要跳过 sm_90 编译；每次重新构建使用新的 `BUILD_ROOT`，预留约 100 MiB 源码快照和数 GiB 依赖/构建空间。

## 预检、试跑与正式提交

先等构建成功。所有 46 个预检成功后，再试跑 2xgy +0（index=2）。脚本从构建快照执行，输出目录唯一，拒绝覆盖旧结果；每次运行建立独立的新 EMAT。Duke 的 46 项冻结预检已通过；H200 仍须用本地新构建重新预检，以核对转移后的 PDB、位点和序列。

```bash
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" \
  --cpus-per-task=4 --mem=16G --time=00:30:00 --array=0-45 \
  --output="$SCRATCH/logs/preflight_%A_%a.out" --error="$SCRATCH/logs/preflight_%A_%a.err" \
  slurm/h200/run.slurm --mode preflight --heap-gib 8 --host-gib 4

sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:2 \
  --cpus-per-task=64 --mem=960G --time=2-00:00:00 \
  --output="$SCRATCH/logs/pilot_%j.out" --error="$SCRATCH/logs/pilot_%j.err" \
  slurm/h200/run.slurm --index 2
```

预检检查请求位点是否被忽略、完整序列数及与 MARK* 的序列顺序；试跑检查实际使用 H200 DP、CPU CCD、完整 39 条序列的 P/L/PL 状态、内存和日志中的数值异常。`manifest.json` 的 `COMPLETED` 才表示所有预期结果成功；`INCOMPLETE_ESTIMATES`、进程失败和 Slurm timeout 都单独保留。

试跑通过后提交全部 46 个（所有 design 独立参与调度，Slurm 按可用资源启动）：

```bash
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:2 \
  --cpus-per-task=64 --mem=960G --time=2-00:00:00 --array=0-45 \
  --output="$SCRATCH/logs/frontier_%A_%a.out" --error="$SCRATCH/logs/frontier_%A_%a.err" \
  slurm/h200/run.slurm
```

默认每任务 2 H200、64 CPU、850 GiB heap、800 GiB DP host budget、每卡 120 GiB GPU budget。CPU CCD 是主要耗时来源之一，不能只分配少量 CPU。若节点没有 960 GiB 可分配内存，所有正式任务统一调整，例如 `--mem=512G`，并在脚本参数中加入 `--heap-gib 440 --host-gib 400`；记录成单独资源配置，不混合汇总。不要为个别慢系统临时增加算法预算。脚本验证 heap 不超过 Slurm 内存以及 GPU 型号、可见卡数和显存预算。

每任务 mapped file 预算最多 512 GiB，另有 EMAT、审计和日志；46 路同时运行时，mapped workspace 的预算合计为 23 TiB，加上其他结果空间。实际占用需在试跑中测量；若 scratch 容量或额度不足，可以加 array 并发上限，例如 `--array=0-45%4`，4 路预算合计为 2 TiB。脚本保留中间文件供检查，归档后只清理明确属于该次运行的 `dp_mmap/`；不要修改冻结包或历史缓存。

38 系统先设置 `INPUT_ROOT="$SCRATCH/pdbs_prepped"`，用相同资源提交 `--array=0-37%4`，脚本参数加 `--cohort baseline38`。先用 `--mode preflight --heap-gib 8 --host-gib 4` 完成其预检。主实验不复用旧 EMAT。

## 8×A5000 最新六系统基线与 2×H200 对照

以下是 Duke 上 build `12626762` 的完整运行，作为 proposal root 修复前的对照保留。每个任务独占一台 GrisMan 节点，使用 104 CPU、8×RTX A5000、约 998 GiB 节点内存、850 GiB Java heap、800 GiB DP host budget、FP64 DP、CPU CCD、CPU sampling、16 个 proposal-fit 线程、seed=42、epsilon=0.683 和 `budget-forward`。工作负载是原始 flex（+0）的全部 P/L/PL 单突变序列；成功数要求一条序列的三个状态全部为 `Estimated`。计时是 `/usr/bin/time` 包住 Java 进程所得 wall，使用私有复制的历史 EMAT，因此不含重新生成 EMAT 的时间。

| 系统 | Duke job/task | 8×A5000 wall (s) | wall (min) | 三状态成功序列 |
|---|---|---:|---:|---:|
| 2xgy | 12626765_0 | 1606.76 | 26.78 | 36/39 |
| 4u3s | 12627258_1 | 2821.55 | 47.03 | 31/39 |
| 1a0r | 12627258_2 | 2758.90 | 45.98 | 13/20 |
| 3ma2 | 12627258_3 | 4358.10 | 72.64 | 36/39 |
| 4wyu | 12627258_4 | 3131.77 | 52.20 | 33/39 |
| 5a6y | 12627258_5 | 3408.17 | 56.80 | 36/39 |

六个结果相对前一 build 的成功集合不变，除时间列外格式化 CSV 完全一致；独立数值 gate `12626763` 也通过。Duke 原始审计位于 `/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A12626765` 和 `A12627258`。

2026-09-18 更新：proposal root 修复的 build `12632056` 已通过 64 项 Slurm 回归检查，包含非零根选择、预览/执行一致性、host 预算拒绝和换根后 triple 分区函数与穷举一致性。六系统重跑 array `12632058` 及汇总 job `12632064` 均已完成，退出码均为 `0:0`。每项仍使用相同的 104 CPU、8×A5000 和上述计时口径，输入校验、序列工作负载及汇总检查的资源参数均一致。H200 使用修复后的代码时，以下表作为最新 A5000 同版本对照，上表保留作旧版记录。

| 系统 | 最新 Duke job/task | 节点 | 最新 wall (s) | 最新 wall (min) | 相对旧版 wall 变化 | 三状态成功序列（旧 → 新） |
|---|---|---|---:|---:|---|---|
| 2xgy | 12632058_0 | fennario-01 | 1634.23 | 27.24 | 增加 1.7% | 36/39 → 36/39 |
| 4u3s | 12632058_1 | fennario-02 | 3057.42 | 50.96 | 增加 8.4% | 31/39 → 32/39 |
| 1a0r | 12632058_2 | fennario-03 | 2745.56 | 45.76 | 减少 0.5% | 13/20 → 13/20 |
| 3ma2 | 12632058_3 | fennario-04 | 4479.84 | 74.66 | 增加 2.8% | 36/39 → 36/39 |
| 4wyu | 12632058_4 | fennario-05 | 3028.79 | 50.48 | 减少 3.3% | 33/39 → 33/39 |
| 5a6y | 12632058_5 | fennario-06 | 3712.09 | 61.87 | 增加 8.9% | 36/39 → 36/39 |

这轮未显示普遍的端到端加速。`4u3s` 新增成功序列 `LEU MET`，没有丢失原成功序列；其他五项成功集合完全相同。因此 `4u3s` 的 wall 变化同时伴随有效结果增加，不能当作相同成功集合下的精确加速比。以上每版各一次运行，尚不能把小幅耗时差异确认为稳定收益或退化。

最新结果位于 `/usr/xtmp/lz280/packstar_proposal_root_20260917/A12632058`，其中 `six_system_comparison.json` 和 `.md` 保存新旧 wall、成功集合、输入校验及资源差异。最终重建 proposal 的日志中，4u3s、3ma2、4wyu、5a6y 分别记录 9、1、1、3 次非零 root split；2xgy、1a0r 没有该类重建记录，不能由此断言其所有 DP 根均为零。

在 H200 上必须用 `baseline38` 的同一原始 flex 配置复现这六项，不能拿扩大 flex 的 frontier design 与上表直接比较。其零基 array index 分别是 1a0r=0、2xgy=12、3ma2=20、4u3s=25、4wyu=29、5a6y=32。若 H200 节点允许，保留相同的 104 CPU、850 GiB heap 和 800 GiB host budget，只把 GPU 配置改为 2×H200：

```bash
export INPUT_ROOT="$SCRATCH/pdbs_prepped"
export RESULT_ROOT="$SCRATCH/results/a5000x8-vs-h200x2"
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:2 \
  --cpus-per-task=104 --mem=960G --time=2-00:00:00 \
  --array=0,12,20,25,29,32%6 \
  --output="$SCRATCH/logs/h200_compare_%A_%a.out" \
  --error="$SCRATCH/logs/h200_compare_%A_%a.err" \
  slurm/h200/run.slurm --cohort baseline38
```

逐系统计算 `A5000 wall / H200 wall`；大于 1 才表示该完整工作负载下 2×H200 更快。H200 runner 强制新建 EMAT，而上表 A5000 wall 使用已有 EMAT，因此还必须并列比较日志中的 search/DP 阶段时间；不能把 H200 的 EMAT 构建开销解释成 GPU 速度。若 H200 节点不能提供 104 CPU，则记录实际 CPU 数并把 wall 结果标成端到端资源配置比较，不能表述为纯 GPU 对比。还要核对序列数、三个状态的成功集合和数值结果；工作负载或成功集合不同的行不计算精确硬件加速比。

## 消融和扩展实验

固定 8 系统：1gwc、2p4a、2xxm、3k3q、2xgy、3ma2、4wyu、5a6y。使用 baseline38 中对应 index、Complex、sequence index=0，seed=42/43/44。每个 arm 用独立 JVM，通过 `--cohort baseline38 --mode pfunc --index N --seed S --arm ARM` 调用同一入口：

| ARM | proposalLearning | tripleEta | selector |
|---|---|---|---|
| no-learning | false | false | 不启用 |
| pair-only | true | false | 不启用 |
| decomposition-cost | true | true | decomposition-cost |
| budget-forward | true | true | budget-forward |

共 8 × 3 × 4 = 96 次。四个 arm 固定相同资源、统计预算和序列；在同一节点依次执行并轮换次序，每次设置相同过程超时。`no-learning` 只关闭 proposal 学习，保留同一统计校准/验证流程；`jointMomentLearning=false` 不是无学习消融。已有 Duke Titan V 消融可作机制参考，和 H200 结果分开报告。

GPU 扩展固定 2xgy、3ma2、3k3q 原始 flex，分别使用 1、2 H200；CPU 数、host budget、seed 和其他参数固定。`sbatch --gres` 与 `--gpus` 必须一致。比较完整流程和 DP 阶段的时间；不能提前假定 2 H200 的端到端时间优于 8 A5000。A5000 历史试验用了 104 CPU、旧 EMAT，需明确资源及计时差别。

## 汇总口径

保留每次 `manifest.json`、完整命令、Git commit / source checksum、硬件/Slurm 信息、`wall.time`、`run.log`、全部 CSV/TSV 和 frequency/severity 审计文件。对杀死的任务以 `sacct` 的状态/Elapsed/Timelimit 为准；遗留 `RUNNING` manifest 不代表任务还活着。

主加速比 = MARK* 完整 wall / PACK* 完整 wall；双方都使用新 EMAT且关闭 stability filter。另列 search-only 时间，不能把一方预处理计入、另一方排除。MARK* 在 14 天超时且 PACK* 完整成功时，只报告加速比 >= 14 天 / PACK* wall，不把超时记成精确运行时间。

结果按 system、flex 位点数列出 wall、搜索时间、成功序列数/总数、状态、峰值 host RSS 和加速比或下界。38 系统报告完整分母和失败原因；不能只统计成功子集而隐去其余系统。排名、区间重叠及覆盖核验仅在同一序列/状态工作负载上进行。保留 PACK* 条件频率/严重度假设标识 `conditional-relative-gauge-S0-20-not-externally-recalibrated`，不要表述成已做外部无条件校准的置信保证。
