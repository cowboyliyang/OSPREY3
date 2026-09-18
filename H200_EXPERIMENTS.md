# PACK* H200 实验执行说明

## 固定版本与实验顺序

先固定本次提交的 Git commit，再编译和提交任务。所有编译、预检、实验和结果分析均通过 Slurm；源码可以放 home，依赖缓存、构建快照、PDB、EMAT、临时文件和结果必须放集群 scratch。以下命令中的账号、partition、GPU 类型名称和 scratch 路径需要按 H200 集群填写。

生产配置：FP64 DP、CPU CCD、CPU sampling、joint moment learning、`budget-forward` triple selection、最多 3 个 triple / 3 条 fill edge / 1,000,000 个 triple assignments。epsilon=0.683，delta=0.05，seed=42。完整参数在 `slurm/h200/production.properties`；不要用旧脚本的参数列表覆盖它。GPU CCD 已在 A5000 上试过，完整流程更慢，本轮保持 CPU CCD。

优先顺序：

1. 编译和 GPU 正确性测试；32 个 design 的位置/序列预检；2xgy +0 完整试跑。
2. 10 个系统、32 个 flex design 的完整 P/L/PL 工作负载，对接正在 Duke 运行的 MARK* 基线。
3. 原始 38 系统完整基准，避免只展示筛选出的高加速比系统。
4. 固定 8 个系统上的学习消融，以及固定系统上的 1/2 H200 扩展测试。先保留 seed=42 的主结果，再用 43/44 检查重复性。

## 32 个 frontier design

`slurm/h200/frontier.tsv` 是与 MARK* 已提交任务完全一致的配置，包括 mutable/flexible 位点、PDB SHA256 和期望序列数。不要重新选择加点顺序。

| 系统 | array index | 相对原始 flex | 序列数/design |
|---|---|---|---:|
| 2xgy | 0–2 | -2, -1, 0 | 39 |
| 4u3s | 3–5 | -2, -1, 0 | 39 |
| 3ma2 | 6–7 | 0, +1 | 39 |
| 4wyu | 8–9 | 0, +1 | 39 |
| 5a6y | 10–12 | 0, +1, +2 | 39 |
| 3k3q | 13–15 | 0, +1, +2 | 39 |
| 4z80 | 16–18 | 0, +1, +2 | 39 |
| 2rl0 | 19–22 | 0, +1, +2, +3 | 20 |
| 2q2a | 23–27 | 0, +1, +2, +3, +4 | 20 |
| 4wyq | 28–31 | 0, +1, +2, +3 | 20 |

2xgy、4u3s 原始 MARK* 已超过一周，只有缩小和原始设计；1a0r 原始已两周超时，未进入本组。该组依据历史 PACK* 表现筛选，是 frontier 实验，不能替代 38 系统总体结果。

MARK* array：`12631847`；每个 design 为 16 CPU、192 GiB 内存、160 GiB Java heap、14 天上限，无 GPU，关闭 stability filter。其固定输入包位于 Duke：

```
/usr/xtmp/lz280/packstar_flex_frontier10_20260917/prep_12631844/package
```

把该目录原样复制到 H200 scratch 的 `INPUT_ROOT`，包括 `designs.tsv`、`structures/`、`preflight/` 和 `SHA256SUMS`。10 个 PDB 共 6,884,974 bytes，其他是小配置和预检记录；不需要迁移历史 EMAT、Java classpath 或整个工作目录。数据包不提交 Git；包内历史说明不作为本仓库的额外文档。传输和校验通过该集群允许的 Slurm 数据搬运任务执行，在包目录运行 `sha256sum -c SHA256SUMS`。运行脚本另外检查配置一致性、PDB checksum 和预检序列顺序。

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
export INPUT_ROOT="$SCRATCH/frontier-package"
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

先等构建成功。所有 32 个预检成功后，再试跑 2xgy +0（index=2）。脚本从构建快照执行，输出目录唯一，拒绝覆盖旧结果；每次运行建立独立的新 EMAT。

```bash
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" \
  --cpus-per-task=4 --mem=16G --time=00:30:00 --array=0-31 \
  --output="$SCRATCH/logs/preflight_%A_%a.out" --error="$SCRATCH/logs/preflight_%A_%a.err" \
  slurm/h200/run.slurm --mode preflight --heap-gib 8 --host-gib 4

sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:2 \
  --cpus-per-task=64 --mem=960G --time=2-00:00:00 \
  --output="$SCRATCH/logs/pilot_%j.out" --error="$SCRATCH/logs/pilot_%j.err" \
  slurm/h200/run.slurm --index 2
```

预检检查请求位点是否被忽略、完整序列数及与 MARK* 的序列顺序；试跑检查实际使用 H200 DP、CPU CCD、完整 39 条序列的 P/L/PL 状态、内存和日志中的数值异常。`manifest.json` 的 `COMPLETED` 才表示所有预期结果成功；`INCOMPLETE_ESTIMATES`、进程失败和 Slurm timeout 都单独保留。

试跑通过后提交全部 32 个（并发 4 只是示例，可按额度放开）：

```bash
sbatch --account="$ACCOUNT" --partition="$GPU_PARTITION" --gres=gpu:h200:2 \
  --cpus-per-task=64 --mem=960G --time=2-00:00:00 --array=0-31%4 \
  --output="$SCRATCH/logs/frontier_%A_%a.out" --error="$SCRATCH/logs/frontier_%A_%a.err" \
  slurm/h200/run.slurm
```

默认每任务 2 H200、64 CPU、850 GiB heap、800 GiB DP host budget、每卡 120 GiB GPU budget。CPU CCD 是主要耗时来源之一，不能只分配少量 CPU。若节点没有 960 GiB 可分配内存，所有正式任务统一调整，例如 `--mem=512G`，并在脚本参数中加入 `--heap-gib 440 --host-gib 400`；记录成单独资源配置，不混合汇总。不要为个别慢系统临时增加算法预算。脚本验证 heap 不超过 Slurm 内存以及 GPU 型号、可见卡数和显存预算。

每任务 mapped file 预算最多 512 GiB，另有 EMAT、审计和日志；4 路并发至少考虑 2 TiB mapped workspace，加上其他结果空间。实际占用需在试跑中测量。脚本保留中间文件供检查，归档后只清理明确属于该次运行的 `dp_mmap/`；不要修改冻结包或历史缓存。

38 系统先设置 `INPUT_ROOT="$SCRATCH/pdbs_prepped"`，用相同资源提交 `--array=0-37%4`，脚本参数加 `--cohort baseline38`。先用 `--mode preflight --heap-gib 8 --host-gib 4` 完成其预检。主实验不复用旧 EMAT。

## 8×A5000 最新六系统基线与 2×H200 对照

以下是 Duke 上最终优化 build `12626762` 的最新完整运行。每个任务独占一台 GrisMan 节点，使用 104 CPU、8×RTX A5000、约 998 GiB 节点内存、850 GiB Java heap、800 GiB DP host budget、FP64 DP、CPU CCD、CPU sampling、16 个 proposal-fit 线程、seed=42、epsilon=0.683 和 `budget-forward`。工作负载是原始 flex（+0）的全部 P/L/PL 单突变序列；成功数要求一条序列的三个状态全部为 `Estimated`。计时是 `/usr/bin/time` 包住 Java 进程所得 wall，使用私有复制的历史 EMAT，因此不含重新生成 EMAT 的时间。

| 系统 | Duke job/task | 8×A5000 wall (s) | wall (min) | 三状态成功序列 |
|---|---|---:|---:|---:|
| 2xgy | 12626765_0 | 1606.76 | 26.78 | 36/39 |
| 4u3s | 12627258_1 | 2821.55 | 47.03 | 31/39 |
| 1a0r | 12627258_2 | 2758.90 | 45.98 | 13/20 |
| 3ma2 | 12627258_3 | 4358.10 | 72.64 | 36/39 |
| 4wyu | 12627258_4 | 3131.77 | 52.20 | 33/39 |
| 5a6y | 12627258_5 | 3408.17 | 56.80 | 36/39 |

六个结果相对前一 build 的成功集合不变，除时间列外格式化 CSV 完全一致；独立数值 gate `12626763` 也通过。Duke 原始审计位于 `/usr/xtmp/lz280/packstar_frontier_fullnode_20260917/A12626765` 和 `A12627258`。

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
