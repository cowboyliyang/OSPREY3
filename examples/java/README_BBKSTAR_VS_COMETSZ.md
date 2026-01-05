# BBK* + MARK* 示例说明

## 重要区别：CometsZ vs BBKStar

### ❌ 之前的误解（CometsZBBKStarMARKStarExample.java）

```java
// 这个例子实际上只是 CometsZ + MARK*
// 并没有真正使用 BBKStar 类！
CometsZ cometsZ = new CometsZ.Builder(objective)...
```

**问题**:
- 使用的是 `CometsZ` 类，不是 `BBKStar`
- 只配置了 minimizing energy matrix，没有 rigid matrix
- 虽然名字叫 BBKStar，但代码里根本没用到 BBKStar！

### ✅ 正确的做法（TrueBBKStarMARKStarExample.java）

```java
// 真正使用 BBKStar 类！
BBKStar bbkstar = new BBKStar(protein, ligand, complex,
                              kstarSettings, bbkstarSettings);

// BBKStar 特有：需要配置 BOTH rigid 和 minimizing
info.confSearchFactoryMinimized = ...  // minimizing 搜索
info.confSearchFactoryRigid = ...      // rigid 搜索（BBKStar特有！）

// MARK* 可以作为 partition function 方法
info.pfuncFactory = rcs -> new MARKStarBound(...);
```

## 关键区别对比

| 特性 | CometsZ + MARK* | BBK* + MARK* |
|-----|----------------|--------------|
| **主类** | `CometsZ` | `BBKStar` |
| **用途** | Multi-state thermodynamic design | Batch-based K* optimization |
| **能量矩阵** | Only minimizing | Both rigid AND minimizing |
| **批处理** | ❌ No | ✅ Yes (`numConfsPerBatch`) |
| **序列优化** | Objective window | K* score ranking |
| **MARK* 集成** | ✅ Yes | ✅ Yes |

## 运行示例

### 方式1: 运行真正的 BBK* + MARK* 例子（推荐）

```bash
cd /home/users/lz280/IdeaProjects/OSPREY3

# 使用 Gradle 任务
./gradlew runTrueBBKStarExample --console=plain

# 或使用脚本
./RUN_TRUE_BBKSTAR_EXAMPLE.sh
```

**特点**:
- ✅ 真正使用 `BBKStar` 类
- ✅ 批处理优化（8 conformations per batch）
- ✅ Rigid 和 Minimizing 能量矩阵
- ✅ MARK* 作为分区函数方法
- ⏱️ 运行时间：几分钟（比 CometsZ 快）

### 方式2: 运行 CometsZ + MARK* 例子

```bash
# 这个例子实际上是 CometsZ + MARK*，不是 BBK*
./gradlew runCometsZExample --console=plain
```

**特点**:
- ❌ 不使用 `BBKStar` 类
- ✅ 使用 `CometsZ` 进行多状态设计
- ✅ MARK* 作为分区函数方法
- ⏱️ 运行时间：10-60 分钟（慢很多）

## 代码模式对比

### 模式1: BBK* + MARK* （正确）

```java
// 来自 TestBBKStar.java

// 1. 创建 BBKStar 实例
BBKStar bbkstar = new BBKStar(protein, ligand, complex,
                              kstarSettings, bbkstarSettings);

// 2. 为每个 conf space 配置
for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

    // 配置 minimizing
    info.confEcalcMinimized = ...
    info.confSearchFactoryMinimized = ...

    // 配置 rigid (BBKStar 特有!)
    EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
        .setIsMinimizing(false)  // <-- 关键：rigid!
        .build();
    info.confSearchFactoryRigid = ...

    // MARK* 作为 pfunc
    info.pfuncFactory = rcs -> {
        MARKStarBound pfunc = new MARKStarBoundFastQueues(
            confSpace,
            rigidEmat,        // rigid matrix
            minimizingEmat,   // minimizing matrix
            ...
        );
        pfunc.setCorrections(...);  // 避免 NullPointerException
        return pfunc;
    };
}

// 3. 运行
List<KStar.ScoredSequence> results = bbkstar.run(ecalc.tasks);
```

### 模式2: CometsZ + MARK* （不同的算法）

```java
// 1. 创建 CometsZ 实例
CometsZ cometsZ = new CometsZ.Builder(objective)
    .setEpsilon(epsilon)
    .build();

// 2. 配置 states
for (CometsZ.State state : cometsZ.states) {
    state.confEcalc = ...

    // 只需要 minimizing，不需要 rigid
    state.pfuncFactory = rcs -> new MARKStarBound(...);
}

// 3. 运行
prepCometsZStates(cometsZ, ffparams, () -> {
    List<CometsZ.SequenceInfo> results =
        cometsZ.findBestSequences(numSequences);
});
```

## 参考代码

### 测试文件
- **TestBBKStar.java** - BBK* 的正确用法（包括 MARK* 集成）
  - 方法: `runBBKStar()` (lines 67-169)
  - 特别看: lines 118-129 (rigid energy setup)

- **TestCometsZWithBBKStarAndMARKStar.java** - CometsZ + MARK*
  - 注意：这个文件名有误导性！实际上是 CometsZ，不是 BBKStar

### 关键类
- `BBKStar` - BBK* 算法实现
- `CometsZ` - Multi-state design with sequence optimization
- `MARKStarBound` - MARK* partition function
- `MARKStarBoundFastQueues` - Fast queue implementation

## 总结

**如果你想要真正的 BBK* + MARK* 集成**:
✅ 使用 `TrueBBKStarMARKStarExample.java`
✅ 运行 `./gradlew runTrueBBKStarExample`

**如果你想要 CometsZ + MARK***:
✅ 使用 `CometsZBBKStarMARKStarExample.java`
✅ 运行 `./gradlew runCometsZExample`

两个都是有效的算法，但它们是**不同的**！
- **BBK***: 批处理 K* 优化
- **CometsZ**: 多状态热力学集成设计

希望这样解释清楚了! 🎯
