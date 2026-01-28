# Triple Energy Context Dependency在MARK*中的表现

## 关键发现

同一个triple `[pos1=RC28, pos2=RC29, pos4=RC8]`在不同的partial conformation中会产生**完全不同的minimized energy**。

### 实验证据

从`TEST-edu.duke.cs.osprey.markstar.TestDPvsOriginal.xml`：

**Case A (Line 148):**
```
[MINIMIZE] conf=32    28    29    32    8     48    28    3     (8 positions)
[TRIPLE_CORRECTION] tuple=Res 1 RC 28 Res 2 RC 29 Res 4 RC 8
    pairwise=-13.964511
    minimized=-13.391970
    correction=0.572541
```

**Case B (Line 2587):**
```
[MINIMIZE] conf=32    28    29    9     8     (5 positions)
[TRIPLE_CORRECTION] tuple=Res 1 RC 28 Res 2 RC 29 Res 4 RC 8
    pairwise=-1.944728
    minimized=-1.813970
    correction=0.130758
```

**差异：**
- Pairwise: -13.964511 vs -1.944728 = **12.02 kcal/mol**
- Minimized: -13.391970 vs -1.813970 = **11.58 kcal/mol**
- 两者都有巨大差异！

---

## 为什么会这样？

### 1. Triple Correction在A* Tree的Partial Conformations中计算

MARK*不是在完整的8-residue conformations中计算triple corrections，而是在**A* tree traversal的中间节点**计算。

#### Case A Context
```
Full conf: [pos0=32, pos1=28, pos2=29, pos3=32, pos4=8, pos5=48, pos6=28, pos7=3]
                                       ^^^^^
Triple在这个**8-position partial assignment**的背景下被计算
```

#### Case B Context
```
Partial conf: [pos0=32, pos1=28, pos2=29, pos3=9, pos4=8]
                                           ^^^^
Triple在这个**5-position partial assignment**的背景下被计算
```

注意pos3的不同：
- Case A: pos3=32
- Case B: pos3=9

### 2. makeMolecule()的行为

当`minimizingEcalc.calcEnergy(tuple)`被调用时：

```java
// In ConfEnergyCalculator.java:285
public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag) {
    return calcEnergy(frag, makeFragInters(frag));
}

// Creates molecule
ParametricMolecule pmol = confSpace.makeMolecule(frag);
```

**关键：`makeMolecule()`只接收triple的RCTuple，不知道calling context！**

```java
// In SimpleConfSpace.java:669-700
public Molecule makeDiscreteMolecule(RCTuple conf) {
    Molecule mol = new Molecule();
    for (Residue res : molTemplate.residues) {
        Position pos = getPositionOrNull(res.getPDBResNumber());
        if (pos != null) {
            int index = conf.pos.indexOf(pos.index);
            if (index >= 0) {
                // Use specified RC
                ResidueConf rc = pos.resConfs.get(conf.RCs.get(index));
                Residue newRes = res.copyToMol(mol, false);
                rc.updateResidue(pos.strand.templateLib, newRes, mutAlignmentCache);
                continue;
            }
        }
        // ⚠️ Otherwise, copy from TEMPLATE molecule
        res.copyToMol(mol, true);
    }
    return mol;
}
```

### 3. 问题所在：Unassigned Positions使用Template Coordinates

对于triple `[pos1=28, pos2=29, pos4=8]`：

**Molecule construction:**
```
Assigned (from RCTuple):
  pos1 → RC 28
  pos2 → RC 29
  pos4 → RC 8

Unassigned (从TEMPLATE复制):
  pos0 → Template coordinates (NOT context RC 32!)
  pos3 → Template coordinates (NOT context RC 32 or 9!)
  pos5 → Template coordinates (NOT context RC 48!)
  pos6 → Template coordinates (NOT context RC 28!)
  pos7 → Template coordinates (NOT context RC 3!)
```

**这导致什么？**

Triple energy的计算**不仅**取决于三个minimized residues，还取决于：
1. **Shell residues** (via `makeFragInters()` → `EnergyPartition.makeFragment()` → `addShell()`)
2. **Shell residues的坐标** = 所有unassigned positions的template coordinates

### 4. 不同Context → 不同的Unassigned Positions → 不同的Shell

**Case A: 8-position conf**
- Assigned: 8 positions
- Unassigned (shell): 其他所有positions，使用template coords

**Case B: 5-position conf**
- Assigned: 5 positions
- Unassigned (shell): 其他所有positions，**包括pos5, pos6, pos7**，使用template coords

**等等，template coords应该是一样的啊！为什么energy不同？**

---

## 真正的原因：Energy Calculation的ResidueInteractions

让我们看`makeFragInters()`到底计算什么：

```java
// In ConfEnergyCalculator.java:218-224
public ResidueInteractions makeFragInters(RCTuple frag) {
    if (addShellInters) {
        return EnergyPartition.makeAll(confSpace, eref, addResEntropy, frag);
    } else {
        return EnergyPartition.makeFragment(confSpace, eref, addResEntropy, frag);
    }
}
```

### makeAll()的内容

```java
// In EnergyPartition.java:324-328
public static ResidueInteractions makeAll(..., RCTuple frag) {
    var inters = makeFragment(confSpace, eref, addResEntropy, frag);
    inters.addAll(makeShell(confSpace));
    return inters;
}

public static ResidueInteractions makeShell(SimpleConfSpace confSpace) {
    return ResInterGen.of(confSpace)
        .addShellIntras()    // ⚠️ Shell internal energies
        .addShellInters()    // ⚠️ Shell-shell pairwise energies
        .make();
}
```

### 关键点：Energy包括哪些Terms？

对于triple `[pos1, pos2, pos4]`，energy calculation包括：

```
E = Intra(pos1) + Intra(pos2) + Intra(pos4)                    // Fragment intras
  + Inter(pos1,pos2) + Inter(pos1,pos4) + Inter(pos2,pos4)    // Fragment inters
  + Shell(pos1) + Shell(pos2) + Shell(pos4)                    // Fragment-shell
  + Intra(shell) + Inter(shell,shell)                          // Shell self-energy
```

**Shell self-energy！这就是关键！**

Shell residues是所有**非flexible**的residues，这包括：
1. 真正的静态shell residues (e.g., 水分子, 配体等)
2. **Unassigned positions in the partial conformation!**

### 在MARK*中，Shell的定义是动态的！

```java
// In SimpleConfSpace.java:437-455
public Set<String> getShellResNumbers() {
    Set<String> shellResNumbers = new LinkedHashSet<>();
    // Add residues within shellDist of flexible positions
    for (Position pos : positions) {
        for (Residue res : pos.strand.mol.residues) {
            double dist = pos.calcDistTo(res);
            if (dist <= shellDist && !isFlexible(res)) {
                shellResNumbers.add(res.getPDBResNumber());
            }
        }
    }
    return shellResNumbers;
}
```

**"非flexible"在partial conformation中的含义：**
- 在8-position conf中，8个positions是flexible
- 在5-position conf中，**只有5个positions是flexible**
- 所以5-position conf有**更多的shell residues**！

---

## 数值分析

### Case A (8-position conf)

**Molecule for triple minimization:**
```
Minimized (3 residues):
  pos1=RC28, pos2=RC29, pos4=RC8

Fixed at template coords (shell):
  All other positions (0, 3, 5, 6, 7, and true shell residues)
```

**Energy calculation:**
```
E_A = Intra(pos1,RC28) + Intra(pos2,RC29) + Intra(pos4,RC8)
    + Inter(pos1-pos2) + Inter(pos1-pos4) + Inter(pos2-pos4)
    + Shell interactions (pos1,pos2,pos4 with all shell)
    + Shell self-energy (intra + inter among shell)
```

### Case B (5-position conf)

**Molecule for triple minimization:**
```
Minimized (3 residues):
  pos1=RC28, pos2=RC29, pos4=RC8

Fixed at template coords (shell):
  All other positions (0, 2, 3, and true shell residues)
```

**Wait, the shell should be the same!**

Let me re-examine. The calling context is:

**Case A:**
```
Full conf: [pos0=32, pos1=28, pos2=29, pos3=32, pos4=8, pos5=48, pos6=28, pos7=3]
```

Triple `[pos1=28, pos2=29, pos4=8]` is calculated.

**Case B:**
```
Partial conf: [pos0=32, pos1=28, pos2=29, pos3=9, pos4=8]
```

Triple `[pos1=28, pos2=29, pos4=8]` is calculated.

**When `makeMolecule([pos1=28, pos2=29, pos4=8])` is called:**
- It creates a molecule with pos1, pos2, pos4 using specified RCs
- ALL other positions (including pos0, pos3, pos5, pos6, pos7) use **template coordinates**

So the molecule should be **identical** in both cases! The positions 1, 2, 4 use the same RCs (28, 29, 8), and all other positions use template coords.

**Then why is the energy different?**

---

## 深入分析：Pairwise值的来源

让我重新看log中的"pairwise"值：

```java
// In MARKStarBound.java:1238-1243
double lowerbound = minimizingEmat.getInternalEnergy(tuple);

System.out.println("[TRIPLE_CORRECTION] tuple=" + tuple.stringListing()
    + " pairwise=" + String.format("%.6f", lowerbound)
    + " minimized=" + String.format("%.6f", tripleEnergy)
    ...
```

`lowerbound`来自`minimizingEmat.getInternalEnergy(tuple)`，这包括：

```java
// In UpdatingEnergyMatrix.java:113-153
public double getInternalEnergy(RCTuple tup) {
    double energy = 0;

    // Add one-body energies
    for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
        energy += getOneBody(posNum,RCNum);
    }

    // Add pairwise energies
    for(int indexInTuple=0; indexInTuple<numPosInTuple; indexInTuple++){
        for(int index2=0; index2<indexInTuple; index2++){
            energy += getPairwise(posNum,RCNum,pos2,rc2);
        }
    }

    // ⚠️ Add higher-order corrections!
    if (useHigherOrderTerms) {
        energy += internalEHigherOrder(tup);
    }

    return energy;
}
```

### 关键发现：Corrections Accumulation!

`getInternalEnergy()`会累加**已有的higher-order corrections**！

这意味着：
1. **第一次**计算triple `[1=28, 2=29, 4=8]`时：
   - `lowerbound = E_1(1,28) + E_1(2,29) + E_1(4,8) + E_2(1,2) + E_2(1,4) + E_2(2,4)`
   - 没有existing corrections
   - 计算minimized energy
   - 保存correction

2. **第二次**计算同一个triple时（在不同context）：
   - `lowerbound = E_1 + E_2 + **EXISTING_CORRECTION**`
   - 如果已经有correction，pairwise值就不同了！

但这还是不能完全解释为什么minimized值也不同...

---

## 等等！我发现问题了！

让我重新看两个contexts的**实际assignments**：

**Case A (Line 148):**
```
conf=32    28    29    32    8     48    28    3
pos: 0     1     2     3     4     5     6     7
```

**Case B (Line 2587):**
```
conf=32    28    29    9     8
pos: 0     1     2     3     4
```

**Triple in both cases:**
```
[pos1=28, pos2=29, pos4=8]
```

但是！在Case B的partial conformation中，**pos3=9**, 而在Case A中，**pos3=32**!

当创建molecule时，即使triple只包含[1,2,4]，**pos3并不在RCTuple中**，所以会使用template coordinates。

**但这意味着两次minimization的molecule应该是完全相同的！**

---

## 真相：不同的Energy Matrix状态

我认为答案在于：**两次计算时，Energy Matrix的状态不同！**

### Scenario

1. **时间点T1 (Line 148 附近):**
   - MARK*正在处理8-position conf
   - Energy matrix已经accumulate了一些corrections
   - 计算triple `[1=28, 2=29, 4=8]`
   - `pairwise = E_1 + E_2 + corrections_so_far = -13.964511`

2. **时间点T2 (Line 2587 附近):**
   - MARK*正在处理5-position conf
   - Energy matrix可能**被reset或更新**了
   - 或者使用**不同的context-specific matrix**
   - 再次计算同一个triple
   - `pairwise = E_1 + E_2 + different_corrections = -1.944728`

### 可能的机制

检查是否有**per-node energy matrix**或**context-dependent corrections**：

需要查看：
1. `MARKStarNode`是否维护自己的correction matrix？
2. Corrections是否是global的还是node-specific的？
3. 不同的partial conformations是否触发不同的matrix查询？

---

## 结论（目前推测）

Triple energy的context dependency可能来自：

1. **Energy Matrix的动态更新**
   - Corrections在MARK*运行期间不断累积
   - `getInternalEnergy()`包括已存在的corrections
   - 后续的triple计算会"看到"之前的corrections

2. **不同的Node Context**
   - 不同的A* nodes可能使用不同的correction sets
   - MARKStarNode可能维护node-specific的corrections

3. **Minimization的初始状态**
   - 虽然molecule构造应该相同，但minimizer的初始state可能不同
   - 可能受到之前minimizations的影响

需要进一步检查：
- MARKStarNode的correction matrix是否是per-node的
- 是否有多个UpdatingEnergyMatrix实例
- Corrections的scope和lifetime是什么
