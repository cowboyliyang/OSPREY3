OSPREY3 能量最小化中的并行机制详解
整体架构
用户代码 (COMETS/CometsZ)
        ↓
ConfEnergyCalculator.calcEnergyAsync()  （异步能量计算）
        ↓
TaskExecutor (ThreadPoolTaskExecutor - 4个线程)
        ↓
线程池 (Java ThreadPoolExecutor)
        ↓
工作线程 (pool-0-0, pool-0-1, pool-0-2, pool-0-3)
        ↓
EnergyCalculator.calcEnergy() → Minimizer → Energy Function
1. 任务提交层
源码：ConfEnergyCalculator.java:317-319
public void calcEnergyAsync(RCTuple frag, ResidueInteractions inters, 
                           TaskListener<EnergyCalculator.EnergiedParametricMolecule> listener) {
    tasks.submit(() -> calcEnergy(frag, inters), listener);
}
工作原理：
每个构象的能量计算作为一个**任务（task）**提交给执行器
任务是一个 lambda 函数，调用 calcEnergy()
立即返回（非阻塞），这样可以继续提交更多任务
2. 线程池任务执行器
源码：ThreadPoolTaskExecutor.java:90-120
public <T> void submit(Task<T> task, TaskListener<T> listener) {
    synchronized (this) {  // 线程安全的提交
        boolean wasAdded = false;
        while (!wasAdded) {
            checkException();
            
            // 尝试将任务加入队列，超时时间 400ms
            wasAdded = threads.submit(400, TimeUnit.MILLISECONDS, () -> {
                try {
                    // 在工作线程上运行任务
                    T result = runTask(task);
                    
                    // 将结果发送回监听器线程
                    threads.submitToListener(() -> {
                        taskSuccess(task, listener, result);
                    });
                } catch (Throwable t) {
                    taskFailure(task, listener, t);
                }
            });
        }
        startedTask();
    }
}
关键特性：
同步提交： 一次只提交一个任务（线程安全）
超时阻塞： 如果所有工作线程都忙，等待最多 400ms 后重试
独立监听器线程： 结果在专用线程上处理，避免阻塞工作线程
3. 线程池实现
源码：Threads.java:19-44
public Threads(int numThreads, int queueSize) {
    // 队列类型取决于 queueSize 参数
    if (queueSize <= 0) {
        queue = new SynchronousQueue<>();  // 无缓冲（默认）
    } else {
        queue = new ArrayBlockingQueue<>(queueSize);  // 缓冲任务
    }
    
    // 创建固定大小的线程池
    pool = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.DAYS, queue, (runnable) -> {
        Thread thread = Executors.defaultThreadFactory().newThread(runnable);
        thread.setDaemon(true);
        thread.setName(String.format("pool-%d-%d", poolId, threadId.getAndIncrement()));
        return thread;
    });
    pool.prestartAllCoreThreads();  // 立即启动所有线程
    
    // 用于回调的独立监听器线程
    listener = new ThreadPoolExecutor(1, 1, 0, TimeUnit.DAYS, new LinkedBlockingQueue<>(), ...);
    listener.prestartAllCoreThreads();
}
线程架构：
工作线程： numThreads 个（测试中是 4 个），命名为 pool-0-0, pool-0-1, pool-0-2, pool-0-3
监听器线程： 1 个线程，命名为 pool-0-listener，用于处理结果
队列： SynchronousQueue（默认）= 无缓冲，主线程会阻塞直到工作线程就绪
4. 能量计算与最小化
源码：EnergyCalculator.java:643-729
public EnergiedParametricMolecule calcEnergy(ParametricMolecule pmol, ResidueInteractions inters, 
                                             ResidueInteractionsApproximator approximator) {
    
    // 从自由度体素的中心开始
    DoubleMatrix1D x = DoubleFactory1D.dense.make(pmol.dofs.size());
    
    // 如果没有自由度或禁用最小化，返回刚性能量
    if (!isMinimizing || pmol.dofBounds.size() <= 0) {
        return new EnergiedParametricMolecule(pmol, efunc.getEnergy(), null);
    }
    
    // 可选：首先用纯 vdW 力场解决冲突
    if (alwaysResolveClashesEnergy != null) {
        Minimizer.Result vdwResult = minimizeWithVdw(pmol, inters, x);
        if (vdwResult.energy >= alwaysResolveClashesEnergy) {
            return new EnergiedParametricMolecule(pmol, Double.POSITIVE_INFINITY, null);
        }
        x = vdwResult.dofValues;
    }
    
    // 使用完整力场进行最小化
    try (EnergyFunction efunc = context.efuncs.make(ffInters, pmol.mol)) {
        MoleculeObjectiveFunction f = new MoleculeObjectiveFunction(pmol, efunc);
        
        try (Minimizer minimizer = context.minimizers.make(f)) {
            Minimizer.Result result = minimizer.minimizeFrom(x);
            
            // 检查是否陷入无限能量井（严重冲突）
            if (isInfiniteWell(result.energy)) {
                // 用 vdW 最小化重试
                Minimizer.Result vdwResult = minimizeWithVdw(pmol, ffInters, x);
                result = minimizer.minimizeFrom(vdwResult.dofValues);
                
                if (isInfiniteWell(result.energy)) {
                    return new EnergiedParametricMolecule(pmol, Double.POSITIVE_INFINITY, null);
                }
            }
            
            return new EnergiedParametricMolecule(pmol, result.energy, result.dofValues);
        }
    }
}
最小化过程（每个构象在工作线程上）：
初始化自由度 到体素中心
可选的冲突解决 使用纯 vdW 力场
主要最小化 使用 CCD（循环坐标下降）和完整力场
回退机制 如果检测到无限能量（严重冲突）
5. 并行实战示例
COMETS 示例，源码：Comets.java:421-455
// 获取下一批构象（批量大小 = 并行核心数）
confs.clear();
for (int i=0; i<state.confEcalc.tasks.getParallelism(); i++) {  // i=0,1,2,3（4核）
    ConfSearch.ScoredConf conf = confTree.nextConf();
    if (conf == null) break;
    confs.add(conf);
}

// 提交所有构象进行并行能量最小化
for (ConfSearch.ScoredConf conf : confs) {
    state.confEcalc.calcEnergyAsync(conf, confTable, econf -> {
        if (minEnergyConf == null || econf.getEnergy() < minEnergyConf.getEnergy()) {
            minEnergyConf = econf;
        }
    });
}

// 等待所有并行最小化完成
state.confEcalc.tasks.waitForFinish();
4核时间线：
时间 →
主线程:     [获取4个构象] [提交4个任务] ──────[等待]──────→ [处理结果]
工作线程0:                      [最小化构象0]────────→
工作线程1:                      [最小化构象1]────────→
工作线程2:                      [最小化构象2]────────→
工作线程3:                      [最小化构象3]────────→
监听器:                                        [回调0,1,2,3]
6. CPU vs GPU 并行
CPU 并行（测试中使用的）：
Parallelism.makeCpu(4)
→ ThreadPoolExecutor 有 4 个线程
→ SimpleCCDMinimizer（基于 CPU 的最小化器）
→ ResidueForcefieldEnergy（CPU 能量函数）
工作流程：
4 个 CPU 线程
每个线程独立最小化一个构象
使用 CCD（循环坐标下降）算法
串行处理每个构象的自由度
GPU 并行（可选）：
Parallelism.make(cpuCores=2, gpus=1, streamsPerGpu=16)
→ GpuStreamPool 有 16 个 CUDA 流
→ CudaCCDMinimizer（基于 GPU 的最小化器）
→ GpuForcefieldEnergy 或 ResidueForcefieldEnergyCuda
GPU 优势：
每个流可以在 GPU 上独立最小化一个构象
16 个流提供比 4 个 CPU 核心更高的吞吐量
GPU 并行处理能量计算（数千个 CUDA 核心）
7. 关键并行机制总结
并行类型
基于任务的并行： 每个构象 = 1 个任务
线程池： 固定大小的工作线程池（测试中 4 个）
异步提交： 非阻塞任务提交，带回调
批处理： 一次处理 N 个构象（N = 并行级别）
独立最小化： 每个构象独立最小化（完美并行）
同步点： waitForFinish() 阻塞直到所有任务完成
为什么是"完美并行"（Embarrassingly Parallel）
✅ 无共享状态： 每个构象的计算完全独立
✅ 无通信： 工作线程之间不需要通信
✅ 无依赖： 任务之间没有依赖关系
✅ 纯计算： 只是并行计算，没有复杂协调
性能特点
扩展性好： 添加更多核心接近线性加速
负载均衡： Java 线程池自动分配任务
资源利用： CPU/GPU 核心保持高利用率
内存友好： 每个构象独立，内存占用可控
测试中的实际配置
CPU 核心数：4
GPU：0（未使用）
队列大小：0（无缓冲，阻塞模式）
最小化算法：SimpleCCDMinimizer（CPU）
能量函数：ResidueForcefieldEnergy（CPU）
8. 实际运行示例
假设要计算 20 个构象的能量，使用 4 个 CPU 核心：
批次 1: [构象1, 构象2, 构象3, 构象4] → 4个工作线程并行处理
        等待全部完成...
        
批次 2: [构象5, 构象6, 构象7, 构象8] → 4个工作线程并行处理
        等待全部完成...
        
批次 3: [构象9, 构象10, 构象11, 构象12] → 4个工作线程并行处理
        等待全部完成...
        
批次 4: [构象13, 构象14, 构象15, 构象16] → 4个工作线程并行处理
        等待全部完成...
        
批次 5: [构象17, 构象18, 构象19, 构象20] → 4个工作线程并行处理
        等待全部完成...
总时间 ≈ (单个构象时间) × (总构象数 / 核心数) 如果单个构象需要 1 秒：
串行：20 秒
4 核并行：5 秒（4倍加速）
16 核并行：1.25 秒（16倍加速）
这就是为什么并行化对 OSPREY 这种计算密集型应用如此重要！