# GPU Monte Carlo backend

本目录保存半解析偏振蒙特卡洛模型的 GPU 计算实现。CPU 与 GPU 使用同一
套云光学参数、Mueller 矩阵和物理定义；二者是同一物理模型的两个计算后端。

## 公共入口

- CPU：`lib/MC/MSMCSim.m`
- GPU：`lib/MC_GPU/MSMCSimGPU.m`
- GPU 公共入口把参数转交给本目录中的 `MSMCSimGPUBackend.m`。
- GPU 不可用或配置不受支持时，可由 `fallbackToCPU` 控制是否回退到
  `MSMCSim`。

项目初始化脚本 `setupMSToolbox.m` 会递归加入 `lib`。

如果不调用初始化脚本，应递归加入完整的 `lib`，不能只加入
`lib/MC_GPU`，因为 GPU 后端仍需使用 `lib/MC_Optics` 中的共享光学
预处理，以及 `lib/MC` 中的 CPU 回退函数：

```matlab
addpath(genpath(fullfile(projectRoot, 'lib')));
```

## CPU 与 GPU 物理过程对应关系

| CPU 实现 | GPU 实现 | 物理过程 |
| --- | --- | --- |
| `freeRunning.m` | `freeRunningGPU.m` | 分层介质中的随机自由传播 |
| `rotSVPhi.m` | `rotSVPhiGPU.m` | Stokes 参考平面旋转 |
| `scaRotAng.m` | `scaRotAngGPU.m` | 散射角与两个参考平面旋转角 |
| `updateDir.m` | `updateDirGPU.m` | 散射后传播方向更新 |
| `sampleScatteringAngle.m` | `sampleScatteringAngleGPU.m` | 散射角分布抽样 |

GPU 函数一次处理成批的 `gpuArray` 光子状态，不能直接替换为面向单光子的


## 内部结构

```text
lib/
├── MC/                          CPU 光子传输
├── MC_Optics/                   CPU/GPU 共享光学预处理
└── MC_GPU/
    ├── MSMCSimGPU.m             GPU 公共入口
    ├── MSMCSimGPUBackend.m      GPU 调度与结果汇总
    ├── simulateLayeredBatchGPU.m
    ├── freeRunningGPU.m
    ├── rotSVPhiGPU.m
    ├── scaRotAngGPU.m
    ├── updateDirGPU.m
    ├── sampleScatteringAngleGPU.m
    └── README.md
```

- `MSMCSimGPUBackend.m`：参数解析、设备选择、光学量上传、批次调度、
  结果汇总、诊断输出和 CPU 回退。
- `simulateLayeredBatchGPU.m`：执行完整的光子散射循环，包括自由传播、
  散射抽样、Mueller 矩阵作用、接收判定和统计累加。该函数支持两种
  光子状态管理方式：
  - `legacy`（默认）：为整批光子保留固定长度的当前状态数组，并使用
    存活标记筛选继续传播的光子，便于与原实现逐项比较和调试。该模式
    不会自动保存每次散射的完整历史轨迹。
  - `packed`：每轮移除已经离开云层、被吸收或终止传播的光子，仅保留
    活跃光子的当前状态，以减少 GPU 索引回写和 CPU-GPU 同步开销。
    该模式用于优化大批量计算效率，实际加速效果取决于 GPU、批大小和
    诊断设置。
  两种方式采用相同的物理模型

完整 GPU 后端需要 MATLAB Parallel Computing Toolbox 和受支持的 GPU。
