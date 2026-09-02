# MSMC optical preprocessing

本目录保存 CPU 与 GPU 半解析偏振蒙特卡洛模型共同使用的云光学预处理。
它位于光子传输之前，把云微物理参数和 Mie 散射数据库转换为两个计算后端
可以直接使用的只读 `optics` 结构体。

## 模块边界

```text
Mie 散射数据库
+ lidarCfg.lambda
+ cloudCfg.clReff、N_d、kappa
          │
          ▼
prepareMSMCOptics
          │
          ├──► MSMCSim       （CPU 光子传输）
          └──► MSMCSimGPU    （GPU 光子传输）
```

本模块不执行光子自由传播、多次散射循环、Stokes 参考平面旋转、回波距离门
累加或微物理参数反演。

仓库自带用于最小运行示例的散射数据库，默认文件为
`data/water_cloud_scattering_phase_function_example.mat`；正式计算也可以通过
`scatMatrixFile` 显式指定其他兼容的高分辨率数据库。

## 文件说明

### `prepareMSMCOptics.m`

调用 `cloudScatteringMatrix` 对给定波长和液滴谱进行体积平均，生成：

- 分层消光系数 `muExt`；
- 分层散射系数 `muSca`；
- 单次散射反照率 `albedo`；
- Mueller 矩阵元素 `S11`、`S12`、`S33`、`S34`；
- 散射角网格 `theta`；
- 散射角抽样使用的 `cdfS11` 和 `cdfS12`。

### `prepareMSMCOpticsCached.m`

为批量正演和重复模拟提供内存/磁盘缓存。固定波长、`Reff/kappa` 廓线和
归一化 `N_d(z)` 形状时，绝对 `N_d` 幅度只线性改变 `muExt` 与
`muSca`；Mueller 矩阵和归一化散射角分布可以复用。

该函数主要用于：

- 大规模 `Reff-N_d` 参数扫描；
- 不同随机种子的重复蒙特卡洛模拟；
- 多个云底高度使用相同云内微物理廓线的批处理。

## 使用方式

推荐通过项目初始化脚本加载全部共享依赖：

```matlab
setupMSToolbox
```

也可以递归加入完整的 `lib`：

```matlab
addpath(genpath(fullfile(projectRoot, 'lib')));
```

普通调用：

```matlab
optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
    'scatMatrixFile', scatteringDatabase);

[SVcpu, diagCPU] = MSMCSim(lidarCfg, cloudCfg, ...
    'preparedOptics', optics);
[SVgpu, diagGPU] = MSMCSimGPU(lidarCfg, cloudCfg, ...
    'preparedOptics', optics);
```

缓存调用：

```matlab
[optics, cacheInfo] = prepareMSMCOpticsCached( ...
    lidarCfg, cloudCfg, ...
    'scatMatrixFile', scatteringDatabase);
```

CPU 和 GPU 必须共享同一份光学定义，不能在两个后端中分别维护消光系数、
Mueller 矩阵或散射角 CDF 的计算公式。
