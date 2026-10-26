# MUltiple Scattering Simulator (MUSS)

MUSS 是基于 MATLAB 的水云偏振激光雷达多次散射模拟程序，结合 Mie 散射与半解析蒙特卡洛方法，计算回波 Stokes 信号、强度和线偏振退偏比。

模型支持分层云配置，以及消光系数或云滴数浓度两种参数驱动方式，提供共用光学预处理的 CPU 和 GPU 实现。

## 快速开始

CPU 计算需要 MATLAB；GPU 计算还需要 Parallel Computing Toolbox 和受 MATLAB 支持的 GPU。

获取仓库的 `master` 分支：

```bash
git clone --branch master https://github.com/ZPYin/MultiScatterSim.git
cd MultiScatterSim
```

在 MATLAB 中将当前文件夹切换到项目根目录，初始化并运行 CPU 示例：

```matlab
setupMSToolbox
run(fullfile('analysis', 'MCSim.m'))
```

使用 GPU 时，运行：

```matlab
run(fullfile('analysis', 'MCSimGPU.m'))
```

两个示例均使用 532 nm 波长、1 万个光子和仓库内置的[示例散射数据库](data/water_cloud_scattering_phase_function_example.mat)。可在示例脚本中修改 `lidarCfg`、`cloudCfg` 和 `nPhoton`。

详细输入参数、输出格式和退偏比公式见[模型原理与参数说明](docs/MUSS基本原理介绍.md)。

## 与 MSCART 的结果对比

提供地基和星载两个案例，将 MUSS 结果与采用第三方 MSCART 模型计算得到的参考数据进行比较。MSCART 的源代码和可执行程序不随本仓库提供。

<details>
<summary>查看地基和星载对比图（英文）</summary>

**地基案例：** 355 nm 波长、5 m 距离门，低/高消光条件，全视场为 0.5 和 2.0 mrad；展示云内累计至 150 m 的线偏振退偏比。

<p align="center">
  <img src="analysis/mscart_comparison/ground/comparison_en.png" alt="MUSS 与 MSCART 地基案例对比" width="800">
</p>

**星载案例：** 532 nm 波长、3 m 距离门，全视场为 0.001、0.13 和 13 mrad；展示云内累计至 300 m 的线偏振退偏比。

<p align="center">
  <img src="analysis/mscart_comparison/satellite/comparison_en.png" alt="MUSS 与 MSCART 星载案例对比" width="1000">
</p>

蓝色实线为 MUSS，黑色虚线与圆点为 MSCART。两个偏振通道先分别对独立重复取平均，再累加取比；阴影表示独立重复的取值范围。

</details>

已有数据、统计定义和绘图步骤见[MSCART 对比说明](analysis/mscart_comparison/README.md)。

## 文档与示例

| 文档 | 内容 |
| --- | --- |
| [模型原理与参数说明](docs/MUSS基本原理介绍.md) | 物理流程、输入输出、调用示例及模块索引 |
| [GPU 模型说明](lib/MC_GPU/README.md) | GPU 依赖、计算结构和状态管理方式 |
| [光学预处理说明](lib/MC_Optics/README.md) | 散射数据库、光学参数预计算与缓存 |
| [MSCART 对比说明](analysis/mscart_comparison/README.md) | 参考数据与对比图绘制 |
| [单层云退偏比示例](analysis/depolarization_ration_cmp.m) | 不同消光系数下退偏比随云内深度的变化 |

## 开发团队

| 姓名 | 单位 | 邮箱 |
| --- | --- | --- |
| Zhenping Yin | Wuhan University | zp.yin@whu.edu.cn |
| Chengwei Luo | Hubei University | cw.luo@stu.hubu.edu.cn |
