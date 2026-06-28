# MUltiple Scattering Simulator (MUSS)
MUSS是一个基于MATLAB的云层多次散射模拟程序，主要用于激光在水云中的传播、散射和偏振演化过程。程序结合Mie散射理论、云滴粒径谱和蒙特卡洛光子传输方法，可模拟激光在云层内发生多次散射后接收到的回波信号等结果。

主要模块                                                      说明
lib/MC/MSMCSim.m                                             蒙特卡洛多次散射主程序
lib/MC/freeRunning.m                                         云层中光子自由程采样
lib/MC/updateDir.m                                           光子传播方向更新
lib/MC/rotSVPhi.m                                            Stokes向量旋转
lib/scattering-phase-function/cloudScatteringMatrix.m        云层体平均Muller散射矩阵计算
lib/scattering-phase-function/water_droplet_scattering_gen.m 水滴Mie散射数据库生成函数
lib/cloud-microphysics/smmgDist.m                            修正Gamma云滴粒径谱
lib/cloud-optics/refIdxWater.m                               水的负折射率计算


## 下载
首先克隆仓库
git clone https://github.com/ZPYin/MultiScatterSim.git
进入项目根目录
cd MultiScatterSim-dev
在MATLAB中打开项目根目录，并运行脚本：setupMSToolbox.m 该脚本会将lib/，include/，analysis/和tests/等目录加入MATLAB搜索路径
（如何下载代码，300字以内）

## 运行
使用示例
参照\MultiScatterSim-dev\analysis\MCSim.m该脚本配置仿真输入参数,调用MSMCSim.m主程序，
该程序默认使用\MultiScatterSim-dev\data\water_cloud_scattering_phase_function_example.mat该Mie散射数据库，
Mie散射数据库根据自己仿真的需求可自行调用lib/scattering-phase-function/water_droplet_scattering_gen.m 该脚本生成。
配上一个已经写好的单层云退偏比测试脚本\MultiScatterSim-dev\analysis\depolarization_ration_cmp.m  


## 云层多次散射结果展示
![单层云线偏振退偏振比对比](img/depolarization_ration_cmp.jpg)
![参考文献](img/参考论文单层云去极化率.png)

![均匀、非均匀、部分均匀介质水云退偏比](img/去极化率图.jpg)
![参考文献](img/参考文献_均匀、非均匀、部分均匀介质退偏比.png)


## 多次散射模型基本原理

多次散射模型的详细原理说明见：[MUSS 基本原理介绍](./docs/MUSS基本原理介绍.md)

## Developer Team

|Name|Department|Email|
|:--:|-----:|:---:|
|Zhenping Yin|Wuhan University|zp.yin@whu.edu.cn|
|Chenwei Luo|Hubei University||cw.luo@stu.hubu.edu.cn|