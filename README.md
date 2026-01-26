# MUltiple Scattering Simulator (MUSS)

## TODO

1 优化拒绝采样方法（参考docs/On_sampling_of_scattering_phase_function.pdf）
2 优化计算效率(参考docs/amt-8-237-2015.pdf)
3 思考反演液水云云顶微物理特性的方法，应用在CALIPSO数据反演中

对于2，采用inverse sampling可以将采样速率提升两个数量级
对于3，使用Variance reduction方法可以让退偏结果快速收敛（约一个量级提升，可以使用更少的光子达到更高信噪比的结果）

## Developer Team

|Name|Department|Email|
|:--:|-----:|:---:|
|Zhenping Yin|Wuhan University|zp.yin@whu.edu.cn|
|Chenwei Luo|Hubei University||
