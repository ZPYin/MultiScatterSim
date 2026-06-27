function [n_r] = smmgDist(r, varargin)
% SMMGDIST 单模态修正伽马分布 - 云滴粒径谱计算函数
%
% 功能描述：
%    该函数基于单模态修正伽马分布（Single-Mode Modified Gamma Distribution）
%    计算云滴的粒径谱分布，用于描述云中水滴的数浓度随半径的变化关系。
%    这是激光雷达云微物理参数反演和蒙特卡罗散射仿真的重要基础函数。
%
% 物理意义：
%    修正伽马分布是描述云滴粒径谱的经典统计模型，能够准确反映自然云中
%    水滴尺寸的分布特征。该分布考虑了云滴的形状参数、有效半径和消光系数
%    等关键微物理参数，为激光雷达反演和散射计算提供真实的云滴分布。
%
% 使用方法：
%    [n_r] = smmgDist(r)                           % 使用默认参数
%    [n_r] = smmgDist(r, 'kappa', 10)             % 指定形状参数
%    [n_r] = smmgDist(r, 'R_eff', 12, 'sigma_e', 0.05)  % 指定多个参数
%
% 输入参数：
%    r: 数值型数组
%        云滴半径，单位：微米 (μm)
%        可以是标量或向量，表示需要计算数浓度的半径值
%
% 可选参数：
%    kappa: 数值型标量
%        形状参数，控制分布的形状特征 (默认值: 7)
%        较大的kappa值对应更窄的分布
%    sigma_e: 数值型标量
%        云的消光系数，单位：m^{-1} (默认值: 6e-2)
%        反映云的光学厚度特性
%    R_eff: 数值型标量
%        有效半径，单位：微米 (μm) (默认值: 9)
%        云滴的平均尺寸特征参数
%
% 输出参数：
%    n_r: 数值型数组
%        给定半径处的云滴数浓度，单位：m^{-3}μm^{-1}
%        表示单位体积、单位半径间隔内的云滴数量
%
% 应用场景：
%    1. 激光雷达云微物理参数反演
%    2. 蒙特卡罗散射仿真中的云滴分布建模
%    3. 云光学特性计算
%    4. 云物理过程数值模拟
%
% 参考文献：
%    [1] Donovan D, Klein Baltink H, Henzing J, De Roode S, Siebesma A. A
%        depolarisation lidar-based method for the determination of liquid-cloud
%        microphysical properties. Atmospheric Measurement Techniques. 2015;8(1)
%        :237--66.
%    [2] Ahmad W, Zhang K, Tong Y, Xiao D, Wu L, Liu D. Water Cloud Detection 
%        with Circular Polarization Lidar: A Semianalytic Monte Carlo Simulation
%        Approach. Sensors. 2022;22(4):1679.
%
% 历史记录：
%    2023-05-22: 首次版本 by Zhenping
% 作者联系方式：zp.yin@whu.edu.cn

% 参数解析器设置
p = inputParser;
p.KeepUnmatched = true;  % 保留未匹配的参数

% 定义必需参数和可选参数
addRequired(p, 'r', @isnumeric);                    % 必需参数：半径数组
addParameter(p, 'kappa', 7, @isnumeric);            % 可选参数：形状参数，默认值7
addParameter(p, 'sigma_e', 6e-2, @isnumeric);       % 可选参数：消光系数，默认值6e-2 m^{-1}
addParameter(p, 'R_eff', 9, @isnumeric);            % 可选参数：有效半径，默认值9 μm

% 解析输入参数
parse(p, r, varargin{:});

% 计算修正伽马分布的关键参数
% 根据参考文献[1]的公式(8)计算平均半径
R_m = p.Results.R_eff ./ (p.Results.kappa + 2);   % Eq.(8) Ref (1)

% 根据参考文献[1]的公式(11)计算k参数
k = (p.Results.kappa + 1) .* p.Results.kappa ./ (p.Results.kappa + 2).^2;   % Eq.(11) Ref (1)

% 根据参考文献[2]的公式(4)计算云滴总数浓度
N_d = 1 / (2 * pi) * p.Results.sigma_e ./ (k .* p.Results.R_eff.^2);   % Eq.(4) Ref(2)

% 计算单模态修正伽马分布的云滴数浓度
% 基于修正伽马分布公式：n(r) = N_d/R_m * (r/R_m)^(κ-1) * exp(-r/R_m) / Γ(κ)
% 其中 Γ(κ) 是伽马函数，对于任意实数κ > 0
n_r = N_d ./ R_m ./ gamma(p.Results.kappa) .* ...
      (r ./ R_m).^(p.Results.kappa - 1) .* exp(- r ./ R_m) * 1e12;

end