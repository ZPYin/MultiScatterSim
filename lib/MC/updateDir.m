function [outDir] = updateDir(inDir, phi, theta)
% UPDATEDIR 更新光子传播方向
%
% 功能描述：
%    根据散射角theta和相对方位角phi，计算光子散射后的新传播方向。
%    该函数实现了球坐标系中的方向向量变换，用于蒙特卡洛光子传输模拟。
%
% 物理原理：
%    在光子散射过程中，散射角theta定义为入射方向与散射方向的夹角，
%    相对方位角phi定义为散射平面相对于入射平面的旋转角。
%    通过球坐标变换公式，可以从入射方向和散射角度计算出散射方向。
%
% 用法：
%    [outDir] = updateDir(inDir, phi, theta)
%
% 输入参数：
%    inDir: 3×1向量
%        入射光子的传播方向单位向量 [x, y, z]
%    phi: 标量
%        相对方位角 (弧度)，范围 [0, 2π]
%    theta: 标量
%        散射角 (弧度)，范围 [0, π]
%
% 输出参数：
%    outDir: 3×1向量
%        散射后光子的传播方向单位向量 [x, y, z]
%
% 算法说明：
%    1. 当入射方向接近垂直（z分量接近±1）时，使用简化公式避免数值不稳定
%    2. 一般情况下，使用完整的球坐标变换公式
%    3. 变换公式基于Rodrigues旋转公式和球面三角学
%
% 参考文献：
%    [1] Hu Y-X, Winker D, Yang P, Baum B, Poole L, Vann L. Identification of 
%        cloud phase from PICASSO-CENA lidar depolarization: a multiple scattering 
%        sensitivity study. Journal of Quantitative Spectroscopy and Radiative 
%        Transfer. 2001;70(4):569-79.
%
% 历史版本：
%    2023-05-31: 初版 by Zhenping
% 作者: zp.yin@whu.edu.cn

% 判断入射方向是否接近垂直方向（z分量接近±1）
% 当 |z| ≈ 1 时，使用特殊处理避免数值不稳定性
if (1 - abs(inDir(3))) < 1e-12
    % 特殊情况：入射方向接近垂直（向上或向下）
    % 修正：保留 inDir 的方位角 φ_in，使散射方位角 φ 成为真正的相对角
    inPhi = atan2(inDir(2), inDir(1));  % 入射方向绝对方位角
    outDir = [sin(theta) * cos(inPhi + phi); ...       % x分量: 绝对方位=φ_in+φ
              sin(theta) * sin(inPhi + phi); ...       % y分量: 绝对方位=φ_in+φ
              cos(theta) * ((inDir(3) >= 0) - 0.5) * 2]; % z分量：保持原方向符号
    % 注：((inDir(3) >= 0) - 0.5) * 2 等价于 sign(inDir(3))，用于保持z方向的符号
else 
    % 一般情况：使用完整的球坐标变换公式
    % temp = sin(θ_in)，其中θ_in是入射方向与z轴的夹角
    temp = sqrt(1 - inDir(3)^2);  % temp = sqrt(x²+y²) = sin(θ_in)
    
    % 使用Rodrigues旋转公式计算散射后的方向向量
    % 该公式考虑了入射方向、散射角和方位角的完整几何关系
    outDir = [sin(theta) * (inDir(1) * inDir(3) * cos(phi) - inDir(2) * sin(phi)) / temp + inDir(1) * cos(theta); ...
              sin(theta) * (inDir(2) * inDir(3) * cos(phi) + inDir(1) * sin(phi)) / temp + inDir(2) * cos(theta); ...
              -sin(theta) * cos(phi) * temp + inDir(3) * cos(theta)];
    
    % 公式解释：
    % x分量 = sin(θ)[(x_in·z_in·cos(φ) - y_in·sin(φ))/sin(θ_in)] + x_in·cos(θ)
    % y分量 = sin(θ)[(y_in·z_in·cos(φ) + x_in·sin(φ))/sin(θ_in)] + y_in·cos(θ)
    % z分量 = -sin(θ)·cos(φ)·sin(θ_in) + z_in·cos(θ)
    %
    % 其中：
    % - 第一项表示散射引起的方向变化
    % - 第二项表示原方向在新方向上的投影
    % - temp = sin(θ_in) 用于归一化处理
end

end
