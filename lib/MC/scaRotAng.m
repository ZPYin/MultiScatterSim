function [scaAng, rotAng1, rotAng2] = scaRotAng(inDir, outDir)
% SCAROTANG 散射角度和旋转角度计算函数
%
% 该函数计算光子散射过程中的散射角和两个旋转角，这些角度是进行
% Stokes向量偏振变换的关键几何参数。在光子散射事件中，需要将
% Stokes向量从入射坐标系变换到散射坐标系，再变换到出射坐标系。
%
% 物理原理：
% 1. 散射角(scaAng)：入射方向与散射方向之间的夹角
% 2. 旋转角1(rotAng1)：入射坐标系到散射平面坐标系的旋转角
% 3. 旋转角2(rotAng2)：散射平面坐标系到出射坐标系的旋转角
%
% 几何定义：
% - 使用球坐标系(θ,φ)表示方向向量
% - θ：天顶角（与z轴夹角）
% - φ：方位角（在xy平面内与x轴夹角）
% - 散射平面：由入射方向和散射方向确定的平面
%
% USAGE:
%    [scaAng, rotAng1, rotAng2] = scaRotAng(inDir, outDir)
%
% INPUTS:
%    inDir: numeric [3×1] - 入射光子传播方向单位向量 [x; y; z]
%    outDir: numeric [3×1] - 散射光子传播方向单位向量 [x; y; z]
%
% OUTPUTS:
%    scaAng: numeric - 散射角 [弧度]
%        范围：[0, π]，0表示前向散射，π表示后向散射
%    rotAng1: numeric - 旋转角1 [弧度]
%        入射坐标系到散射平面的旋转角
%    rotAng2: numeric - 旋转角2 [弧度]  
%        散射平面到出射坐标系的旋转角
%
% 算法步骤：
% 1. 计算入射和出射方向的球坐标角度
% 2. 通过向量点积计算散射角
% 3. 处理特殊情况（平行散射、垂直方向等）
% 4. 计算旋转角的余弦值并进行数值稳定性处理
%
% REFERENCES:
%    [1] Hu Y-X, Winker D, Yang P, Baum B, Poole L, Vann L. Identification 
%        of cloud phase from PICASSO-CENA lidar depolarization: a multiple 
%        scattering sensitivity study. Journal of Quantitative Spectroscopy 
%        and Radiative Transfer. 2001;70(4):569-79.
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

% ========== 计算球坐标角度 ==========
% 计算入射和出射方向的方位角φ和天顶角θ
% 修正为 atan2(y, x) 能够避免除零NaN并且保留正确的几何象限
inPhi = atan2(inDir(2), inDir(1));      % 入射方位角[弧度]
outPhi = atan2(outDir(2), outDir(1));   % 出射方位角[弧度]
inTheta = acos(inDir(3));               % 入射天顶角 [弧度]
outTheta = acos(outDir(3));             % 出射天顶角 [弧度]

% ========== 计算散射角 ==========
% 使用向量点积公式：cos(散射角) = 入射向量 · 出射向量
dotProduct = sum(inDir .* outDir);      % 计算向量点积

% 处理数值精度问题，确保点积在[-1,1]范围内
if (dotProduct >= 1)
    scaAng = 0.0;                       % 完全前向散射
elseif (dotProduct <= -1)
    scaAng = pi;                        % 完全后向散射
else
    scaAng = acos(dotProduct);          % 一般情况的散射角
end

% ========== 特殊情况处理 ==========
% 情况1：散射角接近0或π（平行或反平行散射）
if (abs(sin(scaAng)) < 1e-12)
    rotAng1 = 0;                        % 旋转角设为0
    rotAng2 = 0;
    return;                             % 提前返回
end

% 情况2：入射方向接近垂直（z方向）
if (1 - abs(inDir(3))) <= 1e-12
    % 使用简化公式计算旋转角
    rotAng1 = acos(inDir(3) * cos(outPhi - inPhi));
    rotAng2 = acos(inDir(3));
    return;
end

% 情况3：出射方向接近垂直（z方向）
if (1 - abs(outDir(3))) <= 1e-12
    % 使用简化公式计算旋转角
    rotAng1 = acos(outDir(3) * cos(outPhi - inPhi));
    rotAng2 = acos(outDir(3));
    return;
end

% ========== 一般情况的旋转角计算 ==========
% 使用球面三角学公式计算旋转角的余弦值
% 基于散射平面几何关系的精确计算
rotAng1Cos = (outDir(3) - inDir(3) * cos(scaAng)) / (sin(inTheta) * sin(scaAng));
rotAng2Cos = (inDir(3) - outDir(3) * cos(scaAng)) / (sin(outTheta) * sin(scaAng));

% ========== 数值稳定性处理 ==========
% 限制余弦值在[-1,1]范围内，避免数值误差导致的计算错误
if rotAng1Cos >= 1
    rotAng1Cos = 1;                     % 上限截断
elseif rotAng1Cos <= -1
    rotAng1Cos = -1;                    % 下限截断
end

if rotAng2Cos >= 1
    rotAng2Cos = 1;                     % 上限截断
elseif rotAng2Cos <= -1
    rotAng2Cos = -1;                    % 下限截断
end

% ========== 计算最终旋转角 ==========
rotAng1 = acos(rotAng1Cos);            % 旋转角1 [弧度]
rotAng2 = acos(rotAng2Cos);            % 旋转角2 [弧度]

% ========== 恢复象限符号（找回3、4象限的方向性） ==========
% acos() 只返回 [0, π]，丢失了旋转的顺/逆时针方向。
% 使用 z 轴分量叉乘 (inDir × outDir)_z = inDir(1)*outDir(2) - inDir(2)*outDir(1)
% 判定真实的相对旋转方向，对 Q, U 偏振态积分至关重要。
rotSign = inDir(1) * outDir(2) - inDir(2) * outDir(1);
if rotSign < 0
    rotAng1 = -rotAng1;
    rotAng2 = -rotAng2;
end

end