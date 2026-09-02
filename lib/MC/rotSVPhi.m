function [SVOut] = rotSVPhi(SVIn, phi)
% ROTSVPHI Stokes向量的方位角旋转变换函数
%
% 该函数实现Stokes向量在偏振参考坐标系中的旋转变换。当光的传播方向或
% 观测坐标系发生旋转时，需要对Stokes向量进行相应的坐标变换以保持
% 偏振状态描述的一致性。
%
% 物理原理：
% Stokes向量 [I, Q, U, V] 描述了光的偏振状态：
% - I: 总强度（不受旋转影响）
% - Q: 线偏振分量（0°-90°方向）
% - U: 线偏振分量（45°-135°方向）  
% - V: 圆偏振分量（不受旋转影响）
%
% 当坐标系绕光传播方向旋转角度φ时，线偏振分量Q和U会发生耦合变换，
% 而总强度I和圆偏振V保持不变。这是因为线偏振的方向性与坐标系选择相关。
%
% USAGE:
%    [SVOut] = rotSVPhi(SVIn, phi)
%
% INPUTS:
%    SVIn: vector [4×1] - 输入Stokes向量 [I; Q; U; V]
%        I: 总强度分量
%        Q: 线偏振分量（0°-90°）
%        U: 线偏振分量（45°-135°）
%        V: 圆偏振分量
%    phi: numeric - 旋转角度 [弧度]
%        正值表示逆时针旋转，负值表示顺时针旋转
%
% OUTPUTS:
%    SVOut: vector [4×1] - 输出旋转后的Stokes向量 [I'; Q'; U'; V']
%
% 数学公式：
% 旋转变换矩阵为：
% [I']   [1    0         0        0] [I]
% [Q'] = [0  cos(2φ)  sin(2φ)    0] [Q]
% [U']   [0 -sin(2φ)  cos(2φ)    0] [U]
% [V']   [0    0         0        1] [V]
%
% 注意：旋转角度在变换中被放大2倍，这是Stokes参数的固有特性，
% 因为线偏振方向的π周期性对应Stokes参数的2π周期性。
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

% ========== Stokes向量旋转变换 ==========
% 应用旋转变换矩阵，实现坐标系旋转
SVOut = [SVIn(1); ...                                          % I分量：总强度保持不变
         SVIn(2) * cos(2 * phi) + SVIn(3) * sin(2 * phi); ... % Q分量：线偏振0°-90°变换
         -SVIn(2) * sin(2 * phi) + SVIn(3) * cos(2 * phi); ...% U分量：线偏振45°-135°变换  
         SVIn(4)];                                             % V分量：圆偏振保持不变

end