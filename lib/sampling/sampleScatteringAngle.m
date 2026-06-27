function thetaSample = sampleScatteringAngle(thetaGrid, cdfS11, cdfS12, intensity, qRot)
% SAMPLESCATTERINGANGLE 基于预计算CDF基底的快速散射角采样
%
% 通过线性组合两个预计算CDF基底来采样散射角:
%   CDF = intensity * cumsum(S11*sinθ*dθ) + qRot * cumsum(S12*sinθ*dθ)
%
% 这等价于对偏振相位函数 (S11*I + S12*Q_rot)*sinθ 做逆CDF采样，
% 但避免了每次散射重复计算cumsum。
%
% INPUTS:
%   thetaGrid: [1×N] 散射角网格 (deg)
%   cdfS11:    [1×N] S11*(sinθ*dθ) 的累积和
%   cdfS12:    [1×N] S12*(sinθ*dθ) 的累积和
%   intensity: 标量, Stokes I 分量
%   qRot:      标量, 旋转后的 Q 分量 (= Q*cos2φ + U*sin2φ)
%
% OUTPUT:
%   thetaSample: 标量, 采样得到的散射角 (deg)

cdf = intensity * cdfS11 + qRot * cdfS12;
total = cdf(end);

if ~isfinite(total) || total <= 0
    thetaSample = thetaGrid(end);
    return;
end

u = rand(1, 1) * total;
idx = find(cdf > u, 1, 'first');
if isempty(idx)
    idx = numel(thetaGrid);
end

thetaSample = thetaGrid(idx);

end
