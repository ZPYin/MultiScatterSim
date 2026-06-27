function [pOutPos, isInMedium, ds, outClLayer] = freeRunning(pInPos, pInDir, cloudExt, dClThickness, cloudBaseHeight, cloudBoundX, cloudBoundY)
% FREERUNNING 计算光子在云层中的下一步自由传播
% 
% 该函数实现了蒙特卡洛方法中光子在多层云介质中的自由传播过程。
% 基于Beer-Lambert定律和随机采样，计算光子在发生下一次散射事件前的传播距离和位置。
%
% USAGE:
%    [pOutPos, isInMedium, ds] = freeRunning(pInPos, pInDir, cloudExt, dClThickness, cloudBaseHeight, cloudBoundX, cloudBoundY)
%
% INPUTS:
%    pInPos: numeric [1x3]
%        入射光子的三维位置坐标 [x, y, z] (m)
%    pInDir: numeric [1x3]  
%        光子传播方向的单位向量 [dx, dy, dz]
%    cloudExt: numeric [Nx1]
%        各云层的消光系数数组 (m^-1)，N为云层数量
%    dClThickness: numeric
%        单个云层的垂直厚度 (m)
%    cloudBaseHeight: numeric
%        云底高度 (m)
%    cloudBoundX: numeric [Nx1]
%        各云层在x方向的边界半宽度 (m)
%    cloudBoundY: numeric [Nx1]
%        各云层在y方向的边界半宽度 (m)
%
% OUTPUTS:
%    pOutPos: numeric [1x3]
%        光子传播后的新位置坐标 [x, y, z] (m)
%    isInMedium: logical
%        逻辑变量，判断光子是否仍在云介质内部
%    ds: numeric
%        光子的传播距离 (m)
%    outClLayer: numeric
%        光子传播后所在的云层索引
%
% ALGORITHM:
%    1. 使用指数分布随机采样确定光学厚度 tau0 = -ln(rand)
%    2. 根据光子传播方向分为水平传播和倾斜传播两种情况
%    3. 对于倾斜传播，考虑多层云的光学厚度累积效应
%    4. 计算光子新位置并判断是否超出云边界
%
% HISTORY:
%    2023-06-02: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

% 初始化输出变量
pOutPos = NaN(1, 3);        % 光子输出位置，初始化为NaN
isInMedium = true;          % 光子是否在介质中，初始为true
ds = 0;                     % 传播距离，初始为0

% 使用指数分布随机采样生成光学厚度
% tau0 = -ln(rand) 符合指数分布，用于确定到下次散射事件的光学距离
tau0 = - log(rand(1, 1));

% 计算光子当前所在的云层索引
% 使用ceil函数向上取整，加1e-8避免数值精度问题
% 钳制到有效层数范围，防止云顶边界溢出
inClLayer = ceil((pInPos(3) - cloudBaseHeight) / dClThickness + 1e-8);
inClLayer = max(1, min(inClLayer, length(cloudExt)));

% 判断光子传播方向：水平传播 vs 倾斜传播
if abs(pInDir(3)) < 1e-6
    % ========== 水平传播情况 ==========
    % 当光子垂直方向分量很小时，认为是水平传播
    % 在单一云层内传播，传播距离 = 光学厚度 / 消光系数
    ds = tau0 / cloudExt(inClLayer);
    pOutPos = pInPos + pInDir * ds;     % 计算新位置
    outClLayer = inClLayer;             % 输出云层索引与输入相同

else
    % ========== 倾斜传播情况 ==========
    % 光子可能穿越多个云层，需要考虑累积光学厚度
    
    % 计算各云层的累积光学厚度
    % cumsumTau(i) = 从第1层到第i层的总光学厚度
    cumsumTau = cumsum(cloudExt .* dClThickness);

    % 计算光子当前位置对应的累积光学厚度 tau1
    if inClLayer > 1
        % 如果不在第一层，需要加上前面所有层的光学厚度
        tau1 = (sum(cloudExt(1:(inClLayer - 1)) * dClThickness)) + ...
               (pInPos(3) - cloudBaseHeight - (inClLayer - 1) * dClThickness) * cloudExt(inClLayer);
    else
        % 如果在第一层，只计算当前层内的光学厚度
        tau1 = (pInPos(3) - cloudBaseHeight - (inClLayer - 1) * dClThickness) * cloudExt(inClLayer);
    end
    
    % 计算光子传播后的累积光学厚度 tau2
    % tau0 是沿光子路径的光学厚度，换算到垂直累计光学厚度需要乘以方向余弦
    tau2 = tau1 + tau0 * pInDir(3);
    
    % 查找光子最终所在的云层
    % 找到第一个累积光学厚度大于tau2的云层
    idx = find(tau2 < cumsumTau, 1);
    if isempty(idx)
        % 如果没找到，说明光子已经穿出所有云层
        outClLayer = NaN;
    else
        % 找到目标云层，计算精确的垂直位置
        outClLayer = idx;
        % 根据光学厚度反推垂直位置
        pOutPos(3) = idx * dClThickness + cloudBaseHeight + (tau2 - cumsumTau(idx)) / cloudExt(idx);
        % 根据垂直位移计算总传播距离
        ds = (pOutPos(3) - pInPos(3)) / pInDir(3);
        % 计算完整的三维新位置
        pOutPos = pInPos + pInDir * ds;
    end

end

% ========== 边界条件检查 ==========
% 判断光子是否仍在云介质内部，需要检查多个边界条件
if isnan(outClLayer)
    % 情况1：光子已穿出所有云层（垂直方向）
    isInMedium = false;
elseif (abs(pOutPos(1)) > cloudBoundX(min(outClLayer, length(cloudBoundX)))) || ...
       (abs(pOutPos(2)) > cloudBoundY(min(outClLayer, length(cloudBoundY))))
    % 情况2：光子超出当前云层的水平边界
    % 检查x和y方向是否超出对应云层的边界范围
    isInMedium = false;
elseif (pOutPos(3) > (cloudBaseHeight + dClThickness * length(cloudExt))) || (pOutPos(3) < cloudBaseHeight)
    % 情况3：光子超出整个云系统的垂直边界
    % 高于云顶或低于云底都认为离开介质
    isInMedium = false;
end

end
