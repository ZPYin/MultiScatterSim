function [SVReturn, diag] = MSMCSim(lidarCfg, cloudCfg, varargin)
% MSMCSIM 基于蒙特卡洛方法的激光雷达多次散射模拟
%
% 该函数实现了激光雷达在云层中的多次散射过程的蒙特卡洛模拟，考虑了光子的
% 偏振特性（Stokes向量）和云滴的Mie散射特性。模拟包括光子在云层中的传播、
% 散射、吸收以及最终被激光雷达接收的过程。
%
% 算法原理：
% 1. 光子发射：在云底按激光束截面随机生成光子初始位置
% 2. 自由传播：使用Beer-Lambert定律计算光子到下次散射事件的距离
% 3. 散射事件：基于Mie散射理论和相位函数进行散射角度采样
% 4. 偏振计算：使用Mueller散射矩阵更新Stokes向量
% 5. FOV判断：检查散射光子是否进入激光雷达视场
% 6. 统计累积：将接收到的光子信号按距离分档累积
%
% USAGE:
%    [SVReturn] = MSMCSim(lidarCfg, cloudCfg)
%    [SVReturn] = MSMCSim(lidarCfg, cloudCfg, 'nPhoton', 1e5, 'visible', 'off')
%
% INPUTS:
%    lidarCfg: struct - 激光雷达系统配置参数
%        lambda: numeric - 激光波长 [nm]
%        rBeam: numeric - 激光束在z=0处的半径 [mm]
%        divBeam: numeric - 激光束发散角 [mrad]
%        FOV: numeric - 接收视场角 [mrad]
%        rTel: numeric - 接收望远镜半径 [mm]
%        binWidth: numeric - 距离分辨率 [ns]
%        nBins: numeric - 距离档数
%        distBeamTel: numeric - 发射光轴与接收光轴间距 [mm]
%        SBeam: numeric [4x1] - 入射激光的Stokes向量 [I,Q,U,V]
%    cloudCfg: struct - 云层配置参数
%        clBase: numeric - 云底高度 [m]
%        clDh: numeric - 单层云几何厚度 [m]
%        nClLayers: numeric - 云层数量
%        kappa: numeric - 修正伽马分布形状因子
%        clExt: numeric [Nx1] - 各云层消光系数 [m^-1]
%        clReff: numeric [Nx1] - 各云层有效半径 [μm]
%
% KEYWORDS:
%    nPhoton: numeric - 模拟光子数量 (默认: 1e4)
%    visible: char - 是否显示结果图 'on'|'off' (默认: 'on')
%    debug: logical - 是否输出调试信息 (默认: false)
%
% OUTPUTS:
%    SVReturn: numeric [4 x nBins] - 各距离档的Stokes向量返回信号
%        第1行: I分量（总强度）
%        第2行: Q分量（线偏振度，0°-90°）
%        第3行: U分量（线偏振度，45°-135°）
%        第4行: V分量（圆偏振度）
%
% ALGORITHM DETAILS:
%    采用前向蒙特卡洛方法，结合半解析技术提高计算效率：
%    - 使用指数分布采样确定散射事件间距
%    - 基于云滴粒径分布和Mie理论计算散射特性
%    - 考虑多层云的光学性质差异
%    - 实现偏振光的完整传输过程
%    - 采用强制散射技术增加FOV内的有效光子数
%
% REFERENCES:
%    - Winker & Poole (1995), Applied Physics B, 60(4), 341-344
%    - Mishchenko et al. (2002), Scattering, Absorption, and Emission 
%      of Light by Small Particles
%
% HISTORY:
%    2023-06-05: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

% ========== 输入参数解析 ==========
% 使用inputParser进行参数验证和默认值设置
p = inputParser;
p.KeepUnmatched = true;   % 保留未匹配的参数，便于扩展

% 必需参数验证
addRequired(p, 'lidarCfg', @isstruct);   % 激光雷达配置必须为结构体
addRequired(p, 'cloudCfg', @isstruct);   % 云层配置必须为结构体

% 可选参数及默认值
addParameter(p, 'nPhoton', 1e4, @isnumeric);      % 模拟光子数，默认1万个
addParameter(p, 'visible', 'on', @ischar);        % 结果可视化开关
addParameter(p, 'debug', false, @islogical);      % 调试信息开关
addParameter(p, 'useRussianRoulette', true, @islogical);      % 是否使用轮盘对赌
addParameter(p, 'energyThreshold', 0.01, @isnumeric);         % 轮盘对赌能量阈值
addParameter(p, 'recordHitMap', false, @islogical);
addParameter(p, 'mapNx', 64, @isnumeric);
addParameter(p, 'mapNy', 64, @isnumeric);
addParameter(p, 'mapRadius', NaN, @isnumeric);
addParameter(p, 'nAngleGroups', 5, @isnumeric);
addParameter(p, 'useLocalAtten', false, @islogical); % 调试：使用本层局部透过而非累计垂直光学深度
addParameter(p, 'scatMatrixFile', '', @ischar);      % 可选：显式指定水滴Mie散射查找表
addParameter(p, 'angleSampler', 'fast', @ischar);    % fast|inverse
% 执行参数解析
parse(p, lidarCfg, cloudCfg, varargin{:});

%% ========== 参数定义 ==========
nPhotons = p.Results.nPhoton;   % 获取模拟光子总数
angleSampler = lower(p.Results.angleSampler);
if ~ismember(angleSampler, {'fast', 'inverse'})
    error('angleSampler must be fast or inverse.');
end

diag = struct();
diag.angleSampler = angleSampler;
if p.Results.recordHitMap
    Nx = p.Results.mapNx;
    Ny = p.Results.mapNy;
    mapRadius = p.Results.mapRadius;
    if ~isfinite(mapRadius) || mapRadius <= 0
        mapRadius = lidarCfg.rTel * 1e-3;
    end
    nAngleGroups = round(p.Results.nAngleGroups);
    if nAngleGroups < 1
        nAngleGroups = 1;
    end
    diag.mapNx = Nx;
    diag.mapNy = Ny;
    diag.mapRadius = mapRadius;
    diag.nAngleGroups = nAngleGroups;
    maxAngle = lidarCfg.FOV * 1e-3 / 2;
    diag.angleBinEdges = linspace(0, maxAngle, nAngleGroups + 1);
    diag.angleBinCenters = 0.5 * (diag.angleBinEdges(1:end-1) + diag.angleBinEdges(2:end));

    xVals = linspace(-mapRadius, mapRadius, Nx);
    yVals = linspace(-mapRadius, mapRadius, Ny);
    [Xgrid, Ygrid] = meshgrid(xVals, yVals);
    mask = (Xgrid.^2 + Ygrid.^2) <= mapRadius^2;
    diag.telMapI = zeros(Ny, Nx);
    diag.telMapQ = zeros(Ny, Nx);
    diag.telMapU = zeros(Ny, Nx);
    diag.telMapV = zeros(Ny, Nx);
    diag.telMapI_groups = zeros(Ny, Nx, nAngleGroups);
    diag.telMapQ_groups = zeros(Ny, Nx, nAngleGroups);
    diag.telMapU_groups = zeros(Ny, Nx, nAngleGroups);
    diag.telMapV_groups = zeros(Ny, Nx, nAngleGroups);
    diag.mask = mask;
end

%% ========== 云层散射特性计算 ==========
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

% 添加必要的路径
rootPath = fullfile(projectDir, 'lib');
addpath(fullfile(rootPath, 'sampling'));
addpath(fullfile(rootPath, 'scattering-phase-function'));
addpath(fullfile(rootPath, 'cloud-microphysics'));

% 项目根目录与数据库路径
clScatMatrixFile = fullfile(projectDir, 'tests', 'data', ...
    'water_droplet_scattering_532nm_finer_validation_r_step_0p01um.mat');
if ~isempty(p.Results.scatMatrixFile)
    clScatMatrixFile = p.Results.scatMatrixFile;
end

% 计算云层的散射矩阵和光学特性
% 基于Mie散射理论，考虑云滴粒径分布（修正伽马分布）
dist = 'modified-gamma';
unmatchedArgs = p.Unmatched;  % 复制一份，避免修改只读属性
if isfield(unmatchedArgs, 'cloud_distribution')
    dist = unmatchedArgs.cloud_distribution;
    unmatchedArgs = rmfield(unmatchedArgs, 'cloud_distribution');
end
extraArgs = reshape([fieldnames(unmatchedArgs).'; struct2cell(unmatchedArgs).'], 1, []);
cloudScatMatrix = cloudScatteringMatrix(lidarCfg.lambda, clScatMatrixFile, ...
    'cloud_distribution', dist, ...               % 支持 modified-gamma 或 modified-gamma-paper
    'kappa', cloudCfg.kappa, ...                  % 分布形状参数
    'sigma_e', cloudCfg.clExt, ...                % 消光系数
    'R_eff', cloudCfg.clReff, ...                 % 有效半径
    'b', 1/9, ...                                 % 强制传入给 modified-gamma-paper
    extraArgs{:});                                % 透传其他可选参数（如 b）

% 提取关键光学参数
muExt = cloudScatMatrix.muExt;     % 消光系数 [m^-1]
muSca = cloudScatMatrix.muSca;     % 散射系数 [m^-1]  
albedo = muSca ./ muExt;           % 单次散射反照率 (散射系数/消光系数)

% 提取Mueller散射矩阵元素（用于偏振计算）
S11 = cloudScatMatrix.S11;        % 散射相位函数主对角元素
S12 = cloudScatMatrix.S12;        % 线偏振相关的散射矩阵元素
S33 = cloudScatMatrix.S33;        % 线偏振保持相关元素
S34 = cloudScatMatrix.S34;        % 线偏振-圆偏振转换元素

% 预计算随机散射角采样所需的CDF基底。
% 每次散射的PDF为 (S11*I + S12*Q_rot)*sin(theta)，因此CDF可由两组基底线性组合，省略的cumsum函数的调用次数。
thetaGrid = cloudScatMatrix.theta;
thetaWeights = [0, diff(thetaGrid)] .* sin(thetaGrid * pi / 180);
cdfS11 = cumsum(S11 .* thetaWeights, 2);
cdfS12 = cumsum(S12 .* thetaWeights, 2);

% ========== 云层边界设置 ==========
% 根据Winker & Poole (1995)的建议，云边界对应光学厚度为2
% 这里设置为4，提供更大的模拟空间以捕获边缘效应
clBoundX = 4.0 ./ muExt;   % X方向云层边界半宽度 [m]
clBoundY = 4.0 ./ muExt; % Y方向云层边界半宽度 [m]

%% ========== 蒙特卡洛模拟主循环 ==========
% 计算激光束在云底处的半径
rCl = lidarCfg.rBeam * 1e-3 + lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase * 1e-3;   % 云底激光束半径 [m]

% 初始化输出变量
SVReturn = zeros(4, lidarCfg.nBins);   % 返回信号的Stokes向量 [4×nBins]
diag.SVReturnLayer = zeros(4, lidarCfg.nBins, cloudCfg.nClLayers); % 分层返回累积
len = zeros(1, nPhotons);              % 每个光子的传播距离 [m]
nScas = zeros(1, nPhotons);            % 每个光子的散射次数
layerScatterCount = zeros(1, cloudCfg.nClLayers); % 分层散射事件计数
% 开始光子循环模拟
for iPhoton = 1:nPhotons
    % ========== 光子初始化 ==========
    pPos = [0; 0; 0];                  % 光子位置向量 [x; y; z] [m]
    temp = zeros(1, 5);                % 临时计算变量

    % ========== 光子发射到云底 ==========
    isAlive = true;                    % 光子存活标志
    
    % 在激光束截面内随机生成光子初始位置
    pPos(1) = rand(1, 1) * rCl;       % X方向随机位置
    pPos(2) = rand(1, 1) * rCl;       % Y方向随机位置
    pPos(3) = cloudCfg.clBase;         % Z方向设为云底高度
    
    % 确保光子位置在圆形激光束截面内
    while (norm(pPos(1:2), 2) > rCl)
        % 重新生成位置直到满足圆形边界条件
        pPos(1) = rand(1, 1) * rCl;
        pPos(2) = rand(1, 1) * rCl;
    end
    
    % 计算光子初始传播方向（从激光器指向初始位置）
    pDir = pPos / norm(pPos - [0; 0; 0], 2);
    
    % 根据初始传播方向旋转入射Stokes向量
    % 修正为 atan2，保证与 scaRotAng 中的方位角定义在四个象限匹配，且避免除零
    SVIn = rotSVPhi(lidarCfg.SBeam, atan2(pDir(2), pDir(1)));

    % ========== 光子传播循环 ==========
    while (isAlive)

        % ========== 自由传播步骤 ==========
        % 调用freeRunning函数计算光子下一步位置
        % 输入：当前位置、方向、消光系数、云层厚度、云底高度、边界
        % 输出：新位置、是否仍在介质中、传播距离、散射云层索引
        [pPos, isInMedium, ds, scaClLayer] = freeRunning(pPos, pDir, muExt, cloudCfg.clDh, cloudCfg.clBase, clBoundX, clBoundY);
        len(iPhoton) = len(iPhoton) + ds;   % 累积传播距离

        % ========== 边界检查 ==========
        if (isInMedium)
            % 光子仍在云层内，进行吸收判断
            if (rand(1, 1) > albedo(scaClLayer))
                % 光子被吸收，结束传播
                isAlive = false;
                continue;
            end
            
            % ========== 轮盘对赌(Russian Roulette) ==========
            % 当光子强度低于阈值时激活轮盘对赌，减少低权重光子的计算量
            if p.Results.useRussianRoulette && SVIn(1) < p.Results.energyThreshold
                % 生存概率正比于当前能量与阈值的比值
                survivalProbability = SVIn(1) / p.Results.energyThreshold;
                
                if rand(1, 1) > survivalProbability
                    % 光子被终止
                    isAlive = false;
                    if p.Results.debug && mod(iPhoton, 1000) == 0
                        fprintf('光子 %d: 轮盘对赌终止, 当前能量=%.4f\n', iPhoton, SVIn(1));
                    end
                    break;
                else
                    % 光子存活，按比例增加权重以保持统计正确性
                    SVIn = SVIn / survivalProbability;
                    if p.Results.debug && mod(iPhoton, 10000) == 0
                        fprintf('光子 %d: 轮盘对赌存活, 权重增加到 %.4f\n', iPhoton, SVIn(1));
                    end
                end
            end

            % ========== FOV内散射检查和半解析方法 ==========
            % 判断光子是否在激光雷达视场角(FOV)内
            % FOV条件：光子到望远镜光轴的距离 < FOV半径
            % FOV半径 = (FOV角度 × 高度 + 望远镜孔径) / 2
            if (sqrt((pPos(1) - lidarCfg.distBeamTel * 1e-3)^2 + pPos(2)^2) < ((lidarCfg.FOV * 1e-3 * pPos(3) + lidarCfg.rTel * 1e-3 * 2) / 2))
                
                % ========== 强制散射到FOV方向（半解析方法）==========                                                      
                % 计算从当前位置指向望远镜的单位方向向量
                pScaDir = [lidarCfg.distBeamTel * 1e-3; 0; 0] - pPos;
                pScaDir = pScaDir / norm(pScaDir, 2);

                % ========== 计算散射角度 ==========
                % 计算散射角和旋转角，用于Stokes向量变换
                [scaAng,rotAng1, rotAng2] = scaRotAng(pDir, pScaDir); %散射角scaAng就是pDir入射方向与pScaDir散射方向之间的夹角  rotAng1和rotAng2是两个旋转角  用于Stokes向量变换

                % 记录层内散射计数（FOV强制散射）
                layerScatterCount(scaClLayer) = layerScatterCount(scaClLayer) + 1;

                % ========== 计算进入FOV的概率 ==========
                % 根据散射相位函数查找对应角度的索引
                [~, iAng] = min(abs(cloudScatMatrix.theta - scaAng / pi * 180));  %查找对应散射矩阵中散射角度最近的索引  方便后续根据这个索引计算对应角度的散射相位函数值 例如S11 S12等
                
                % 计算垂直光学厚度（考虑多层云累积消光）
                % 1. 下方完整层的光学厚度
                opticalDepthBelow = 0;
                if scaClLayer > 1
                    opticalDepthBelow = sum(muExt(1:scaClLayer-1)) * cloudCfg.clDh;
                end
                
                % 2. 当前层光子位置到底部的光学厚度
                layerBottom = cloudCfg.clBase + (scaClLayer - 1) * cloudCfg.clDh;
                opticalDepthCurrent = (pPos(3) - layerBottom) * muExt(scaClLayer);
                
                % 3. 总垂直光学厚度
                verticalOD = opticalDepthBelow + opticalDepthCurrent;
                % 选择透过光学厚度：默认用累计垂直OD；调试时可仅用本层局部OD，减弱跨层衰减
                attenOD = verticalOD;
                if p.Results.useLocalAtten
                    attenOD = opticalDepthCurrent;
                end

                if cloudScatMatrix.theta(iAng) < 90
                    % 前向散射情况，概率设为0（简化处理）
                    probEnter = 0;
                else
                    % 后向散射情况，计算详细的进入概率
                    % 考虑：单次散射反照率、立体角、传输损失、相位函数
                    probEnter = albedo(scaClLayer) * ...
                        pi * (lidarCfg.rTel * 1e-3)^2 / norm(pPos - [lidarCfg.distBeamTel * 1e-3; 0; 0], 2).^2 * ...
                        exp(-attenOD / abs(pScaDir(3))) * ...
                        S11(scaClLayer, iAng);
                end

                % ========== 计算散射后的Stokes向量 ==========
                % 第一步：根据散射几何旋转入射Stokes向量
                SVTemp1 = rotSVPhi(SVIn, rotAng1);
                
                % 第二步：应用Mueller散射矩阵进行偏振变换
                SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                           S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                           0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                           0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                        SVTemp1 / S11(scaClLayer, iAng);

                % ========== 计算返回时间和距离仓位 ==========
                % 计算光子从散射点返回到激光雷达的总时间
                % 包括：已传播距离 + 返回距离（考虑垂直分量）
                timeBin = (len(iPhoton) + (pPos(3) - cloudCfg.clBase) / ...
                    max(abs(pScaDir(3)), eps)) / ...
                    3e8 / (lidarCfg.binWidth * 1e-9);
                iBin = ceil(timeBin);
                
                % ========== 累积返回信号 ==========
                if isfinite(iBin) && iBin >= 1 && iBin <= lidarCfg.nBins
                    % 第三步：旋转到最终参考系并累积到对应时间仓
                    % 权重包括：进入概率 × 多次散射增强因子
                    SV_final = rotSVPhi(rotSVPhi(SVTemp2, -rotAng2), atan2(pScaDir(2), pScaDir(1))) * probEnter* 2^(nScas(iPhoton) - 1);
                    SVReturn(:, iBin) = SVReturn(:, iBin) + SV_final;
                    diag.SVReturnLayer(:, iBin, scaClLayer) = diag.SVReturnLayer(:, iBin, scaClLayer) + SV_final;
                end

                % ========== 更新剩余光子强度 ==========
                % 减去已经"逃逸"到FOV的光子部分
                SVIn = SVIn * (1 - probEnter);

            end

            % ========== 随机散射事件处理 ==========
            % 对于不在FOV内或未被半解析方法处理的光子，进行随机散射
            
            % ========== 随机采样散射角度 ==========
            % 第一步：随机生成方位角（0到2π均匀分布）
            rotAng1 = rand(1, 1) * 2 * pi;
            
            % 第二步：根据偏振相位函数进行散射角逆采样
            % 考虑入射Stokes向量的偏振状态对散射概率的影响
            % sin(θ)甲醛修正:CDF积分需包含立体角权重，否则大粒子前向衍射峰主导采样。
            qRot = SVIn(2) * cos(2 * rotAng1) + SVIn(3) * sin(2 * rotAng1);  %坐标转换到散射主平面后的线偏振分量
            if strcmp(angleSampler, 'inverse')
                yPhase = S11(scaClLayer, :) * SVIn(1) + S12(scaClLayer, :) * qRot;
                scaAng = inverseSampling(thetaGrid, ...
                    yPhase .* sin(thetaGrid * pi / 180), ...
                    1) / 180 * pi;
            else
                scaAng = sampleScatteringAngle(thetaGrid, ...
                    cdfS11(scaClLayer, :), ...
                    cdfS12(scaClLayer, :), ...
                    SVIn(1), qRot) / 180 * pi;
            end
            % 记录层内散射计数（随机散射）
            layerScatterCount(scaClLayer) = layerScatterCount(scaClLayer) + 1;

            % ========== 计算散射后的Stokes向量 ==========
            % 根据散射角查找对应的散射矩阵元素索引
            [~, iAng] = min(abs(cloudScatMatrix.theta - scaAng / pi * 180));
            
            % 计算散射后的新传播方向
            pScaDir = updateDir(pDir, rotAng1, scaAng);
            
            % ========== Stokes向量的偏振变换 ==========
            % 第一步：旋转到散射平面坐标系
            SVTemp1 = rotSVPhi(SVIn, rotAng1);
            
            % 第二步：应用Mueller散射矩阵
            SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                       S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                       0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                       0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                      SVTemp1 / S11(scaClLayer, iAng);
            
            % ========== 计算子午面旋转角度 ==========
            % 用于将Stokes向量从散射平面旋转到新的传播方向参考系
            temp(1) = sqrt((1 - cos(scaAng)^2) * (1 - pScaDir(3)^2));
            if abs(temp(1)) < 1e-8
                % 特殊情况：散射角接近0或π时的处理
                temp(2) = pi / 2;
            else
                % 一般情况：计算子午面旋转角
                % 考虑方位角的象限和散射几何
                temp(3) = (double((rotAng1 > pi) && (rotAng1 < 2*pi)) - 0.5) * 2 * (pScaDir(3) * cos(scaAng) - pDir(3)) / temp(1);
                
                % 限制余弦值在[-1, 1]范围内，避免数值误差
                temp(3) = double(temp(3) <= -1) * -1 + ...
                                   double(temp(3) >= 1) * 1 + ...
                                   double((temp(3) > -1) && (temp(3) < 1)) * temp(3);
                temp(2) = acos(temp(3));
            end

            % ========== 更新光子状态 ==========
            % 第三步：旋转到新的子午面参考系
            SVIn = rotSVPhi(SVTemp2, -temp(2));   
            
            % 更新散射计数和传播方向
            nScas(iPhoton) = nScas(iPhoton) + 1;
            pDir = pScaDir;
        else
            % ========== 光子离开云层 ==========
            % 光子已经离开云层边界，结束传播
            isAlive = false;

            % 调试信息输出
            if p.Results.debug
                fprintf('Finished %6.2f%%: Photon %d is dead\n', iPhoton/ nPhotons * 100, iPhoton);
            end
        end

    end
end

%% ========== 额外诊断统计（调试用） ==========
% 各层返回I总和，便于分层能量占比分析
diag.layerI = squeeze(sum(diag.SVReturnLayer(1, :, :), 2));
% 散射次数直方图（边界与计数），便于观察散射阶分布
diag.nScasHistEdges = 0:(max(nScas) + 1);
diag.nScasHist = histcounts(nScas, diag.nScasHistEdges);
% 分层散射事件计数
diag.layerScatterCount = layerScatterCount;

%% ========== 结果显示 ==========
% 根据用户设置决定是否显示模拟结果
if strcmp(p.Results.visible, 'on')
    % 创建白色背景的图形窗口
    figure('Color', 'w');

    % ========== 子图1：散射次数分布直方图 ==========
    subplot(311);
    histogram(nScas);
    title('光子散射次数分布', 'FontSize', 12);
    xlabel('散射次数', 'FontSize', 10);
    ylabel('光子数量', 'FontSize', 10);

    % ========== 子图2：Stokes向量分量随距离变化 ==========
    subplot(312);
    hold on;
    
    % 计算距离轴（单位：km）
    distanceAxis = (1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15;
    
    % 绘制四个Stokes向量分量
    plot(distanceAxis, SVReturn(1, :), 'DisplayName', 'I (总强度)');
    plot(distanceAxis, SVReturn(2, :), 'DisplayName', 'Q (线偏振)');
    plot(distanceAxis, SVReturn(3, :), 'DisplayName', 'U (线偏振45°)');
    plot(distanceAxis, SVReturn(4, :), 'DisplayName', 'V (圆偏振)');
    
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
    set(gca, 'YScale', 'log');  % 使用对数坐标显示强度
    title('Stokes向量分量', 'FontSize', 12);
    xlabel('距离 (km)', 'FontSize', 10);
    ylabel('信号强度 (对数尺度)', 'FontSize', 10);
    legend('show');

    % ========== 子图3：退偏振比随距离变化 ==========
    subplot(313);
    
    % 计算线偏振退偏振比：δ = (I - Q) / (I + Q)
    % 这是激光雷达中常用的偏振参数，用于区分球形和非球形粒子
    depRatio = (SVReturn(1, :) - SVReturn(2, :)) ./ (SVReturn(1, :) + SVReturn(2, :));
    
    % 备选计算方法（圆偏振退偏振比）：
    % depRatio = (SVReturn(1, :) + SVReturn(4, :)) ./ (SVReturn(1, :) - SVReturn(4, :));
    
    plot(distanceAxis, depRatio, 'LineWidth', 1.5);
    
    ylim([0, 1]);  % 退偏振比范围：0-1
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
    title('线偏振退偏振比', 'FontSize', 12);
    xlabel('距离 (km)', 'FontSize', 10);
    ylabel('退偏振比', 'FontSize', 10);
    grid on;
end

end