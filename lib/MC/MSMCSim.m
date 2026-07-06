function [SVReturn, diag] = MSMCSim(lidarCfg, cloudCfg, varargin)
% MSMCSIM 基于蒙特卡洛方法的偏振激光雷达多次散射模拟
%
% 该函数追踪光子在云层中的传播、散射、偏振态演化和接收过程，
% 最终统计不同距离门内的 Stokes 回波信号。云滴散射特性由
% cloudScatteringMatrix 根据 Mie 散射数据库和粒径分布计算。
%
% 用法:
%    SVReturn = MSMCSim(lidarCfg, cloudCfg)
%    [SVReturn, diag] = MSMCSim(lidarCfg, cloudCfg, 'nPhoton', 1e5, 'visible', 'off')
%
% lidarCfg 字段:
%    lambda       激光波长 [nm]
%    rBeam        发射光束半径 [mm]
%    divBeam      发射光束发散角 [mrad]
%    FOV          接收视场角 [mrad]
%    rTel         望远镜半径 [mm]
%    binWidth     距离门宽度 [ns]
%    nBins        距离门数量
%    distBeamTel  发射光轴与接收光轴间距 [mm]
%    SBeam        入射 Stokes 向量 [I,Q,U,V]
%
% cloudCfg 字段:
%    clBase       云底高度 [m]
%    clDh         单层云厚度 [m]
%    nClLayers    云层数
%    kappa        修正伽马分布形状参数
%    clReff       有效半径 [um]
%    clExt        消光系数 [m^-1]，旧模式使用
%    N_d          云滴数浓度 [m^-3]，若提供则优先使用 N_d/Reff 模式
%
% 可选参数:
%    nPhoton              模拟光子数，默认 1e4
%    visible              是否显示诊断图，'on' 或 'off'
%    debug                是否输出调试信息
%    useRussianRoulette   是否使用轮盘赌终止低权重光子
%    energyThreshold      轮盘赌权重阈值
%    scatMatrixFile       可选 Mie 散射数据库文件
%    angleSampler         散射角采样方法，'fast' 或 'inverse'
%
% 输出:
%    SVReturn  4 x nBins 的 Stokes 回波信号
%    diag      分层回波、散射次数等诊断信息
%
% References:
%    Winker and Poole (1995), Applied Physics B, 60(4), 341-344
%    Mishchenko et al. (2002), Scattering, Absorption, and Emission of Light by Small Particles

% ========== 输入参数解析 ==========

p = inputParser;
p.KeepUnmatched = true;

% 必需参数验证
addRequired(p, 'lidarCfg', @isstruct);
addRequired(p, 'cloudCfg', @isstruct);   % 云层配置必须为结构体


addParameter(p, 'nPhoton', 1e4, @isnumeric);      % 模拟光子数，默认1万个
addParameter(p, 'visible', 'on', @ischar);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'useRussianRoulette', true, @islogical);      % 是否使用轮盘对赌
addParameter(p, 'energyThreshold', 0.01, @isnumeric);
addParameter(p, 'useLocalAtten', false, @islogical); % 调试：使用本层局部透过而非累计垂直光学深度
addParameter(p, 'scatMatrixFile', '', @ischar);
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



addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'sampling'));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scattering-phase-function'));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cloud-microphysics'));


projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
clScatMatrixFile = fullfile(projectDir, 'tests', 'data', ...
    'water_droplet_scattering_532nm_finer_validation_r_step_0p01um.mat');
if ~isempty(p.Results.scatMatrixFile)
    clScatMatrixFile = p.Results.scatMatrixFile;
end

    %'water_droplet_scattering_532nm_finer_validation_r_step_0p01um.mat');  % 可选：高精度数据库


% 基于Mie散射理论，考虑云滴粒径分布（修正伽马分布）
dist = 'modified-gamma';
unmatchedArgs = p.Unmatched;
if isfield(unmatchedArgs, 'cloud_distribution')
    dist = unmatchedArgs.cloud_distribution;
    unmatchedArgs = rmfield(unmatchedArgs, 'cloud_distribution');
end
extraArgs = reshape([fieldnames(unmatchedArgs).'; struct2cell(unmatchedArgs).'], 1, []);


cloudArgs = {
    'cloud_distribution', dist, ...
    'kappa', cloudCfg.kappa, ...
    'R_eff', cloudCfg.clReff, ...
    'b', 1/9};
cloudArgs = [cloudArgs, extraArgs];

if isfield(cloudCfg, 'N_d') && ~isempty(cloudCfg.N_d)
    cloudArgs = [cloudArgs, {'N_d', cloudCfg.N_d}];
else
    cloudArgs = [cloudArgs, {'sigma_e', cloudCfg.clExt}];
end

cloudScatMatrix = cloudScatteringMatrix(lidarCfg.lambda, clScatMatrixFile, cloudArgs{:});

% 提取关键光学参数
muExt = cloudScatMatrix.muExt;     % 消光系数 [m^-1]
muSca = cloudScatMatrix.muSca;     % 散射系数 [m^-1]  
albedo = muSca ./ muExt;

% 提取Mueller散射矩阵元素（用于偏振计算）
S11 = cloudScatMatrix.S11;
S12 = cloudScatMatrix.S12;        % 线偏振相关的散射矩阵元素
S33 = cloudScatMatrix.S33;
S34 = cloudScatMatrix.S34;



thetaGrid = cloudScatMatrix.theta;
thetaWeights = [0, diff(thetaGrid)] .* sin(thetaGrid * pi / 180);
cdfS11 = cumsum(S11 .* thetaWeights, 2);
cdfS12 = cumsum(S12 .* thetaWeights, 2);

% ========== 云层边界设置 ==========
% 根据Winker & Poole (1995)的建议，云边界对应光学厚度为2

clBoundX = 4.0 ./ muExt;
clBoundY = 4.0 ./ muExt;





rCl = lidarCfg.rBeam * 1e-3 + lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase * 1e-3;   % 云底激光束半径 [m]


SVReturn = zeros(4, lidarCfg.nBins);   % 返回信号的Stokes向量 [4×nBins]
diag.SVReturnLayer = zeros(4, lidarCfg.nBins, cloudCfg.nClLayers); % 分层返回累积
len = zeros(1, nPhotons);
nScas = zeros(1, nPhotons);
layerScatterCount = zeros(1, cloudCfg.nClLayers); % 分层散射事件计数

for iPhoton = 1:nPhotons

    pPos = [0; 0; 0];                  % 光子位置向量 [x; y; z] [m]
    temp = zeros(1, 5);                % 临时计算变量


    isAlive = true;                    % 光子存活标志
    

    pPos(1) = rand(1, 1) * rCl;       % X方向随机位置
    pPos(2) = rand(1, 1) * rCl;       % Y方向随机位置
    pPos(3) = cloudCfg.clBase;         % Z方向设为云底高度
    

    while (norm(pPos(1:2), 2) > rCl)
        % 重新生成位置直到满足圆形边界条件
        pPos(1) = rand(1, 1) * rCl;
        pPos(2) = rand(1, 1) * rCl;
    end
    

    pDir = pPos / norm(pPos - [0; 0; 0], 2);
    
    % 根据初始传播方向旋转入射Stokes向量

    SVIn = rotSVPhi(lidarCfg.SBeam, atan2(pDir(2), pDir(1)));

    % ========== 光子传播循环 ==========
    while (isAlive)

        % ========== 自由传播步骤 ==========



        [pPos, isInMedium, ds, scaClLayer] = freeRunning(pPos, pDir, muExt, cloudCfg.clDh, cloudCfg.clBase, clBoundX, clBoundY);
        len(iPhoton) = len(iPhoton) + ds;   % 累积传播距离


        if (isInMedium)
            % 光子仍在云层内，进行吸收判断
            if (rand(1, 1) > albedo(scaClLayer))
                % 光子被吸收，结束传播
                isAlive = false;
                continue;
            end
            
            % ========== 轮盘对赌(Russian Roulette) ==========

            if p.Results.useRussianRoulette && SVIn(1) < p.Results.energyThreshold

                survivalProbability = SVIn(1) / p.Results.energyThreshold;
                
                if rand(1, 1) > survivalProbability

                    if p.Results.debug && mod(iPhoton, 1000) == 0
                        fprintf('光子 %d: 轮盘对赌终止, 当前能量=%.4f\n', iPhoton, SVIn(1));
                    end
                    break;
                else

                    SVIn = SVIn / survivalProbability;
                    if p.Results.debug && mod(iPhoton, 10000) == 0
                        fprintf('Photon %d: Russian roulette survived, weight=%.4f\n', iPhoton, SVIn(1));
                    end
                end
            end



            % FOV条件：光子到望远镜光轴的距离 < FOV半径

            if (sqrt((pPos(1) - lidarCfg.distBeamTel * 1e-3)^2 + pPos(2)^2) < ((lidarCfg.FOV * 1e-3 * pPos(3) + lidarCfg.rTel * 1e-3 * 2) / 2))
                


                pScaDir = [lidarCfg.distBeamTel * 1e-3; 0; 0] - pPos;
                pScaDir = pScaDir / norm(pScaDir, 2);

                % ========== 计算散射角度 ==========
                % 计算散射角和旋转角，用于Stokes向量变换
                [scaAng,rotAng1, rotAng2] = scaRotAng(pDir, pScaDir);


                layerScatterCount(scaClLayer) = layerScatterCount(scaClLayer) + 1;



                [~, iAng] = min(abs(cloudScatMatrix.theta - scaAng / pi * 180));
                
                % 计算垂直光学厚度（考虑多层云累积消光）
                % 1. 下方完整层的光学厚度
                opticalDepthBelow = 0;
                if scaClLayer > 1
                    opticalDepthBelow = sum(muExt(1:scaClLayer-1)) * cloudCfg.clDh;
                end
                

                layerBottom = cloudCfg.clBase + (scaClLayer - 1) * cloudCfg.clDh;
                opticalDepthCurrent = (pPos(3) - layerBottom) * muExt(scaClLayer);
                

                verticalOD = opticalDepthBelow + opticalDepthCurrent;

                attenOD = verticalOD;
                if p.Results.useLocalAtten
                    attenOD = opticalDepthCurrent;
                end

                if cloudScatMatrix.theta(iAng) < 90

                    probEnter = 0;
                else
                    % 后向散射情况，计算详细的进入概率

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




                timeBin = (len(iPhoton) + (pPos(3) - cloudCfg.clBase) / max(abs(pScaDir(3)), eps)) / 3e8 / (lidarCfg.binWidth * 1e-9);
                iBin = ceil(timeBin);
                
                % ========== 累积返回信号 ==========
                if isfinite(iBin) && iBin >= 1 && iBin <= lidarCfg.nBins


                    SV_final = rotSVPhi(rotSVPhi(SVTemp2, -rotAng2), atan2(pScaDir(2), pScaDir(1))) * probEnter* 2^(nScas(iPhoton) - 1);
                    SVReturn(:, iBin) = SVReturn(:, iBin) + SV_final;
                    diag.SVReturnLayer(:, iBin, scaClLayer) = diag.SVReturnLayer(:, iBin, scaClLayer) + SV_final;
                end

                % ========== 更新剩余光子强度 ==========

                SVIn = SVIn * (1 - probEnter);

            end

            % ========== 随机散射事件处理 ==========

            
            % ========== 随机采样散射角度 ==========

            rotAng1 = rand(1, 1) * 2 * pi;
            



            qRot = SVIn(2) * cos(2 * rotAng1) + SVIn(3) * sin(2 * rotAng1);
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
            

            pScaDir = updateDir(pDir, rotAng1, scaAng);
            

            % 第一步：旋转到散射平面坐标系
            SVTemp1 = rotSVPhi(SVIn, rotAng1);
            
            % 第二步：应用Mueller散射矩阵
            SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                                   S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                                   0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                                   0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                            SVTemp1 / S11(scaClLayer, iAng);
            

            % 用于将Stokes向量从散射平面旋转到新的传播方向参考系
            temp(1) = sqrt((1 - cos(scaAng)^2) * (1 - pScaDir(3)^2));
            if abs(temp(1)) < 1e-8
                % 特殊情况：散射角接近0或π时的处理
                temp(2) = pi / 2;
            else
                % 一般情况：计算子午面旋转角

                temp(3) = (double((rotAng1 > pi) && (rotAng1 < 2*pi)) - 0.5) * 2 * (pScaDir(3) * cos(scaAng) - pDir(3)) / temp(1);
                

                temp(3) = double(temp(3) <= -1) * -1 + ...
                                   double(temp(3) >= 1) * 1 + ...
                                   double((temp(3) > -1) && (temp(3) < 1)) * temp(3);
                temp(2) = acos(temp(3));
            end


            % 第三步：旋转到新的子午面参考系
            SVIn = rotSVPhi(SVTemp2, -temp(2));   
            

            nScas(iPhoton) = nScas(iPhoton) + 1;
            pDir = pScaDir;
        else
            % ========== 光子离开云层 ==========

            isAlive = false;

            % 调试信息输出
            if p.Results.debug
                fprintf('Finished %6.2f%%: Photon %d is dead\n', iPhoton/ nPhotons * 100, iPhoton);
            end
        end

    end
end



diag.layerI = squeeze(sum(diag.SVReturnLayer(1, :, :), 2));
% 散射次数直方图（边界与计数），便于观察散射阶分布
diag.nScasHistEdges = 0:(max(nScas) + 1);
diag.nScasHist = histcounts(nScas, diag.nScasHistEdges);
% 分层散射事件计数
diag.layerScatterCount = layerScatterCount;

%% ========== Result display ==========
if strcmp(p.Results.visible, 'on')
    figure('Color', 'w');

    subplot(311);
    histogram(nScas);
    title('Photon scattering order distribution', 'FontSize', 12);
    xlabel('Scattering order', 'FontSize', 10);
    ylabel('Photon count', 'FontSize', 10);

    subplot(312);
    hold on;
    distanceAxis = (1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15;
    plot(distanceAxis, SVReturn(1, :), 'DisplayName', 'I');
    plot(distanceAxis, SVReturn(2, :), 'DisplayName', 'Q');
    plot(distanceAxis, SVReturn(3, :), 'DisplayName', 'U');
    plot(distanceAxis, SVReturn(4, :), 'DisplayName', 'V');
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
    set(gca, 'YScale', 'log');
    title('Returned Stokes components', 'FontSize', 12);
    xlabel('Distance (km)', 'FontSize', 10);
    ylabel('Signal intensity', 'FontSize', 10);
    legend('show');

    subplot(313);
    depRatio = (SVReturn(1, :) - SVReturn(2, :)) ./ (SVReturn(1, :) + SVReturn(2, :));
    plot(distanceAxis, depRatio, 'LineWidth', 1.5);
    ylim([0, 1]);
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
    title('Linear depolarization ratio', 'FontSize', 12);
    xlabel('Distance (km)', 'FontSize', 10);
    ylabel('Depolarization ratio', 'FontSize', 10);
    grid on;
end

end

