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
addParameter(p, 'preparedOptics', [], @(x) isempty(x) || isstruct(x));
addParameter(p,'forcedReceiverRotation','vectorbasis',@ischar);
addParameter(p,'forcedReceiverPathMode','allangles',@ischar);
addParameter(p,'propagationWeightNormalization','p11',@ischar);
addParameter(p,'proposalBackwardMix',.2,@isnumeric);
addParameter(p,'proposalTailDegrees',170,@isnumeric);
addParameter(p,'geometryMode','native',@ischar);
addParameter(p,'exactTimeWindowStop',false,@islogical);
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
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'MC_Optics'));


if isempty(p.Results.preparedOptics)
    unmatchedArgs = p.Unmatched;
    extraArgs = reshape([fieldnames(unmatchedArgs).'; struct2cell(unmatchedArgs).'], 1, []);
    optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
        'scatMatrixFile', p.Results.scatMatrixFile, extraArgs{:});
    diag.opticsSource = 'computed';
else
    optics = p.Results.preparedOptics;
    validatePreparedOptics(optics, lidarCfg, cloudCfg);
    diag.opticsSource = 'prepared';
end

weightMode=validatestring(lower(p.Results.propagationWeightNormalization),{'p11','legacy','samplingpdf'});
receiverRotation=validatestring(lower(p.Results.forcedReceiverRotation),{'vectorbasis','legacy'});
receiverPaths=validatestring(lower(p.Results.forcedReceiverPathMode),{'allangles','reciprocity'});
geometryMode=validatestring(lower(p.Results.geometryMode),{'native','mscart'});
if strcmp(weightMode,'p11')
    optics=prepScaSamp(optics,p.Results.proposalBackwardMix,p.Results.proposalTailDegrees);
end
diag.propagationWeightNormalization=weightMode;
diag.forcedReceiverRotation=receiverRotation;
diag.forcedReceiverPathMode=receiverPaths;
diag.geometryMode=geometryMode;
diag.modelRevision='mscart_validated_20260929';
cloudScatMatrix = optics.cloudScatMatrix;

% 提取关键光学参数
muExt = optics.muExt;     % 消光系数 [m^-1]
muSca = optics.muSca;     % 散射系数 [m^-1]
albedo = optics.albedo;

% Publish the exact optical properties used by photon transport so that
% full-profile composition uses the same cloud optical depth.
diag.muExt = muExt;
diag.muSca = muSca;
diag.albedo = albedo;

% 提取Mueller散射矩阵元素（用于偏振计算）
S11 = optics.S11;
S12 = optics.S12;        % 线偏振相关的散射矩阵元素
S33 = optics.S33;
S34 = optics.S34;



thetaGrid = optics.theta;
cdfS11 = optics.cdfS11;
cdfS12 = optics.cdfS12;

% ========== 云层边界设置 ==========
% 根据Winker & Poole (1995)的建议，云边界对应光学厚度为2

clBoundX = 4.0 ./ muExt;
clBoundY = 4.0 ./ muExt;
if strcmp(geometryMode,'mscart'),clBoundX(:)=inf;clBoundY(:)=inf;end





rCl = lidarCfg.rBeam * 1e-3 + lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase;   % 云底激光束半径 [m]


SVReturn = zeros(4, lidarCfg.nBins);   % 返回信号的Stokes向量 [4×nBins]
diag.SVReturnLayer = zeros(4, lidarCfg.nBins, cloudCfg.nClLayers); % 分层返回累积
len = zeros(1, nPhotons);
nScas = zeros(1, nPhotons);
layerScatterCount = zeros(1, cloudCfg.nClLayers); % 分层散射事件计数

% 散射阶数诊断: 按阶数累积返回信号
maxOrderEstimate = 50;
diag.SVByOrder = zeros(4, lidarCfg.nBins, maxOrderEstimate);
diag.orderPhotonCount = zeros(1, maxOrderEstimate);

for iPhoton = 1:nPhotons

    pPos = [0; 0; 0];                  % 光子位置向量 [x; y; z] [m]
    temp = zeros(1, 5);                % 临时计算变量


    isAlive = true;                    % 光子存活标志
    

    % 在云底完整圆形光束截面内均匀抽样：r=R*sqrt(U)，phi∈[0, 2*pi)
    launchRadius = rCl * sqrt(rand(1, 1));
    launchAzimuth = 2 * pi * rand(1, 1);
    pPos(1) = launchRadius * cos(launchAzimuth);
    pPos(2) = launchRadius * sin(launchAzimuth);
    pPos(3) = cloudCfg.clBase;         % Z方向设为云底高度
    

    pDir = pPos / norm(pPos - [0; 0; 0], 2);
    
    % 根据初始传播方向旋转入射Stokes向量

    SVIn = rotSVPhi(lidarCfg.SBeam, atan2(pDir(2), pDir(1)));

    % ========== 光子传播循环 ==========
    while (isAlive)

        % ========== 自由传播步骤 ==========



        [pPos, isInMedium, ds, scaClLayer] = freeRunning(pPos, pDir, muExt, cloudCfg.clDh, cloudCfg.clBase, clBoundX, clBoundY);
        len(iPhoton) = len(iPhoton) + ds;   % 累积传播距离


        if p.Results.exactTimeWindowStop && len(iPhoton)>=3e8*lidarCfg.binWidth*1e-9*lidarCfg.nBins
            isAlive=false;continue;
        end
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

            fovRadius=(lidarCfg.FOV*1e-3*pPos(3)+lidarCfg.rTel*1e-3*2)/2;
            if strcmp(geometryMode,'mscart'),fovRadius=tan(lidarCfg.FOV*1e-3/2)*pPos(3);end
            if sqrt((pPos(1)-lidarCfg.distBeamTel*1e-3)^2+pPos(2)^2)<fovRadius
                


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

                if strcmp(receiverPaths,'reciprocity') && cloudScatMatrix.theta(iAng)<90

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
                receiverOutRotation=-rotAng2;
                if strcmp(receiverRotation,'vectorbasis')
                    [rotAng1,receiverOutRotation]=recRotAng(pDir.',pScaDir.',pDir.');
                end
                SVTemp1 = rotSVPhi(SVIn, rotAng1);
                
                % 第二步：应用Mueller散射矩阵进行偏振变换
                SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                                       S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                                       0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                                       0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                                    SVTemp1 / S11(scaClLayer, iAng);




                timeBin = (len(iPhoton) + (pPos(3) - cloudCfg.clBase) / max(abs(pScaDir(3)), eps)) / 3e8 / (lidarCfg.binWidth * 1e-9);
                if strcmp(geometryMode,'mscart')
                    receiverDistance=norm(pPos-[lidarCfg.distBeamTel*1e-3;0;0]);
                    timeBin=(len(iPhoton)+receiverDistance-cloudCfg.clBase)/3e8/(lidarCfg.binWidth*1e-9);
                end
                iBin = ceil(timeBin);
                
                % ========== 累积返回信号 ==========
                if isfinite(iBin) && iBin >= 1 && iBin <= lidarCfg.nBins


                    if strcmp(receiverRotation,'vectorbasis')
                        receiverSV=rotSVPhi(SVTemp2,receiverOutRotation);
                    else
                        receiverSV=rotSVPhi(rotSVPhi(SVTemp2,-rotAng2),atan2(pScaDir(2),pScaDir(1)));
                    end
                    pathFactor=2^(nScas(iPhoton)-1);
                    if strcmp(receiverPaths,'allangles'),pathFactor=.5;end
                    SV_final=receiverSV*probEnter*pathFactor;
                    SVReturn(:, iBin) = SVReturn(:, iBin) + SV_final;
                    diag.SVReturnLayer(:, iBin, scaClLayer) = diag.SVReturnLayer(:, iBin, scaClLayer) + SV_final;

                    % 散射阶数诊断: 接收阶数 = 已完成随机散射数 + 本次强制散射
                    order = nScas(iPhoton) + 1;
                    if order <= maxOrderEstimate
                        diag.SVByOrder(:, iBin, order) = diag.SVByOrder(:, iBin, order) + SV_final;
                        diag.orderPhotonCount(order) = diag.orderPhotonCount(order) + 1;
                    end
                end

                % ========== 更新剩余光子强度 ==========

                SVIn = SVIn * (1 - probEnter);

            end

            % ========== 随机散射事件处理 ==========

            
            % ========== 随机采样散射角度 ==========

            rotAng1 = rand(1, 1) * 2 * pi;
            



            qRot = SVIn(2) * cos(2 * rotAng1) + SVIn(3) * sin(2 * rotAng1);
            if strcmp(weightMode,'p11')
                % Exact inverse CDF of the discrete scalar proposal in either sampler mode.
                scaAng=sampleScatteringAngle(thetaGrid,cdfS11(scaClLayer,:),cdfS12(scaClLayer,:),SVIn(1),qRot)/180*pi;
            elseif strcmp(angleSampler, 'inverse')
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
            

            if strcmp(weightMode,'p11')
                SVTemp2=SVTemp2*optics.proposalCorrection(scaClLayer,iAng);
            elseif strcmp(weightMode,'samplingpdf')
                density=S11(scaClLayer,iAng)*SVTemp1(1)+S12(scaClLayer,iAng)*SVTemp1(2);
                total=SVTemp1(1)*cdfS11(scaClLayer,end)+SVTemp1(2)*cdfS12(scaClLayer,end);
                denominator=max(density*cdfS11(scaClLayer,end)/max(total,eps),eps);
                SVTemp2=SVTemp2*S11(scaClLayer,iAng)/denominator;
            end
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

% ========== 散射阶数分解的偏振诊断 ==========
diag.maxOrderEstimate = maxOrderEstimate;
I_by_order = reshape(sum(diag.SVByOrder(1, :, :), 2), [], 1);  % [maxOrder x 1]
Q_by_order = reshape(sum(diag.SVByOrder(2, :, :), 2), [], 1);
diag.I_by_order = I_by_order;
diag.Q_by_order = Q_by_order;
diag.Pperp_by_order = (I_by_order - Q_by_order) / 2;
diag.Ppar_by_order = (I_by_order + Q_by_order) / 2;
diag.delta_by_order = diag.Pperp_by_order ./ max(diag.Ppar_by_order, eps);
diag.f_by_order = I_by_order ./ max(sum(I_by_order), eps);
diag.lastActiveOrder = find(I_by_order > 0, 1, 'last');
if isempty(diag.lastActiveOrder), diag.lastActiveOrder = 0; end

%% ========== Result display ==========
if strcmp(p.Results.visible, 'on')
    % CPU and GPU use the same diagnostic quantities and axis conventions.
    counts = double(diag.nScasHist(:).');
    lastNonzero = find(counts > 0, 1, 'last');
    if isempty(lastNonzero)
        counts = 0;
    else
        counts = counts(1:lastNonzero);
    end
    orders = 0:(numel(counts) - 1);

    distanceAxis = (1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15;
    IPlot = abs(double(SVReturn(1, :)));
    QPlot = abs(double(SVReturn(2, :)));
    IPlot(IPlot == 0) = NaN;
    QPlot(QPlot == 0) = NaN;

    denominator = double(SVReturn(1, :) + SVReturn(2, :));
    signalFloor = max(max(abs(denominator)) * 1e-12, eps);
    validDelta = isfinite(denominator) & denominator > signalFloor;
    delta = NaN(1, lidarCfg.nBins);
    delta(validDelta) = double(SVReturn(1, validDelta) - ...
        SVReturn(2, validDelta)) ./ denominator(validDelta);

    figure('Color', 'w', 'Name', 'MSMC diagnostic results', ...
        'NumberTitle', 'off');

    subplot(3, 1, 1);
    bar(orders, counts, 1, 'FaceColor', [0.25, 0.60, 0.82], ...
        'EdgeColor', [0.15, 0.15, 0.15]);
    xlim([-0.5, max(10, orders(end) + 1)]);
    title('Photon scattering order distribution', 'FontSize', 12);
    xlabel('Scattering order', 'FontSize', 10);
    ylabel('Photon count', 'FontSize', 10);
    grid on;

    subplot(3, 1, 2);
    semilogy(distanceAxis, IPlot, 'LineWidth', 1.0, ...
        'DisplayName', 'I');
    hold on;
    semilogy(distanceAxis, QPlot, 'LineWidth', 1.0, ...
        'DisplayName', '|Q|');
    hold off;
    xlim([0, distanceAxis(end)]);
    title('Returned Stokes components', 'FontSize', 12);
    xlabel('Distance (m)', 'FontSize', 10);
    ylabel('Signal intensity', 'FontSize', 10);
    legend('show', 'Location', 'best');
    grid on;

    subplot(3, 1, 3);
    plot(distanceAxis, delta, 'LineWidth', 1.5);
    xlim([0, distanceAxis(end)]);
    ylim([0, 1]);
    title('Linear depolarization ratio', 'FontSize', 12);
    xlabel('Distance (m)', 'FontSize', 10);
    ylabel('Depolarization ratio', 'FontSize', 10);
    grid on;
end

end

function validatePreparedOptics(optics, lidarCfg, cloudCfg)
required = {'cloudScatMatrix', 'muExt', 'muSca', 'albedo', ...
    'S11', 'S12', 'S33', 'S34', 'theta', 'cdfS11', 'cdfS12'};
for iField = 1:numel(required)
    if ~isfield(optics, required{iField})
        error('MSMCSim:InvalidPreparedOptics', ...
            'preparedOptics is missing field "%s".', required{iField});
    end
end
if isfield(optics, 'lambda') && optics.lambda ~= lidarCfg.lambda
    error('MSMCSim:PreparedOpticsWavelength', ...
        'preparedOptics wavelength (%g nm) does not match lidarCfg.lambda (%g nm).', ...
        optics.lambda, lidarCfg.lambda);
end
if numel(optics.muExt) ~= cloudCfg.nClLayers
    error('MSMCSim:PreparedOpticsLayers', ...
        'preparedOptics has %d layers but cloudCfg.nClLayers is %d.', ...
        numel(optics.muExt), cloudCfg.nClLayers);
end
end

