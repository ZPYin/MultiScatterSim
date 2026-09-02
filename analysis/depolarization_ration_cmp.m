function depolarization_ration_cmp()
% 单层云线偏振退偏振比示例
% 固定有效半径 reff=6 um，对比不同消光系数下的退偏振比随深度变化。

clear; clc; close all;



%% 激光雷达与云层参数
lidarCfg.lambda = 532;             % 波长 [nm]
lidarCfg.rBeam = 5.0;              % 光束半径 [mm]
lidarCfg.divBeam = 0.05;           % 发散角 [mrad]
lidarCfg.distBeamTel = 0;          % 发射/接收光轴间距 [mm]
lidarCfg.FOV = 6;                  % 接收视场角 [mrad]
lidarCfg.rTel = 50.0;              % 望远镜半径 [mm]
lidarCfg.binWidth = 10;            % 距离门宽度 [ns]
lidarCfg.nBins = 200;              % 距离门数量
lidarCfg.SBeam = [1; 1; 0; 0];     % 入射线偏振 Stokes 向量

cloudCfg.clBase = 1000;            % 云底高度 [m]
cloudCfg.clDh = 200;               % 云层厚度 [m]
cloudCfg.nClLayers = 1;
cloudCfg.kappa = 9;
cloudCfg.clReff = 6;               % 有效半径 [um]

extCoeffs = [0.01, 0.03, 0.05];    % 消光系数 [m^-1]
nPhotons = 1e4;                    % 示例用光子数；正式结果可增大到1e6
depth = (1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15;

%% 蒙特卡洛仿真与退偏振比计算
figure('Color', 'w', 'Position', [200, 200, 720, 460]);
hold on;

lineStyles = {'k-', 'r--', 'b:'};
for iExt = 1:numel(extCoeffs)
    cloudCfg.clExt = extCoeffs(iExt);
    fprintf('Running ext = %.2f m^-1 ...\n', cloudCfg.clExt);

    SVReturn = MSMCSim(lidarCfg, cloudCfg, ...
        'nPhoton', nPhotons, ...
        'visible', 'off', ...
        'cloud_distribution', 'modified-gamma');

    I = movmean(SVReturn(1, :), 5);
    Q = movmean(SVReturn(2, :), 5);
    delta = (I - Q) ./ (I + Q);     % delta = P_perp / P_parallel

    valid = isfinite(delta) & delta >= 0 & I > max(I) * 1e-5;
    plot(depth(valid), delta(valid), lineStyles{iExt}, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('\\sigma_{ext}=%.2f m^{-1}', cloudCfg.clExt));
end

xlabel('Depth / m');
ylabel('Linear depolarization ratio');
title('Single-layer cloud depolarization ratio, r_{eff}=6 \mum');
legend('Location', 'northeast');
xlim([0, cloudCfg.clDh]);
ylim([0, 1]);
grid on;
box on;

end
