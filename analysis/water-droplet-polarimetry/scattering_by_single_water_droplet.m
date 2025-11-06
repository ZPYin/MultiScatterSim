clc;
close all;
global MS_ENVS;

%% Load Scattering Matrix
scaMatrix = load(fullfile(MS_ENVS.RootPath, 'data', 'water-droplet-scattering-4000.mat'));
dTel = 200e-3;   % diameter of the telescope. (mm)
FOV_tel = 1e-3 / pi * 180;   % field of view of the telescope. (degree)

%% Parameter Definition
r = scaMatrix.r;   % radius of water droplet. [um]
theta = scaMatrix.theta;
inWL = scaMatrix.inWL;

%% Visualization
figure('Position', [0, 0, 400, 500], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.15, 0.1, 0.77, 0.82], 2, 1, 0, 0.1);

% phase function
subplot('Position', subfig(1, :), 'Units', 'Normalized');

lineInstance = [];
lineColor = colormap('jet');

hold on;

for iLine = 1:length(r)
    p = plot(180 - theta, scaMatrix.P11(iLine, :), 'Color', lineColor(ceil(size(lineColor, 1) / (length(r) + 2) * iLine), :), 'DisplayName', sprintf('r=%4d \\mum', r(iLine)));
    lineInstance = cat(2, lineInstance, p);
end

hold off;

xlabel('Scattering Angle (180 - \theta\circ)');
ylabel('Phase Func. (a.u.)');

xlim([180 - max(theta), 180 - min(theta)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'YScale', 'log');

legend(lineInstance, 'Location', 'best');

% polarizing
subplot('Position', subfig(2, :), 'Units', 'Normalized');

lineInstance = [];
lineColor = colormap('jet');

hold on;

for iLine = 1:length(r)
    p = plot(180 - theta, scaMatrix.P_c(iLine, :) ./ scaMatrix.P_p(iLine, :), 'Color', lineColor(ceil(size(lineColor, 1) / (length(r) + 2) * iLine), :), 'DisplayName', sprintf('r=%4d \\mum', r(iLine)));
    lineInstance = cat(2, lineInstance, p);
end

hold off;

xlabel('Scattering Angle (180 - \theta\circ)');
ylabel('\delta_v');

xlim([180 - max(theta), 180 - min(theta)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

% receiving field of view along height
figure('Position', [0, 0, 350, 400], 'Units', 'Pixels', 'Color', 'w');

height = 1:2000;
FOV = 0.2 ./ height / pi * 180;

plot(FOV, height, '-k'); hold on;

xlabel('FOV (\circ)');
ylabel('Height (m)');
title('Receiving FOV along height');

xlim([0, 0.5]);
ylim([0, 2000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

% polarimetry along height
depol = NaN(length(height), length(r));
depol_mean = NaN(1, length(height));
for iH = 1:length(height)
    for iR = 1:length(r)
        FOVLim = dTel / height(iH) / pi * 180 / 2;
        isInFOV = (theta > (180 - FOVLim)) & (theta > (180 - FOV_tel));
        depol(iH, iR) = trapz(theta(isInFOV), scaMatrix.P_c(iR, isInFOV)) / trapz(theta(isInFOV), scaMatrix.P_p(iR, isInFOV));
    end

    depol_mean(1, iH) = nansum(trapz(theta(isInFOV), scaMatrix.P_c(:, isInFOV), 2)) / nansum(trapz(theta(isInFOV), scaMatrix.P_p(:, isInFOV), 2));
end

figure('Position', [0, 0, 350, 400], 'Units', 'Pixels', 'Color', 'w');

p = pcolor(r, height, depol); hold on;
p.EdgeColor = 'None';

ylabel('Height (m)');
xlabel('radius (\mum)');
title('Polarimetry of single water droplet');

ylim([0, 2000]);
xlim([min(r), max(r)]);
caxis([0, 0.4]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickDir', 'out');

colorbar();

% ensemble effects
figure('Position', [0, 0, 350, 400], 'Units', 'Pixels', 'Color', 'w');

hold on;

p = plot(depol_mean, height', 'Color', 'k');

hold off;

ylabel('Height (m)');
xlabel('\delta_v');

ylim([0, 2000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);