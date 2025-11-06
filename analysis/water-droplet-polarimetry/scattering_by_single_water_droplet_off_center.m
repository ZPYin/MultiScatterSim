clc;
close all;
global MS_ENVS;

% Parameter Definition
d = 200e-3;   % distance between laser beam and optical axis of the telescope. (mm)
dTel = 200e-3;   % diameter of the telescope. (mm)
FOV_tel = 1e-3 / pi * 180;   % field of view of the telescope. (degree)
height = 1:2000;
sizeBin = [300, 1000];
nWD = [1, 1];
iWL = 1;   % 1-> 355; 2-> 532; 3-> 1064;

%% Load Scattering Matrix
scaMatrix = load(fullfile(MS_ENVS.RootPath, 'data', 'raindrop-scattering-matrix.mat'));
scaMatrix.P11 = 1 / 2 * (abs(scaMatrix.S1).^2 + abs(scaMatrix.S2).^2);
scaMatrix.P12 = 1 / 2 * (abs(scaMatrix.S2).^2 - abs(scaMatrix.S1).^2);
scaMatrix.P33 = 1 / 2 * (conj(scaMatrix.S2) .* scaMatrix.S1 + scaMatrix.S2 .* conj(scaMatrix.S1));

%% Calculation

% phase function
P11 = NaN(1, length(scaMatrix.theta));
P12 = NaN(1, length(scaMatrix.theta));
P33 = NaN(1, length(scaMatrix.theta));
nr = interp1(sizeBin, nWD, scaMatrix.r);
isValidNR = ~isnan(nr);
for iTheta = 1:length(scaMatrix.theta)
    P11(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P11(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
    P12(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P12(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
    P33(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P33(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
end

L1 = (P11 + P12);
L2 = (P11 - P12);
L3 = P33;
u = cos(scaMatrix.theta / 180 * pi);
P_p = 1 / 8 * (3 * L1 + 3 * L2 .* u.^2 + 2 * L3 .* u);
P_c = 1 / 8 * (L1 + L2 .* u.^2 - 2 * L3 .* u);

%% Visualization
figure('Position', [0, 0, 400, 500], 'Units', 'Pixels', 'Color', 'w');

subfig = subfigPos([0.15, 0.1, 0.77, 0.82], 2, 1, 0, 0.1);

% phase function
subplot('Position', subfig(1, :), 'Units', 'Normalized');

hold on;
p = plot(180 - scaMatrix.theta, P11, 'Color', 'k');

hold off;

xlabel('Scattering Angle (180 - \theta\circ)');
ylabel('Phase Func. (a.u.)');

xlim([180 - max(scaMatrix.theta), 180 - min(scaMatrix.theta)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02], 'YScale', 'log');

% polarizing
subplot('Position', subfig(2, :), 'Units', 'Normalized');

hold on;
plot(180 - scaMatrix.theta, P_c ./ P_p, 'Color', 'k');

hold off;

xlabel('Scattering Angle (180 - \theta\circ)');
ylabel('\delta_v');

xlim([180 - max(scaMatrix.theta), 180 - min(scaMatrix.theta)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

% receiving field of view along height
figure('Position', [0, 0, 350, 400], 'Units', 'Pixels', 'Color', 'w');

FOVMin = atan((d - dTel / 2) ./ height) / pi * 180;
FOVMax = atan((d + dTel / 2) ./ height) / pi * 180;

patch([FOVMin, fliplr(FOVMax)], [height, fliplr(height)], 'k'); hold on;

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
    FOVMin = atan((d - dTel / 2) / height(iH)) / pi * 180;
    FOVMax = atan((d + dTel / 2) / height(iH)) / pi * 180;
    isInFOV = (scaMatrix.theta > (180 - FOVMax)) & (scaMatrix.theta < (180 - FOVMin)) & (scaMatrix.theta > (180 - FOV_tel));
    
    if sum(isInFOV) < 2
        continue;
    end

    depol_mean(1, iH) = nansum(trapz(scaMatrix.theta(isInFOV), P_c(isInFOV), 2)) / nansum(trapz(scaMatrix.theta(isInFOV), P_p(isInFOV), 2));
end

% ensemble effects
figure('Position', [0, 0, 350, 400], 'Units', 'Pixels', 'Color', 'w');

hold on;

p = plot(depol_mean, height', 'Color', 'k');

hold off;

ylabel('Height (m)');
xlabel('\delta_v');

ylim([0, 2000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

