clc;
% close all;
global MS_ENVS;

% Parameter Definition
height = 1:6000;
d = 100e-3;   % distance between laser beam and optical axis of the telescope. (mm)
dTel = 150e-3;   % diameter of the telescope. (mm)
FOV_tel = 0.56e-3 / pi * 180;   % field of view of the telescope. (degree)
iWL = 3;   % 1-> 355; 2-> 532; 3-> 1064;
timeRange = [datenum(2019, 11, 18, 06, 15, 0), datenum(2019, 11, 18, 06, 20, 0)];

%% Read Data
rWD = [0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.250,3.750,4.250,4.750,5.500,6.500,7.500,8.500,9.500,11.00,13.00,15.00,17.00,19.00,21.50,24.50] / 2 * 1e3;   % radius of water droplet measured by OTT. (um)
OTTFile = 'OTT-20191118.h5';
WDDist = h5read(fullfile(MS_ENVS.RootPath, 'data', OTTFile), '/N_Di');
time = datenum(2019, 11, 18, 0, 0:1439, 0);
sizeBin = rWD(1:end);
isInTRange = (time >= timeRange(1)) & (time <= timeRange(2));
nWD = nanmean(WDDist(isInTRange, 1:end), 1);

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
scaCS = NaN(1, length(scaMatrix.theta));
isValidNR = ~isnan(nr);
for iTheta = 1:length(scaMatrix.theta)
    P11(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P11(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
    P12(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P12(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
    P33(iTheta) = trapz(scaMatrix.r(isValidNR), scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2 .* scaMatrix.P33(iWL, isValidNR, iTheta)) / trapz(scaMatrix.Qsca(iWL, isValidNR, iTheta) .* nr(isValidNR) .* pi .* scaMatrix.r(isValidNR).^2);
    scaCS(iTheta) = trapz(scaMatrix.r(isValidNR) * 1e-6, nr(isValidNR) .* scaMatrix.Qsca(iWL, isValidNR, iTheta) .* pi .* (scaMatrix.r(isValidNR) * 1e-6).^2 .* scaMatrix.P11(iWL, isValidNR, iTheta) * (scaMatrix.theta(1) - scaMatrix.theta(2)) / 180 * pi);
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
ylim([0, 4000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

% polarimetry along height
depol_mean = NaN(1, length(height));
bscAll = NaN(1, length(height));
extAll = NaN(1, length(height));
for iH = 1:length(height)
    FOVMin = atan((d - dTel / 2) / height(iH)) / pi * 180;
    FOVMax = atan((d + dTel / 2) / height(iH)) / pi * 180;
    isInFOV = (scaMatrix.theta > (180 - FOVMax)) & (scaMatrix.theta < (180 - FOVMin)) & (scaMatrix.theta > (180 - FOV_tel));
    
    if sum(isInFOV) < 2
        continue;
    end

    depol_mean(1, iH) = nansum(trapz(scaMatrix.theta(isInFOV), P_c(isInFOV), 2)) / nansum(trapz(scaMatrix.theta(isInFOV), P_p(isInFOV), 2));
    bscAll(iH) = nanmean(scaCS(isInFOV));
end

% Rayleigh Scattering
[temp, pres] = read_meteordata(datenum(0, 1, 0), height, 'meteor_data', 'standard_atmosphere');
[mBsc, mExt] = rayleigh_scattering(scaMatrix.inWL(iWL), pres, temp + 273.14, 380, 80);
mdr = 0.009;
vdr = (bscAll .* (depol_mean) ./ (1 + depol_mean) + mBsc .* mdr / (1 + mdr)) ./ (bscAll ./ (1 + depol_mean) + mBsc ./ (1 + mdr));

% ensemble effects
figure('Position', [0, 0, 600, 400], 'Units', 'Pixels', 'Color', 'w');

subplot(121);
hold on;

p1 = plot(depol_mean, height', 'Color', 'k', 'DisplayName', '\delta_{p}');
p2 = plot(vdr, height', 'Color', 'b', 'DisplayName', '\delta_{v}');

hold off;

ylabel('Height (m)');
xlabel('\delta_v');

ylim([0, 4000]);
xlim([0, 0.6]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

legend([p1, p2], 'Location', 'NorthEast');

subplot(122);
hold on;

p1 = plot(bscAll, height, 'Color', 'k', 'DisplayName', 'drizzle');
p2 = plot(mBsc, height, 'Color', 'b', 'DisplayName', 'mol.');

hold off;

ylabel('Height (m)');
xlabel('backscatter (m-1sr-1)');

ylim([0, 4000]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'TickLen', [0.02, 0.02]);

legend([p1, p2], 'Location', 'NorthEast');
