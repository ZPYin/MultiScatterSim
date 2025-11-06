clc;
close all;
global MS_ENVS;

% Parameter Definition
height = 1:6000;
d = 200e-3;   % distance between laser beam and optical axis of the telescope. (mm)
dTel = 200e-3;   % diameter of the telescope. (mm)
FOV_tel = 1e-3 / pi * 180 / 2;   % field of view of the telescope. (degree)
iWL = 1;   % 1-> 355; 2-> 532; 3-> 1064;
mdr = 0.01;   % molecule depolarization ratio
% tIdxArr = [953, 962, 990, 991, 999, 1010, 1024, 1035, 1041, 1122];
tIdxArr = 900:5:1200;

%% Read Data
rWD = [0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.250,3.750,4.250,4.750,5.500,6.500,7.500,8.500,9.500,11.00,13.00,15.00,17.00,19.00,21.50,24.50] / 2 * 1e3;   % radius of water droplet measured by OTT. (um)
OTTFile = 'OTT-Fitted-20200602.h5';
WDDist = h5read(fullfile(MS_ENVS.RootPath, 'data', OTTFile), '/N_Di_fit');
time = datenum(OTTFile(5:12), 'yyyymmdd') + datenum(0, 1, 0, 0, 0:1439, 0);
sizeBin = rWD(1:end);
WDDist(:, 1:3) = WDDist(:, 1:3)*10;

%% Load Scattering Matrix
scaMatrix = load(fullfile(MS_ENVS.RootPath, 'data', 'raindrop-scattering-matrix.mat'));
scaMatrix.P11 = 1 / 2 * (abs(scaMatrix.S1).^2 + abs(scaMatrix.S2).^2);
scaMatrix.P12 = 1 / 2 * (abs(scaMatrix.S2).^2 - abs(scaMatrix.S1).^2);
scaMatrix.P33 = 1 / 2 * (conj(scaMatrix.S2) .* scaMatrix.S1 + scaMatrix.S2 .* conj(scaMatrix.S1));

%% Rayleigh Scattering
[temp, pres] = read_meteordata(datenum(0, 1, 0), height, 'meteor_data', 'standard_atmosphere');
[mBsc, mExt] = rayleigh_scattering(scaMatrix.inWL(iWL), pres, temp + 273.14, 380, 80);

%% Calculate VDR and rain amount
rainAmountFullArr = NaN(1, length(tIdxArr));
rainAmount1mmArr = NaN(1, length(tIdxArr));
vdrArr = NaN(1, length(tIdxArr));
secOrMmtFullArr = NaN(1, length(tIdxArr));
secOrMmt1mmArr = NaN(1, length(tIdxArr));
depol_mean = NaN(1, length(height));
bscAll = NaN(1, length(height));
extAll = NaN(1, length(height));
for iT = 1:length(tIdxArr)
    isInTRange = ((0:1439) == tIdxArr(iT));
    nWD = nanmean(WDDist(isInTRange, 1:end), 1);

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

    % polarimetry along height
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

    isNoRain = abs(bscAll) < 1e-8;
    vdr(~ isNoRain) = (bscAll(~ isNoRain) .* (depol_mean(~ isNoRain)) ./ (1 + depol_mean(~ isNoRain)) + mBsc(~ isNoRain) .* mdr / (1 + mdr)) ./ (bscAll(~ isNoRain) ./ (1 + depol_mean(~ isNoRain)) + mBsc(~ isNoRain) ./ (1 + mdr));
    vdr(isNoRain) = mdr;

    rainAmountFullArr(iT) = 4 * pi / 2 * sum(rWD(1:end).^3 .* nWD) * 1e-12;
    isInRRange = (rWD >= 0.5/2*1e3) & (rWD <= 1/2*1e3);
    rainAmount1mmArr(iT) = 4 * pi / 2 * sum(rWD(isInRRange).^3 .* nanmean(WDDist(isInTRange, isInRRange), 1)) * 1e-12;
    vdrArr(iT) = max(vdr((height >= 500) & (height <= 700)));
    secOrMmt1mmArr(iT) = sum(rWD(isInRRange).^2 .* nanmean(WDDist(isInTRange, isInRRange), 1));
    secOrMmtFullArr(iT) = sum(rWD(1:end).^2 .* nanmean(WDDist(isInTRange, 1:end), 1));
end

%% Visualization

% rain amount vs vdr
figure('Position', [0, 0, 800, 600], 'Units', 'Pixels', 'Color', 'w');

subplot(221);
scatter(vdrArr, rainAmountFullArr, 15, 'Marker', 'o', 'MarkerFaceColor', 'green'); hold on;

xlabel('Max \delta');
ylabel('Rain Amount (full)');

xlim([0, max(vdrArr)]);
ylim([0, max(rainAmountFullArr)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

% rain amount (d<1mm) vs vdr
subplot(222);

scatter(vdrArr, rainAmount1mmArr, 15, 'Marker', 'o', 'MarkerFaceColor', 'green'); hold on;

xlabel('Max \delta');
ylabel('Rain Amount (d<1mm)');

xlim([0, max(vdrArr)]);
ylim([0, max(rainAmount1mmArr)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

% second moment vs vdr
subplot(223);

scatter(vdrArr, secOrMmtFullArr, 15, 'Marker', 'o', 'MarkerFaceColor', 'green'); hold on;

xlabel('Max \delta');
ylabel('Second Moment (full)');

xlim([0, max(vdrArr)]);
ylim([0, max(secOrMmtFullArr)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

% second moment vs vdr
subplot(224);

scatter(vdrArr, secOrMmt1mmArr, 15, 'Marker', 'o', 'MarkerFaceColor', 'green'); hold on;

xlabel('Max \delta');
ylabel('Second Moment (d<1mm)');

xlim([0, max(vdrArr)]);
ylim([0, max(secOrMmt1mmArr)]);

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');

%% Save Data
a = cell(6, (length(vdrArr) + 1));
a{1, 1} = 'Time index';
a{2, 1} = 'vdr (max)';
a{3, 1} = 'Rain Amount (full)';
a{4, 1} = 'Rain Amount (d<1mm)';
a{5, 1} = 'Second Moment (full)';
a{6, 1} = 'Second Moment (d<1mm)';

for iT = 1:length(vdrArr)
    a{1, iT + 1} = tIdxArr(iT);
    a{2, iT + 1} = vdrArr(iT);
    a{3, iT + 1} = rainAmountFullArr(iT)';
    a{4, iT + 1} = rainAmount1mmArr(iT);
    a{5, iT + 1} = secOrMmtFullArr(iT);
    a{6, iT + 1} = secOrMmt1mmArr(iT);
end

xlswrite('simulation-outputs.xls', transpose(a));