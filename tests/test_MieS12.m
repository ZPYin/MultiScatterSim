clc;
global MS_ENVS;

%% Parameter Definition
wl = 532;
r = 5;   % radius. [um]
m = 1.377 + 1.4992e-09i;
S1CmpFile = 'Mie_S1Data_for_compare.txt';
S2CmpFile = 'Mie_S2Data_for_compare.txt';

%% Read Comparison Data
fid1 = fopen(fullfile(MS_ENVS.RootPath, 'data', S1CmpFile), 'r');
fgetl(fid1);
fgetl(fid1);
wlStr = fgetl(fid1);
wlCmp = str2double(strsplit(wlStr(11:end), '\t'));
isWL = abs(wlCmp - wl) < 1e-2;
if ~ any(isWL)
    error('scattering properties at given wavelength was not simulated');
end

data = textscan(fid1, repmat('%f', 1, length(wlCmp) * 2 - 1), 'delimiter', '\t');
angS1 = data{1};
S1Cmp = data{find(isWL) * 2} + data{find(isWL) * 2 + 1} * 1i;

fclose(fid1);

fid2 = fopen(fullfile(MS_ENVS.RootPath, 'data', S2CmpFile), 'r');
fgetl(fid2);
fgetl(fid2);
wlStr = fgetl(fid2);
wlCmp = str2double(strsplit(wlStr(11:end), '\t'));
isWL = abs(wlCmp - wl) < 1e-2;
if ~ any(isWL)
    error('scattering properties at given wavelength was not simulated');
end

data = textscan(fid2, repmat('%f', 1, length(wlCmp) * 2 - 1), 'delimiter', '\t');
angS2 = data{1};
S2Cmp = data{find(isWL) * 2} + data{find(isWL) * 2 + 1} * 1i;

fclose(fid2);

%% Mie Calculation
S1 = NaN(1, length(angS1));
S2 = NaN(1, length(angS1));
for iT = 1:length(angS1)
    [S1(iT), S2(iT)] = MieCalc(m, 2 * pi * r / (wl * 1e-3), angS1(iT));
end

%% Display
figure('color', 'w');

subplot(211);

hold on;
p1 = semilogy(angS1, abs(S1Cmp), 'r-', 'DisplayName', 'Comparison');
p2 = semilogy(angS1, abs(S1), 'b--', 'DisplayName', 'matlab');

xlabel('Scattering Angle');
ylabel('|S1|');

ylim([1e-1, 1e5]);
xlim([0, 180]);

legend([p1, p2]);

set(gca, 'YScale', 'log');

subplot(212);

hold on;
p1 = semilogy(angS1, abs(S2Cmp), 'r-', 'DisplayName', 'Comparison');
p2 = semilogy(angS1, abs(S2), 'b--', 'DisplayName', 'matlab');

xlabel('Scattering Angle');
ylabel('|S2|');
ylim([1e-1, 1e5]);
xlim([0, 180]);

legend([p1, p2]);

set(gca, 'YScale', 'log');