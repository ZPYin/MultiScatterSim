clc;
close all;

%% Parameter Definition
inWL = 905;   % incident wavelength. (nm)
r = logspace(log10(0.05), log10(1000), 1000);   % radius of water droplet. (um)

%% Mie Calculation
Qbsc = NaN(1, length(r));
Qsca = NaN(1, length(r));
Qext = NaN(1, length(r));

[m1, m2] = refIdxWater(inWL);
for iR = 1:length(r)
    x = 2 * pi * r(iR) / inWL * 1e3;

    [~, ~, Qext(iR), Qsca(iR), Qbsc(iR)] = MieCalc(m1 + 1i * m2, x, 180);
end

%% Display
figure('Position', [0, 0, 400, 400], 'Units', 'Pixels', 'Color', 'w');

hold on;

p1 = plot(r * 2, Qext, 'color', [99, 128, 120] / 255, 'DisplayName', 'Extinction Eff.');
p2 = plot(r * 2, Qbsc, 'color', [255, 217, 141] / 255, 'DisplayName', 'Backscatter Eff.');

hold off;

grid on;

xlim([0.1, 2000]);
ylim([1e-4, 1e2]);

xlabel('Diameter (\mum)');
ylabel('Efficienty');
title('Extinction/Backscatter efficienty for water droplets');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XScale', 'log', 'YScale', 'log', 'Box', 'on');

legend([p1, p2], 'Location', 'SouthEast');

%% Lidar Ratio Calculation
% normalized gamma distribution
muArr = [2, 10];   % shape factor
MVD = logspace(log10(1), log10(2000), 1000);
lr = NaN(length(muArr), length(MVD));
for iMu = 1:length(muArr)
    for iMVD = 1:length(MVD)
        nD = @(x) ((x / MVD(iMVD)).^muArr(iMu) .* exp(-(3.67 + muArr(iMu)) .* x / MVD(iMVD)));
        ext = trapz(r, nD(2 * r) .* Qext .* pi .* (r * 1e-6).^2);
        bsc = trapz(r, nD(2 * r) .* Qbsc .* pi .* (r * 1e-6).^2) / (4 * pi);

        lr(iMu, iMVD) = ext / bsc;
    end
end

%% Display
figure('Position', [0, 0, 400, 400], 'Units', 'Pixels', 'Color', 'w');

hold on;

p1 = plot(MVD, lr(1, :), 'color', [99, 128, 120] / 255, 'DisplayName', sprintf('\\mu=%d', muArr(1)));
p2 = plot(MVD, lr(2, :), 'color', [255, 217, 141] / 255, 'DisplayName', sprintf('\\mu=%d', muArr(2)));

hold off;

grid on;

xlim([1, 2000]);
ylim([1, 200]);

xlabel('Median volume diameter, D_0 (\mum)');
ylabel('lidar ratio (sr)');
title('lidar ratio for water droplets with different sizes');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XScale', 'log', 'YScale', 'log', 'Box', 'on');

legend([p1, p2], 'Location', 'SouthWest');