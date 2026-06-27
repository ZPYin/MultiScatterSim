clc;

%% Parameter Definition
r = 0.01:0.5:20;
theta = 175:0.1:180;
kappa = 7;
sigma_e = 4e-2;
R_eff = [4, 8, 13, 18];
inWL = 532;

%% Mie Calculation
nTheta = length(theta);
S1 = NaN(length(r), nTheta);
S2 = NaN(length(r), nTheta);
[m1, m2] = refIdxWater(inWL);
for iR = 1:length(r)
    x = 2 * pi * r(iR) / inWL * 1e3;

    for iTheta = 1:length(theta)
        [S1(iR, iTheta), S2(iR, iTheta)] = MieCalc(m1 + 1i * m2, x, theta(iTheta));
    end
end

%% Mueller Matrix Calculation
P11 = 1 / 2 * (abs(S1).^2 + abs(S2).^2);
factor = trapz(theta / 180 * pi, P11, 2);   % 将概率归一化
P12 = 1 / 2 * (abs(S2).^2 - abs(S1).^2);
P33 = 1 / 2 * (conj(S2) .* S1 + S2 .* conj(S1));
P11 = P11 ./ repmat(factor, 1, nTheta);
P12 = P12 ./ repmat(factor, 1, nTheta);
P33 = P33 ./ repmat(factor, 1, nTheta);

L1 = (P11 + P12);
L2 = (P11 - P12);
L3 = P33;

u = cos(repmat(theta, length(r), 1) / 180 * pi);
P_p = 1 / 8 * (3 * L1 + 3 * L2 .* u.^2 + 2 * L3 .* u);
P_c = 1 / 8 * (L1 + L2 .* u.^2 - 2 * L3 .* u);

P_p_int = NaN(length(R_eff), nTheta);
P_c_int = NaN(length(R_eff), nTheta);
for iReff = 1:length(R_eff)
    dis = smmgDist(r, 'kappa', kappa, 'sigma_e', sigma_e, 'R_eff', R_eff(iReff));
    P_p_int(iReff, :) = trapz(r, repmat(dis, nTheta, 1) .* P_p', 2) / ...
                        trapz(r, dis);
    P_c_int(iReff, :) = trapz(r, repmat(dis, nTheta, 1) .* P_c', 2) / ...
                        trapz(r, dis);
end

%% Display
figure('Position', [0, 0, 400, 1200], 'Units', 'Pixels', 'Color', 'w');

subPos = subfigPos([0.15, 0.1, 0.81, 0.87], 3, 1, 0, 0.04);

subplot('Position', subPos(1, :), 'Units', 'Normalized');

hold on;
p = [];
for iReff = 1:length(R_eff)
    p1 = plot(theta, P_p_int(iReff, :), 'LineWidth', 2, ...
              'DisplayName', sprintf('R_{eff}=%d\\mum', R_eff(iReff)));
    p = cat(1, p, p1);
end

xlabel('');
ylabel('Phase Func. P_{||}');

xlim([0, 180]);
ylim([1e-4, 1e4]);
grid on;

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', ...
         'YScale', 'log', 'XTickLabel', '');
legend(p, 'Location', 'NorthEast');

subplot('Position', subPos(2, :), 'Units', 'Normalized');
hold on;
p = [];
for iReff = 1:length(R_eff)
    p1 = plot(theta, P_c_int(iReff, :), ...
                'LineWidth', 2, ...
                'DisplayName', sprintf('R_{eff}=%d\\mum', R_eff(iReff)));
    p = cat(1, p, p1);
end

xlabel('');
ylabel('Phase Func. P_\perp');

xlim([0, 180]);
ylim([1e-5, 1e-2])
grid on;

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', ...
'YScale', 'log', 'XTickLabel', '');

subplot('Position', subPos(3, :), 'Units', 'Normalized');

hold on;
p = [];
for iReff = 1:length(R_eff)
    p1 = plot(theta, P_c_int(iReff, :) ./ P_p_int(iReff, :), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('R_{eff}=%d\\mum', R_eff(iReff)));
    p = cat(1, p, p1);
end

xlabel('散射角 (\circ)');
ylabel('退偏比');

xlim([177, 180]);
ylim([0, 1])
grid on;

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on', 'FontSize', 12);