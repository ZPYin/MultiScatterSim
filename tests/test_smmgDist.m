clc;

%% Parameter Definition
Re = [4, 8, 13, 18];   % [¦Ìm]
sigma_e = [7e-2, 3e-1, 6e-1, 1];
kappa = 7;

%% Display
r = 0:0.001:20;   % [¦Ìm]

figure('Position', [0, 0, 500, 350], 'Units', 'Pixels', 'Color', 'w');

p1 = plot(r, smmgDist(r, 'kappa', kappa, 'R_eff', Re(1), 'sigma_e', sigma_e(1)) * 1e-6, '-k', 'LineWidth', 2, 'DisplayName', sprintf('R_e=%d\\mum', Re(1)));
hold on;
p2 = plot(r, smmgDist(r, 'kappa', kappa, 'R_eff', Re(2), 'sigma_e', sigma_e(2)) * 1e-6, '-b', 'LineWidth', 2, 'DisplayName', sprintf('R_e=%d\\mum', Re(2)));
p3 = plot(r, smmgDist(r, 'kappa', kappa, 'R_eff', Re(3), 'sigma_e', sigma_e(3)) * 1e-6, '-r', 'LineWidth', 2, 'DisplayName', sprintf('R_e=%d\\mum', Re(3)));
p4 = plot(r, smmgDist(r, 'kappa', kappa, 'R_eff', Re(4), 'sigma_e', sigma_e(4)) * 1e-6, '-g', 'LineWidth', 2, 'DisplayName', sprintf('R_e=%d\\mum', Re(4)));

xlabel('Radius (\mum)');
ylabel('Gamma Dist. (\mum^{-1}cm^{-3})');

xlim([0, 20]);
ylim([0, 400]);
legend([p1, p2, p3, p4], 'Location', 'NorthEast');

set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'Box', 'on');