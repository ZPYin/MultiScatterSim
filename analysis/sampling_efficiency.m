clc;

%% Parameter Definition
matFile = 'C:\Users\ZPYin\Documents\Coding\Matlab\MultiScatterSim\data\water_cloud_scattering_phase_function_example.mat';

%% Load Data
S = load(matFile);

%% Generate Probability Distribution
theta = S.theta;
r = S.r;
S1 = S.S1_Real + 1i * S.S1_Imag;
S2 = S.S2_Real + 1i * S.S2_Imag;
PhaseFunc = 1 / 2 * (abs(S1).^2 + abs(S2).^2);
PhaseFunc = transpose(PhaseFunc) * transpose(smmgDist(r, 'kappa', 7, 'sigma_e', 3e-2, 'R_eff', 10)) * (r(2) - r(1));
PhaseFunc = PhaseFunc / trapz(theta, PhaseFunc);

% test probability function
% PhaseFunc1 = (theta - 90).^2 * 3;
% PhaseFunc1 = PhaseFunc1 / trapz(theta, PhaseFunc1);

%% Rejection Sampling
tStart1 = now;
pfSamples = rejectSampling(theta, PhaseFunc, 1e7);
tEnd1 = now;

tStart2 = now;
pf1Samples = inverseSampling(theta, PhaseFunc, 1e7);
tEnd2 = now;

%% Display
figure('Position', [0, 0, 700, 600], 'Units', 'Pixels', 'Color', 'w');

subplot(211);
hold on;
p1 = histogram(pfSamples, 1000, 'Normalization', 'pdf');
p2 = plot(theta, PhaseFunc);
p1.EdgeColor = 'None';

xlabel('Scattering Angle (\circ)');
title(sprintf('Rejection Sampling (time usage: %f s)', (tEnd1 - tStart1) * 24 * 3600));

xlim([0, 180]);

set(gca, 'YScale', 'log');

subplot(212);

hold on;
p1 = histogram(pf1Samples, 1000, 'Normalization', 'pdf');
p2 = plot(theta, PhaseFunc);
p1.EdgeColor = 'None';

xlabel('Scattering Angle (\circ)');
title(sprintf('Inverse Sampling (time usage: %f s)', (tEnd2 - tStart2) * 24 * 3600));

xlim([0, 180]);

set(gca, 'YScale', 'log');