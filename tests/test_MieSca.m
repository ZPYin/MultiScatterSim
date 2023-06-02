clc;

%% Parameter Definition
theta = 0:1:180;
R = 4;   % radius of water droplet. [¦Ìm]
inWL = 532;   % [nm]

%% Mie Calculation
[m1, m2] = refIdxWater(inWL);
P11 = NaN(1, length(theta));
P12 = NaN(1, length(theta));
P33 = NaN(1, length(theta));
P34 = NaN(1, length(theta));

for iTheta = 1:length(theta)
    x = 2 * pi * R / inWL * 1e3;
    [S1, S2] = MieCalc(m1 + i * m2, x, theta(iTheta));

    P11(iTheta) = 1 / 2 * (abs(S2).^2 + abs(S1).^2);
    P12(iTheta) = 1 / 2 * (abs(S2).^2 - abs(S1).^2);
    P33(iTheta) = 1 / 2 * (conj(S2) * S1 + S2 * conj(S1));
    P34(iTheta) = -i / 2 * (conj(S2) * S1 - S2 * conj(S1));
end

%% Display
figure();

p1 = polar(theta / 180 * pi, P11);