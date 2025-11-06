function [S1, S2] = MieSim4RainDrops(r, theta, inWL, varargin)
% MIESIM4RAINDROPS Mie Scattering simulation for certain shperes.
%
% USAGE:
%    [S1, S2] = MieSim4RainDrops(r, theta, inWL)
%
% INPUTS:
%    r: array
%        radius of rain drops. (um)
%    theta: array
%        scattering angles. (degree)
%    inWL: array
%        incident wavelengths. (nm)
%
% KEYWORDS:
%    saveFile
%
% OUTPUTS:
%    S1: 3-D matrix (wavelength x radius x theta)
%        S1 component from Mie Scattering.
%    S2: 3-D matrix (wavelength x radius x theta)
%        S2 component from Mie Scattering.
%    Qext: 3-D matrix (wavelength x radius x theta)
%        extinction efficiency.
%    Qsca: 3-D matrix (wavelength x radius x theta)
%        scattering efficiency.
%
% HISTORY:
%    2023-10-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'r', @isnumeric);
addRequired(p, 'theta', @isnumeric);
addRequired(p, 'inWL', @isnumeric);
addParameter(p, 'saveFile', '', @ischar);
addParameter(p, 'debug', false, @islogical);

parse(p, r, theta, inWL, varargin{:});

%% Parameter Definition
% r = [0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.250,3.750,4.250,4.750,5.500] / 2 * 1e3;   % radius of water droplet. [um]
% r = (4000:0.1:4001) / 2;   % radius of water droplet. [um]
% theta = 180 - (0:0.001:0.1);
% inWL = 355;

%% Mie Scattering
S1 = NaN(length(inWL), length(r), length(theta));
S2 = NaN(length(inWL), length(r), length(theta));
Qext = NaN(length(inWL), length(r), length(theta));
Qsca = NaN(length(inWL), length(r), length(theta));

for iWL = 1:length(inWL)
    [m1, m2] = refIdxWater(inWL(iWL));

    for iR = 1:length(r)
        x = 2 * pi * r(iR) / inWL(iWL) * 1e3;

        for iTheta = 1:length(theta)
            if (p.Results.debug)
                fprintf('Current sim: wavelength=%f nm; r=%f um; theta=%f\deg\n', inWL(iWL), r(iR), theta(iTheta));
            end

            [S1(iWL, iR, iTheta), S2(iWL, iR, iTheta), Qext(iWL, iR, iTheta), Qsca(iWL, iR, iTheta), ] = MieCalc(m1 + 1i * m2, x, theta(iTheta));
        end
    end
end

% P11 = 1 / 2 * (abs(S1).^2 + abs(S2).^2);
% factor = trapz(theta / 180 * pi, P11, 2);   % 将概率归一化
% P12 = 1 / 2 * (abs(S2).^2 - abs(S1).^2);
% P33 = 1 / 2 * (conj(S2) .* S1 + S2 .* conj(S1));
% P11 = P11 ./ repmat(factor, 1, length(theta));
% P12 = P12 ./ repmat(factor, 1, length(theta));
% P33 = P33 ./ repmat(factor, 1, length(theta));

% L1 = (P11 + P12);
% L2 = (P11 - P12);
% L3 = P33;

% u = cos(repmat(theta, length(r), 1) / 180 * pi);
% P_p = 1 / 8 * (3 * L1 + 3 * L2 .* u.^2 + 2 * L3 .* u);
% P_c = 1 / 8 * (L1 + L2 .* u.^2 - 2 * L3 .* u);

% save(fullfile(MS_ENVS.RootPath, 'data', saveFile), 'theta', 'r', 'inWL', 'P_p', 'P_c', 'P11', 'P12', 'P33');
if ~ isempty(p.Results.saveFile)
    save(p.Results.saveFile, 'theta', 'r', 'inWL', 'S1', 'S2', 'Qext', 'Qsca');
end

end