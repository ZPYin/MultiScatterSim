function [SVReturn] = MSMCSim(lidarCfg, cloudCfg, varargin)
% MSMCSIM Multiple-Scattering with Monte-Carlo method.
%
% USAGE:
%    [SVReturn] = MSMCSim(lidarCfg, cloudCfg)
%
% INPUTS:
%    lidarCfg: struct
%        lambda: incident wavelength. [nm]
%        rBeam: radius of laser beam at z=0. [mm]
%        divBeam: divergence of laser beam. [mrad]
%        FOV: field of view. [mrad]
%        rTel: radius of receiving telescope. [mm]
%        binWidth: bin width. [ns]
%        nBins: number of range bins.
%        distBeamTel: distance between laser axis and optics axis of receiving telescope. [mm]
%        SBeam: Stokes vector of the incident laser beam.
%    cloudCfg: struct
%        clBase: cloud base height. [m]
%        clDh: geometrical thickness of each cloud layer. [m]
%        nClLayers: number of cloud layers.
%        kappa: shape factor of modified gamma distribution.
%        clExt: cloud extinction. [m-1]
%        clReff: cloud droplet effective radius. [um]
%
% KEYWORDS:
%    nPhoton: numeric
%        number of photons to be simulated. (default: 1e4)
%    visible: char
%        'on' | 'off'
%    debug: char
%
% OUTPUTS:
%    SVReturn:
%        4 * nBins [I, Q, U, V] at each range gate.
%
% HISTORY:
%    2023-06-05: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'lidarCfg', @isstruct);
addRequired(p, 'cloudCfg', @isstruct);
addParameter(p, 'nPhoton', 1e4, @isnumeric);
addParameter(p, 'visible', 'on', @ischar);
addParameter(p, 'debug', false, @islogical);

parse(p, lidarCfg, cloudCfg, varargin{:});

%% Parameter Definition
nPhotons = p.Results.nPhoton;

%% Cloud Scattering Properties
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
clScatMatrixFile = fullfile(projectDir, 'data', 'water_cloud_scattering_phase_function_example.mat');
cloudScatMatrix = cloudScatteringMatrix(lidarCfg.lambda, clScatMatrixFile, 'cloud_distribution', 'modified-gamma', 'kappa', cloudCfg.kappa, 'sigma_e', cloudCfg.clExt, 'R_eff', cloudCfg.clReff);
muExt = cloudScatMatrix.muExt;
muSca = cloudScatMatrix.muSca;
albedo = muSca ./ muExt;
S11 = cloudScatMatrix.S11;
S12 = cloudScatMatrix.S12;
S33 = cloudScatMatrix.S33;
S34 = cloudScatMatrix.S34;

% Cloud Boundary
% Usually bound corresponding to an optical depth of 2, according to DM Winker and LR Poole, "Monte-Carlo calculations of cloud returns for ground-based and space-based lidars," Applied Physics B 60 (4), 341-344 (1995).   
% Note: we set 4.
clBoundX = 4.0 ./ muExt;
clBoundY = 4.0 ./ muExt;

%% Monte-Carlo Simulation
rCl = lidarCfg.rBeam * 1e-3 + lidarCfg.divBeam * 1e-3 / 2 * cloudCfg.clBase * 1e-3;   % radius of laser beam at cloud base
SVReturn = zeros(4, lidarCfg.nBins);   % Stokes Vector of return signal
len = zeros(1, nPhotons);   % transmit length
nScas = zeros(1, nPhotons);   % number of scattering events for each photon
for iPhoton = 1:nPhotons
    pPos = [0; 0; 0];
    temp = zeros(1, 5);

    % launch to cloud base
    isAlive = true;
    pPos(1) = rand(1, 1) * rCl;
    pPos(2) = rand(1, 1) * rCl;
    pPos(3) = cloudCfg.clBase;
    while (norm(pPos(1:2), 2) > rCl)
        % generate photon according to initial distribution of photons
        pPos(1) = rand(1, 1) * rCl;
        pPos(2) = rand(1, 1) * rCl;
    end
    pDir = pPos / norm(pPos - [0; 0; 0], 2);
    SVIn = rotSVPhi(lidarCfg.SBeam, atan(pDir(2) / pDir(1)));

    while (isAlive)

        % move to the next point
        [pPos, isInMedium, ds, scaClLayer] = freeRunning(pPos, pDir, muExt, cloudCfg.clDh, cloudCfg.clBase, clBoundX, clBoundY);
        len(iPhoton) = len(iPhoton) + ds;

        % still in the cloud?
        if (isInMedium)

            if (rand(1, 1) > albedo(scaClLayer))
                % if the photon was absorbed
                isAlive = false;
                break;
            end

            % whether the photon was still in the FOV
            if (sqrt((pPos(1) - lidarCfg.distBeamTel * 1e-3)^2 + pPos(2)^2) < ((lidarCfg.FOV * 1e-3 * pPos(3) + lidarCfg.rTel * 1e-3 * 2) / 2))
                % Force scattering into the FOV (semi-analytic method)
                pScaDir = [lidarCfg.distBeamTel * 1e-3; 0; 0] - pPos;
                pScaDir = pScaDir / norm(pScaDir, 2);

                % scattering rotation angle
                [scaAng, rotAng1, rotAng2] = scaRotAng(pDir, pScaDir);

                % probability without scattering again until entering lidar FOV
                [~, iAng] = min(abs(cloudScatMatrix.theta - scaAng / pi * 180));
                if cloudScatMatrix.theta(iAng) < 90
                    probEnter = 0;
                else
                    probEnter = albedo(scaClLayer) * ...
                        pi * (lidarCfg.rTel * 1e-3)^2 / norm(pPos - [lidarCfg.distBeamTel * 1e-3; 0; 0], 2).^2 * ...
                        exp(- ((pPos(3) - cloudCfg.clBase) * muExt(scaClLayer)) / abs(pScaDir(3))) * ...
                        S11(scaClLayer, iAng);
                end

                % Stokes vector after scattering
                SVTemp1 = rotSVPhi(SVIn, rotAng1);
                SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                                       S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                                       0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                                       0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                                    SVTemp1 / S11(scaClLayer, iAng);

                iBin = ceil((len(iPhoton) + (pPos(3) - cloudCfg.clBase) / abs(pScaDir(3))) / 3e8 / (lidarCfg.binWidth * 1e-9));
                if iBin <= lidarCfg.nBins
                    SVReturn(:, iBin) = SVReturn(:, iBin) + rotSVPhi(rotSVPhi(SVTemp2, -rotAng2), atan(pScaDir(2) / pScaDir(1)) * probEnter * 2^(nScas(iPhoton) - 1));
                end

                SVIn = SVIn * (1 - probEnter);

            end

            % Scattering
            rotAng1 = rand(1, 1) * 2 * pi;
            scaAng = inverseSampling(cloudScatMatrix.theta, ...
                S11(scaClLayer, :) * SVIn(1) + S12(scaClLayer, :) * (SVIn(2) * cos(2 * rotAng1) + SVIn(3) * sin(2 * rotAng1)), ...
                1) / 180 * pi;

            % Stokes vector of scattered nPhotons
            [~, iAng] = min(abs(cloudScatMatrix.theta - scaAng / pi * 180));
            pScaDir = updateDir(pDir, rotAng1, scaAng);
            SVTemp1 = rotSVPhi(SVIn, rotAng1);
            SVTemp2 = [S11(scaClLayer, iAng), S12(scaClLayer, iAng), 0, 0; ...
                                   S12(scaClLayer, iAng), S11(scaClLayer, iAng), 0, 0; ...
                                   0, 0, S33(scaClLayer, iAng), S34(scaClLayer, iAng); ...
                                   0, 0, -S34(scaClLayer, iAng), S33(scaClLayer, iAng)] * ...
                            SVTemp1 / S11(scaClLayer, iAng);
            temp(1) = sqrt((1 - cos(scaAng)^2) * (1 - pScaDir(3)^2));
            if abs(temp(1)) < 1e-8
                temp(2) = pi / 2;
            else
                temp(3) = (double((rotAng1 > pi) && (rotAng1 < 2*pi)) - 0.5) * 2 * (pScaDir(3) * cos(scaAng) - pDir(3)) / temp(1);
                temp(3) = double(temp(3) <= -1) * -1 + ...
                                   double(temp(3) >= 1) * 1 + ...
                                   double((temp(3) > -1) && (temp(3) < 1)) * temp(3);
                temp(2) = acos(temp(3));
            end

            SVIn = rotSVPhi(SVTemp2, -temp(2));   % rotate meridian plane
            nScas(iPhoton) = nScas(iPhoton) + 1;
            pDir = pScaDir;
        else
            isAlive = false;

            if p.Results.debug
                fprintf('Finished %6.2f%%: Photon %d is dead\n', iPhoton/ nPhotons * 100, iPhoton);
            end
        end

    end
end

%% Display
if strcmp(p.Results.visible, 'on')
    figure('Color', 'w');

    subplot(311);

    histogram(nScas);

    subplot(312);

    hold on;
    plot((1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15, SVReturn(1, :));
    plot((1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15, SVReturn(2, :));
    plot((1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15, SVReturn(3, :));
    plot((1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15, SVReturn(4, :));
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
    set(gca, 'YScale', 'log');

    subplot(313);
    depRatio = (SVReturn(1, :) - SVReturn(2, :)) ./ (SVReturn(1, :) + SVReturn(2, :));
    %depRatio = (SVReturn(1, :) + SVReturn(4, :)) ./ (SVReturn(1, :) - SVReturn(4, :));
    plot((1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15, depRatio);

    ylim([0, 1]);
    xlim([0, lidarCfg.nBins * lidarCfg.binWidth * 0.15]);
end

end