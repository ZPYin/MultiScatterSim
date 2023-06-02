clc;
global MS_ENVS;

%% Parameter Definition
% lidar
lambda = 1064.0;   % incident wavelength. [nm]
rBeam = 5.0;   % radius of laser beam at z=0. [mm]
divBeam = 0.05;   % divergence of laser beam. [mrad]
FOV = 1;   % field of view. [mrad]
rTel = 150;   % radius of receiving telescope. [mm]
binWidth = 20.0;   % bin width. [ns]
nBins = 100;   % number of range bins.
distBeamTel = 50;   % distance between laser axis and optics axis of receiving telescope. [mm]
SBeam = [1, 1, 0, 0];   % Stokes vector of the incident laser beam.

% simulation
nPhotons = 1e4;

% medium
clBase = 1000;   % cloud base height. [m]
kappa = 9;   % shape factor of modified gamma distribution.
clExt = 2e-2;   % cloud extinction. [m-1]
clReff = 9;   % cloud droplet effective radius. [um]
nClLayers = 30;   % number of cloud layers.
clDh = 10;   % geometrical thickness of each cloud layer. [m]
clScatMatrixFile = fullfile(MS_ENVS.RootPath, 'data', 'water_cloud_scattering_phase_function_example.mat');

%% Load Cloud Scattering Properties
mieScaInfo = load(clScatMatrixFile);

% Calculate Scattering Properties of Water Cloud
isWLinDB = any(abs(mieScaInfo.inWL - lambda) < 1);
if ~ isWLinDB
    error('Scattering properties at given wavelength (%f) is not in database.', lambda);
end
wlIdx = abs(mieScaInfo.inWL - lambda) < 1;
muExt = pi * trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff) .* (mieScaInfo.r * 1e-6).^2 .* squeeze(mieScaInfo.Qext(wlIdx, :)));
muSca = pi * trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff) .* (mieScaInfo.r * 1e-6).^2 .* squeeze(mieScaInfo.Qsca(wlIdx, :)));
albedo = muSca / muExt;
S1 = squeeze(mieScaInfo.S1_Real(wlIdx, :, :)) + 1i * squeeze(mieScaInfo.S1_Imag(wlIdx, :, :));
S2 = squeeze(mieScaInfo.S2_Real(wlIdx, :, :)) + 1i * squeeze(mieScaInfo.S2_Imag(wlIdx, :, :));
S11 = 0.5 * (abs(S1).^2 + abs(S2).^2);
S12 = 0.5 * (abs(S2).^2 - abs(S1).^2);
S33 = 0.5 * (conj(S2) .* S1 + S2 .* conj(S1));
S34 = -0.5i * (S1 .* conj(S2) - S2 .* conj(S1));

factor = trapz(mieScaInfo.theta / 180 * pi, S11, 2);
S11 = S11 ./ repmat(factor, 1, length(mieScaInfo.theta));
S12 = S12 ./ repmat(factor, 1, length(mieScaInfo.theta));
S33 = S33 ./ repmat(factor, 1, length(mieScaInfo.theta));
S34 = S34 ./ repmat(factor, 1, length(mieScaInfo.theta));
S11 = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff)) .* S11, 1) / ...
      trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff));
S12 = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff)) .* S12, 1) / ...
      trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff));
S33 = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff)) .* S33, 1) / ...
      trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff));
S34 = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff)) .* S34, 1) / ...
      trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', kappa, 'sigma_e', clExt, 'R_eff', clReff));

% Cloud Boundary
% Usually bound corresponding to an optical depth of 2, according to DM Winker and LR Poole, "Monte-Carlo calculations of cloud returns for ground-based and space-based lidars," Applied Physics B 60 (4), 341-344 (1995).   
% Note: we set 4.
clBoundX = 4.0 / muExt;
clBoundY = 4.0 / muExt;
clTop = clBase + clDh * nClLayers;

%% MC
pPos = [0, 0, 0];   % position of photon
pDir = [0, 0, 0];   % transmit direction of photon
len = zeros(1, nPhotons);   % transmit length
SVReturn = zeros(4, nBins);   % Stokes Vector of return signal
nScas = zeros(1, nPhotons);   % number of scattering events for each photon
for iP = 1:nPhotons

    % launch to cloud base
    isAlive = true;
    rCl = rBeam * 1e-3 + divBeam * 1e-3 / 2 * clBase * 1e-3;
    pPos(1) = rand(1, 1) * rCl;
    pPos(2) = rand(1, 1) * rCl;
    pPos(3) = clBase;
    while (norm(pPos(1:2), 2) > rCl)
        % generate photon according to initial distribution of photons
        pPos(1) = rand(1, 1) * rCl;
        pPos(2) = rand(1, 1) * rCl;
    end
    pDir = pPos / norm(pPos - [0, 0, 0], 2);
    SVIn = rotSVPhi(SBeam, atan(pDir(2) / pDir(1)));
    iLayer = 1;

    while (isAlive)

        % move to the next point
        rnd = rand(1, 1) + 1e-10;
        ds = -log(rnd) / muExt;
        pPos = pPos + pDir * ds;
        len(iP) = len(iP) + ds;
        iLayer = ceil((pPos(3) - clBase) / clDh);

        % still in the cloud?
        if (isInMedium(pPos, [clBase, clTop, clBoundX, clBoundY]))

            isAbsorbed = (rand(1, 1) > albedo);
            if isAbsorbed
                isAlive = false;
                break;
            end

            diaH = FOV * 1e-3 * pPos(3) + rTel * 1e-3 * 2;   % diameter of FOV at height H
            isInFOV = sqrt((pPos(1) - distBeamTel * 1e-3)^2 + pPos(2)^2) < (diaH / 2);
            if (isInFOV)
                % Force scattering into the FOV (semi-analytic method)

                pInDir = pDir;
                pOutDir = [distBeamTel * 1e-3, 0, 0] - pPos;
                pOutDir = pOutDir / norm(pOutDir, 2);

                % scattering rotation angle
                [scaAng, rotAng1, rotAng2] = scaRotAng(pInDir, pOutDir);

                % solid angle
                solAng = pi * (rTel * 1e-3)^2 / norm(pPos - [distBeamTel * 1e-3, 0, 0], 2).^2;

                % probability without scattering again until entering lidar FOV
                [~, iAng] = min(abs(mieScaInfo.theta - scaAng / pi * 180));
                if mieScaInfo.theta(iAng) < 90
                    probEnter = 0;
                else
                    probEnter = albedo * solAng * exp(- ((pPos(3) - clBase) * muExt) / abs(pOutDir(3))) * S11(iAng);
                end

                % Stokes vector after scattering
                SVTemp1 = rotSVPhi(SVIn, rotAng1);
                MuMatrix = [S11(iAng), S12(iAng), 0, 0; ...
                            S12(iAng), S11(iAng), 0, 0; ...
                            0, 0, S33(iAng), S34(iAng); ...
                            0, 0, -S34(iAng), S33(iAng)];
                SVTemp2 = MuMatrix * transpose(SVTemp1) / S11(iAng);
                SVSca = rotSVPhi(SVTemp2, -rotAng2);

                tReturn = (len(iP) + (pPos(3) - clBase) / abs(pOutDir(3))) / 3e8;
                iBin = ceil(tReturn / (binWidth * 1e-9));
                if iBin <= nBins
                    SVReturn(:, iBin) = SVReturn(:, iBin) + transpose(rotSVPhi(SVSca, atan(pOutDir(2) / pOutDir(1))) * probEnter);
                end

                SVIn = SVIn * (1 - probEnter);

            end

            % Scattering
            phi = rand(1, 1) * 2 * pi;
            phaseFunc = S11 * SVIn(1) + S12 * (SVIn(2) * cos(2 * phi) + SVIn(3) * sin(2 * phi));
            theta = inverseSampling(mieScaInfo.theta, phaseFunc, 1) / 180 * pi;

            % Stokes vector of scattered nPhotons
            [~, iAng] = min(abs(mieScaInfo.theta - theta / pi * 180));
            pInDir = pDir;
            pOutDir = updateDir(pInDir, phi, theta);
            SVTemp1 = rotSVPhi(SVIn, phi);
            MuMatrix = [S11(iAng), S12(iAng), 0, 0; ...
                        S12(iAng), S11(iAng), 0, 0; ...
                        0, 0, S33(iAng), S34(iAng); ...
                        0, 0, -S34(iAng), S33(iAng)];
            SVTemp2 = MuMatrix * transpose(SVTemp1) / S11(iAng);
            temp = sqrt((1 - cos(theta)^2) * (1 - pOutDir(3)^2));
            if abs(temp) < 1e-8
                gammaAng = pi / 2;
            else
                gammaAngCos = (double((phi > pi) && (phi < 2*pi)) - 0.5) * 2 * (pOutDir(3) * cos(theta) - pInDir(3)) / temp;
                gammaAngCos = double(gammaAngCos <= -1) * -1 + ...
                              double(gammaAngCos >= 1) * 1 + ...
                              double((gammaAngCos > -1) && (gammaAngCos < 1)) * gammaAngCos;
                gammaAng = acos(gammaAngCos);
            end

            SVIn = rotSVPhi(SVTemp2, -gammaAng);   % rotate meridian plane
            nScas(iP) = nScas(iP) + 1;
            pDir = pOutDir;
        else
            isAlive = false;
            fprintf('Photon %d is dead\n', iP);
        end

    end
end

%% Display
figure('Color', 'w');

subplot(311);

histogram(nScas);

subplot(312);

hold on;
plot((1:nBins) * binWidth * 0.15, SVReturn(1, :));
plot((1:nBins) * binWidth * 0.15, SVReturn(2, :));
plot((1:nBins) * binWidth * 0.15, SVReturn(3, :));
plot((1:nBins) * binWidth * 0.15, SVReturn(4, :));
xlim([0, nBins * binWidth * 0.15]);
set(gca, 'YScale', 'log');

subplot(313);
depRatio = (SVReturn(1, :) - SVReturn(2, :)) ./ (SVReturn(1, :) + SVReturn(2, :));
%depRatio = (SVReturn(1, :) + SVReturn(4, :)) ./ (SVReturn(1, :) - SVReturn(4, :));
plot((1:nBins) * binWidth * 0.15, depRatio);

ylim([0, 1]);
xlim([0, nBins * binWidth * 0.15]);