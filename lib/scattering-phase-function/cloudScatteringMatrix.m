function [cloudScatMatrix] = cloudScatteringMatrix(lambda, wdScatFile, varargin)
% CLOUDSCATTERINGMATRIX Calculate Cloud Scattering Matrix.
%
% USAGE:
%    [cloudScatMatrix] = cloudScatteringMatrix(lambda, wdScatFile)
%
% INPUTS:
%    lambda: numeric
%        incident wavelength. (nm)
%    wdScatFile: char
%        absolute path of water droplets scattering lookup table. (see ./lib/scattering-phase-function/water_droplet_scattering_gen.m)
%
% KEYWORDS:
%    cloud_distribution: char
%        'modifed-gamma' or 'log-normal' (not supported yet)
%    kappa: numeric
%        shape factor of modified-gamma distribution.
%    sigma_e: numeric
%        cloud extinction. (m-1)
%    R_eff: numeric
%        cloud effective radius. (um)
%
% OUTPUTS:
%    cloudScatMatrix: struct
%        muExt (m-1)
%        muSca (m-1)
%        S11
%        S12
%        S33
%        S34
%
% HISTORY:
%    2023-06-02: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'lambda', @isnumeric);
addRequired(p, 'wdScatFile', @ischar);
addParameter(p, 'cloud_distribution', 'modified-gamma', @ischar);
addParameter(p, 'kappa', 7, @isnumeric);
addParameter(p, 'sigma_e', 2e-2, @isnumeric);
addParameter(p, 'R_eff', 5, @isnumeric);

parse(p, lambda, wdScatFile, varargin{:});

cloudScatMatrix = struct();

mieScaInfo = load(wdScatFile);

switch lower(p.Results.cloud_distribution)
case 'modified-gamma'
    nClLayer = length(p.Results.sigma_e);
    nTheta = length(mieScaInfo.theta);

    isWLinDB = any(abs(mieScaInfo.inWL - lambda) < 1);
    if ~ isWLinDB
        error('Scattering properties at given wavelength (%f) is not in database.', lambda);
    end
    wlIdx = abs(mieScaInfo.inWL - lambda) < 1;

    muExt = NaN(1, nClLayer);
    muSca = NaN(1, nClLayer);
    S11_PD = NaN(nClLayer, nTheta);
    S12_PD = NaN(nClLayer, nTheta);
    S33_PD = NaN(nClLayer, nTheta);
    S34_PD = NaN(nClLayer, nTheta);

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

    for iLayer = 1:nClLayer
        muExt(iLayer) = pi * trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)) .* (mieScaInfo.r * 1e-6).^2 .* squeeze(mieScaInfo.Qext(wlIdx, :)));
        muSca(iLayer) = pi * trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)) .* (mieScaInfo.r * 1e-6).^2 .* squeeze(mieScaInfo.Qsca(wlIdx, :)));

        S11_PD(iLayer, :) = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer))) .* S11, 1) / trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)));
        S12_PD(iLayer, :) = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer))) .* S12, 1) / trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)));
        S33_PD(iLayer, :) = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer))) .* S33, 1) / trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)));
        S34_PD(iLayer, :) = trapz(mieScaInfo.r, transpose(smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer))) .* S34, 1) / trapz(mieScaInfo.r, smmgDist(mieScaInfo.r, 'kappa', p.Results.kappa(iLayer), 'sigma_e', p.Results.sigma_e(iLayer), 'R_eff', p.Results.R_eff(iLayer)));
    end

    cloudScatMatrix.muExt = muExt;
    cloudScatMatrix.muSca = muSca;
    cloudScatMatrix.S11 = S11_PD;
    cloudScatMatrix.S12 = S12_PD;
    cloudScatMatrix.S33 = S33_PD;
    cloudScatMatrix.S34 = S34_PD;
    cloudScatMatrix.theta = mieScaInfo.theta;

case 'log-normal'
otherwise
    error('Unknown cloud microphysical distribution: %s', p.Results.cloud_distribution);
end

end