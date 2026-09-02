function optics = prepareMSMCOptics(lidarCfg, cloudCfg, varargin)
%PREPAREMSMCOPTICS Prepare optical data shared by CPU and GPU MSMC backends.
%
% optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...)
%
% This function separates the expensive Mie database reduction from photon
% transport.  The returned struct can be reused by MSMCSim and MSMCSimGPU
% through the 'preparedOptics' option.  Reusing it is especially important
% for repeated polarization states or Monte Carlo repeats of the same cloud.

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'lidarCfg', @isstruct);
addRequired(p, 'cloudCfg', @isstruct);
addParameter(p, 'scatMatrixFile', '', @ischar);
addParameter(p, 'cloud_distribution', 'modified-gamma', @ischar);
parse(p, lidarCfg, cloudCfg, varargin{:});

projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
scatMatrixFile = fullfile(projectDir, 'data', ...
    'water_droplet_scattering_532nm_finer_validation_r_step_0p01um.mat');
if ~isempty(p.Results.scatMatrixFile)
    scatMatrixFile = p.Results.scatMatrixFile;
end

extraArgs = reshape([fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'], 1, []);
cloudArgs = { ...
    'cloud_distribution', p.Results.cloud_distribution, ...
    'kappa', cloudCfg.kappa, ...
    'R_eff', cloudCfg.clReff, ...
    'b', 1/9};
cloudArgs = [cloudArgs, extraArgs];

if isfield(cloudCfg, 'N_d') && ~isempty(cloudCfg.N_d)
    cloudArgs = [cloudArgs, {'N_d', cloudCfg.N_d}];
else
    cloudArgs = [cloudArgs, {'sigma_e', cloudCfg.clExt}];
end

cloudScatMatrix = cloudScatteringMatrix( ...
    lidarCfg.lambda, scatMatrixFile, cloudArgs{:});

thetaGrid = cloudScatMatrix.theta(:).';
thetaWeights = [0, diff(thetaGrid)] .* sin(thetaGrid * pi / 180);

optics = struct();
optics.formatVersion = 1;
optics.lambda = lidarCfg.lambda;
optics.nLayers = numel(cloudScatMatrix.muExt);
optics.scatMatrixFile = scatMatrixFile;
optics.cloudScatMatrix = cloudScatMatrix;
optics.muExt = cloudScatMatrix.muExt(:);
optics.muSca = cloudScatMatrix.muSca(:);
optics.albedo = optics.muSca ./ optics.muExt;
optics.S11 = cloudScatMatrix.S11;
optics.S12 = cloudScatMatrix.S12;
optics.S33 = cloudScatMatrix.S33;
optics.S34 = cloudScatMatrix.S34;
optics.theta = thetaGrid;
optics.cdfS11 = cumsum(optics.S11 .* thetaWeights, 2);
optics.cdfS12 = cumsum(optics.S12 .* thetaWeights, 2);
end
