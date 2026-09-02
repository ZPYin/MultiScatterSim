%% Setup
run(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'setupMSToolbox.m'));

%% Parameter Definition
% lidar
lidarCfg.lambda = 532.0;   % incident wavelength. [nm]
lidarCfg.rBeam = 5.0;   % radius of laser beam at z=0. [mm]
lidarCfg.divBeam = 0.05;   % divergence of laser beam. [mrad]
lidarCfg.FOV = 1;   % field of view. [mrad]
lidarCfg.rTel = 150;   % radius of receiving telescope. [mm]
lidarCfg.binWidth = 20.0;   % bin width. [ns]
lidarCfg.nBins = 100;   % number of range bins.
lidarCfg.distBeamTel = 50;   % distance between laser axis and optics axis of receiving telescope. [mm]
lidarCfg.SBeam = [1; 1; 0; 0];   % Stokes vector of the incident laser beam.

% medium
cloudCfg.clBase = 1000;   % cloud base height. [m]
cloudCfg.clDh = 10;   % geometrical thickness of each cloud layer. [m]
cloudCfg.nClLayers = 30;   % number of cloud layers.
cloudCfg.kappa = 9 * ones(1, cloudCfg.nClLayers);   % shape factor of modified gamma distribution.
cloudCfg.clExt = 2e-2 * ones(1, cloudCfg.nClLayers);   % cloud extinction. [m-1]
cloudCfg.clReff = 9 * ones(1, cloudCfg.nClLayers);   % cloud droplet effective radius. [um]

%% CPU MC
rng(1);
MSMCSim(lidarCfg, cloudCfg, 'visible', 'on', 'debug', false, 'nPhoton', 1e4);
