function [SVReturn, diag] = MSMCSimGPUBackend(lidarCfg, cloudCfg, varargin)
%MSMCSIMGPUBackend Internal batched GPU transport implementation in MC_GPU.
%
% [SVReturn, diag] = MSMCSimGPU(lidarCfg, cloudCfg, ...)
%
% The implementation supports vertically layered clouds with one optical
% state per layer.  Unsupported configurations fall back to MSMCSim by
% default.  Photon states are advanced in batches with an active-photon
% mask, while the scattering CDF is searched with a per-photon binary search
% that does not materialize an NPhoton-by-NTheta matrix.
%
% GPU-specific options:
%   batchSize                  photons per GPU batch (default 131072)
%   precision                  'single' or 'double' (default 'double')
%   seed                       GPU RNG seed (default random)
%   maxScatterOrder            safety limit per photon (default 256)
%   fallbackToCPU              use MSMCSim if GPU path unsupported (true)
%   enableCUDAForwardCompatibility  support newer GPUs on older MATLAB
%   preparedOptics             output of prepareMSMCOptics
%   diagnosticsLevel           'full' (default) or 'minimal'
%   stateLayout                'legacy' (default) or 'packed'

tFunction = tic;

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'lidarCfg', @isstruct);
addRequired(p, 'cloudCfg', @isstruct);
addParameter(p, 'nPhoton', 1e4, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'visible', 'on', @ischar);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'useRussianRoulette', true, @islogical);
addParameter(p, 'energyThreshold', 0.01, @isnumeric);
addParameter(p, 'useLocalAtten', false, @islogical);
addParameter(p, 'scatMatrixFile', '', @ischar);
addParameter(p, 'angleSampler', 'fast', @ischar);
addParameter(p, 'preparedOptics', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'batchSize', 131072, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'precision', 'double', @ischar);
addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'maxScatterOrder', 256, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'fallbackToCPU', true, @islogical);
addParameter(p, 'enableCUDAForwardCompatibility', true, @islogical);
addParameter(p, 'diagnosticsLevel', 'full', @ischar);
addParameter(p, 'stateLayout', 'legacy', @ischar);
parse(p, lidarCfg, cloudCfg, varargin{:});

precision = validatestring(lower(p.Results.precision), {'single', 'double'});
angleSampler = validatestring(lower(p.Results.angleSampler), {'fast', 'inverse'});
diagnosticsLevelRequested = validatestring( ...
    lower(p.Results.diagnosticsLevel), {'full', 'minimal'});
stateLayoutRequested = validatestring( ...
    lower(p.Results.stateLayout), {'legacy', 'packed'});
stateLayout = stateLayoutRequested;
diagnosticsLevel = diagnosticsLevelRequested;
if strcmp(diagnosticsLevel, 'minimal') && strcmp(p.Results.visible, 'on')
    warning('MSMCSimGPU:MinimalDiagnosticsVisible', ...
        ['diagnosticsLevel="minimal" cannot build the diagnostic plots; ' ...
         'collecting full diagnostics because visible="on".']);
    diagnosticsLevel = 'full';
end
collectFullDiagnostics = strcmp(diagnosticsLevel, 'full');
nPhotons = floor(double(p.Results.nPhoton));
batchSize = min(nPhotons, floor(double(p.Results.batchSize)));
maxScatterOrder = floor(double(p.Results.maxScatterOrder));
maxOrderEstimate = 50;

unsupportedReason = '';
if numel(cloudCfg.clDh) ~= 1
    unsupportedReason = 'GPU prototype currently requires scalar clDh.';
elseif ~strcmp(angleSampler, 'fast')
    unsupportedReason = 'GPU prototype currently supports angleSampler="fast".';
end

if ~isempty(unsupportedReason)
    [SVReturn, diag] = fallbackOrError(unsupportedReason);
    return;
end

tGPUSelect = tic;
[gpuOK, gpuReason, gpuInfo] = selectGPU(p.Results.enableCUDAForwardCompatibility);
gpuSelectSeconds = toc(tGPUSelect);
if ~gpuOK
    [SVReturn, diag] = fallbackOrError(gpuReason);
    return;
end

unmatchedArgs = p.Unmatched;
extraArgs = reshape([fieldnames(unmatchedArgs).'; struct2cell(unmatchedArgs).'], 1, []);

tPrepare = tic;
if isempty(p.Results.preparedOptics)
    optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
        'scatMatrixFile', p.Results.scatMatrixFile, extraArgs{:});
    opticsSource = 'computed';
else
    optics = p.Results.preparedOptics;
    validatePreparedOptics(optics, lidarCfg, cloudCfg);
    opticsSource = 'prepared';
end
prepareSeconds = toc(tPrepare);

if numel(optics.muExt) ~= cloudCfg.nClLayers
    [SVReturn, diag] = fallbackOrError( ...
        'Prepared optical data layer count does not match cloudCfg.', optics);
    return;
end

if isempty(p.Results.seed)
    seed = randi([1, intmax('int32')], 1, 1);
else
    seed = mod(floor(double(p.Results.seed)), double(intmax('int32')));
    if seed < 1, seed = 1; end
end
gpurng(seed);

tUpload = tic;
G = moveOpticsToGPU(optics, precision);
G.SBeam = gpuArray(cast(lidarCfg.SBeam(:).', precision));
G.cloudBase = gpuArray(cast(cloudCfg.clBase, precision));
G.layerThickness = gpuArray(cast(cloudCfg.clDh, precision));
G.cloudTop = gpuArray(cast( ...
    cloudCfg.clBase + cloudCfg.nClLayers * cloudCfg.clDh, precision));
G.cumulativeTau = cumsum(G.muExt .* G.layerThickness);
G.lowerTau = zeros(size(G.cumulativeTau), 'like', G.cumulativeTau);
if G.nLayers > 1
    G.lowerTau(2:end) = G.cumulativeTau(1:end-1);
end
wait(gpuInfo);
uploadSeconds = toc(tUpload);

SVReturn = zeros(4, lidarCfg.nBins);
if collectFullDiagnostics
    SVReturnLayer = zeros(4, lidarCfg.nBins, cloudCfg.nClLayers);
    SVByOrder = zeros(4, lidarCfg.nBins, maxOrderEstimate);
    orderPhotonCount = zeros(1, maxOrderEstimate);
    layerScatterCount = zeros(1, cloudCfg.nClLayers);
    nScasHist = zeros(1, maxScatterOrder + 1);
else
    SVReturnLayer = [];
    SVByOrder = [];
    orderPhotonCount = [];
    layerScatterCount = [];
    nScasHist = [];
end
truncatedPhotons = 0;

batchContext = struct();
batchContext.lidarCfg = lidarCfg;
batchContext.cloudCfg = cloudCfg;
batchContext.options = p.Results;
batchContext.stateLayout = stateLayout;
batchContext.collectFullDiagnostics = collectFullDiagnostics;
batchContext.maxScatterOrder = maxScatterOrder;
batchContext.maxOrderEstimate = maxOrderEstimate;
batchContext.precision = precision;

nBatches = ceil(nPhotons / batchSize);
batchSeconds = zeros(1, nBatches);
batchCollectSeconds = zeros(1, nBatches);
tTransport = tic;
for iBatch = 1:nBatches
    firstPhoton = (iBatch - 1) * batchSize + 1;
    nThis = min(batchSize, nPhotons - firstPhoton + 1);
    tBatch = tic;
    batch = simulateLayeredBatchGPU(nThis, G, batchContext);
    wait(gpuInfo);
    batchSeconds(iBatch) = toc(tBatch);

    tCollect = tic;
    SVReturn = SVReturn + double(gather(batch.SVReturn));
    truncatedPhotons = truncatedPhotons + double(gather(batch.truncatedPhotons));
    if collectFullDiagnostics
        SVReturnLayer = SVReturnLayer + double(gather(batch.SVReturnLayer));
        SVByOrder = SVByOrder + double(gather(batch.SVByOrder));
        orderPhotonCount = orderPhotonCount + ...
            double(gather(batch.orderPhotonCount));
        layerScatterCount = layerScatterCount + ...
            double(gather(batch.layerScatterCount)).';

        if isfield(batch, 'nScasHist')
            nScasHist = nScasHist + double(gather(batch.nScasHist)).';
        else
            nScasCPU = double(gather(batch.nScas));
            h = histcounts(nScasCPU, -0.5:1:(maxScatterOrder + 0.5));
            nScasHist = nScasHist + h;
        end
    end
    batchCollectSeconds(iBatch) = toc(tCollect);

    if p.Results.debug
        fprintf('GPU batch %d/%d: %d photons, %.3f s\n', ...
            iBatch, nBatches, nThis, batchSeconds(iBatch));
    end
end
transportSeconds = toc(tTransport);

tPostprocess = tic;
diag = struct();
diag.backend = 'gpu';
diag.gpuName = gpuInfo.Name;
diag.computeCapability = gpuInfo.ComputeCapability;
diag.precision = precision;
diag.seed = seed;
diag.batchSize = batchSize;
diag.nBatches = nBatches;
diag.batchSeconds = batchSeconds;
diag.batchCollectSeconds = batchCollectSeconds;
diag.gpuComputeSeconds = sum(batchSeconds);
diag.hostCollectSeconds = sum(batchCollectSeconds);
diag.gpuSelectSeconds = gpuSelectSeconds;
diag.prepareSeconds = prepareSeconds;
diag.uploadSeconds = uploadSeconds;
diag.transportSeconds = transportSeconds;
diag.opticsSource = opticsSource;
diag.angleSampler = angleSampler;
diag.diagnosticsLevelRequested = diagnosticsLevelRequested;
diag.diagnosticsLevel = diagnosticsLevel;
diag.stateLayoutRequested = stateLayoutRequested;
diag.stateLayout = stateLayout;
diag.muExt = optics.muExt;
diag.muSca = optics.muSca;
diag.albedo = optics.albedo;
diag.maxOrderEstimate = maxOrderEstimate;
diag.truncatedPhotons = truncatedPhotons;

if collectFullDiagnostics
    diag.SVReturnLayer = SVReturnLayer;
    diag.layerI = squeeze(sum(SVReturnLayer(1, :, :), 2));
    diag.layerScatterCount = layerScatterCount;
    diag.SVByOrder = SVByOrder;
    diag.orderPhotonCount = orderPhotonCount;
    diag.nScasHistEdges = 0:(maxScatterOrder + 1);
    diag.nScasHist = nScasHist;

    IByOrder = squeeze(sum(SVByOrder(1, :, :), 2));
    QByOrder = squeeze(sum(SVByOrder(2, :, :), 2));
    diag.I_by_order = IByOrder;
    diag.Q_by_order = QByOrder;
    diag.Pperp_by_order = (IByOrder - QByOrder) / 2;
    diag.Ppar_by_order = (IByOrder + QByOrder) / 2;
    diag.delta_by_order = diag.Pperp_by_order ./ ...
        max(diag.Ppar_by_order, eps);
    diag.f_by_order = IByOrder ./ max(sum(IByOrder), eps);
    diag.lastActiveOrder = find(IByOrder > 0, 1, 'last');
    if isempty(diag.lastActiveOrder), diag.lastActiveOrder = 0; end
else
    diag.SVReturnLayer = [];
    diag.layerI = [];
    diag.layerScatterCount = [];
    diag.SVByOrder = [];
    diag.orderPhotonCount = [];
    diag.nScasHistEdges = [];
    diag.nScasHist = [];
    diag.I_by_order = [];
    diag.Q_by_order = [];
    diag.Pperp_by_order = [];
    diag.Ppar_by_order = [];
    diag.delta_by_order = [];
    diag.f_by_order = [];
    diag.lastActiveOrder = [];
end

diag.postprocessSeconds = toc(tPostprocess);
diag.totalSeconds = toc(tFunction);
diag.otherSeconds = max(0, diag.totalSeconds - ...
    (diag.gpuSelectSeconds + diag.prepareSeconds + ...
     diag.uploadSeconds + diag.transportSeconds + diag.postprocessSeconds));

if strcmp(p.Results.visible, 'on')
    % CPU and GPU use the same diagnostic quantities and axis conventions.
    counts = double(diag.nScasHist(:).');
    lastNonzero = find(counts > 0, 1, 'last');
    if isempty(lastNonzero)
        counts = 0;
    else
        counts = counts(1:lastNonzero);
    end
    orders = 0:(numel(counts) - 1);

    distanceAxis = (1:lidarCfg.nBins) * lidarCfg.binWidth * 0.15;
    IPlot = abs(double(SVReturn(1, :)));
    QPlot = abs(double(SVReturn(2, :)));
    IPlot(IPlot == 0) = NaN;
    QPlot(QPlot == 0) = NaN;

    denominator = double(SVReturn(1, :) + SVReturn(2, :));
    signalFloor = max(max(abs(denominator)) * 1e-12, eps);
    validDelta = isfinite(denominator) & denominator > signalFloor;
    delta = NaN(1, lidarCfg.nBins);
    delta(validDelta) = double(SVReturn(1, validDelta) - ...
        SVReturn(2, validDelta)) ./ denominator(validDelta);

    figure('Color', 'w', 'Name', 'MSMC diagnostic results', ...
        'NumberTitle', 'off');

    subplot(3, 1, 1);
    bar(orders, counts, 1, 'FaceColor', [0.25, 0.60, 0.82], ...
        'EdgeColor', [0.15, 0.15, 0.15]);
    xlim([-0.5, max(10, orders(end) + 1)]);
    title('Photon scattering order distribution', 'FontSize', 12);
    xlabel('Scattering order', 'FontSize', 10);
    ylabel('Photon count', 'FontSize', 10);
    grid on;

    subplot(3, 1, 2);
    semilogy(distanceAxis, IPlot, 'LineWidth', 1.0, ...
        'DisplayName', 'I');
    hold on;
    semilogy(distanceAxis, QPlot, 'LineWidth', 1.0, ...
        'DisplayName', '|Q|');
    hold off;
    xlim([0, distanceAxis(end)]);
    title('Returned Stokes components', 'FontSize', 12);
    xlabel('Distance (m)', 'FontSize', 10);
    ylabel('Signal intensity', 'FontSize', 10);
    legend('show', 'Location', 'best');
    grid on;

    subplot(3, 1, 3);
    plot(distanceAxis, delta, 'LineWidth', 1.5);
    xlim([0, distanceAxis(end)]);
    ylim([0, 1]);
    title('Linear depolarization ratio', 'FontSize', 12);
    xlabel('Distance (m)', 'FontSize', 10);
    ylabel('Depolarization ratio', 'FontSize', 10);
    grid on;
end

    function [SV, d] = fallbackOrError(reason, prepared)
        if nargin < 2
            prepared = p.Results.preparedOptics;
        end
        if ~p.Results.fallbackToCPU
            error('MSMCSimGPU:Unsupported', '%s', reason);
        end
        cpuArgs = { ...
            'nPhoton', p.Results.nPhoton, ...
            'visible', p.Results.visible, ...
            'debug', p.Results.debug, ...
            'useRussianRoulette', p.Results.useRussianRoulette, ...
            'energyThreshold', p.Results.energyThreshold, ...
            'useLocalAtten', p.Results.useLocalAtten, ...
            'scatMatrixFile', p.Results.scatMatrixFile, ...
            'angleSampler', p.Results.angleSampler};
        if ~isempty(prepared)
            cpuArgs = [cpuArgs, {'preparedOptics', prepared}];
        end
        unmatched = reshape([fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'], 1, []);
        [SV, d] = MSMCSim(lidarCfg, cloudCfg, cpuArgs{:}, unmatched{:});
        d.backend = 'cpu-fallback';
        d.fallbackReason = reason;
    end
end

function G = moveOpticsToGPU(optics, precision)
toGPU = @(x) gpuArray(cast(x, precision));
G.theta = toGPU(optics.theta(:));
G.S11 = toGPU(optics.S11);
G.S12 = toGPU(optics.S12);
G.S33 = toGPU(optics.S33);
G.S34 = toGPU(optics.S34);
G.cdfS11 = toGPU(optics.cdfS11);
G.cdfS12 = toGPU(optics.cdfS12);
G.muExt = toGPU(optics.muExt(:));
G.albedo = toGPU(optics.albedo(:));
G.nLayers = numel(optics.muExt);
end

function [ok, reason, device] = selectGPU(enableForwardCompatibility)
ok = false;
reason = '';
device = [];
try
    if gpuDeviceCount('available') == 0 && enableForwardCompatibility
        parallel.gpu.enableCUDAForwardCompatibility(true);
    end
    if gpuDeviceCount('available') == 0
        reason = 'No MATLAB-compatible GPU is available.';
        return;
    end
    device = gpuDevice;
    ok = true;
catch ME
    reason = sprintf('GPU initialization failed: %s', ME.message);
end
end

function validatePreparedOptics(optics, lidarCfg, cloudCfg)
required = {'cloudScatMatrix', 'muExt', 'muSca', 'albedo', ...
    'S11', 'S12', 'S33', 'S34', 'theta', 'cdfS11', 'cdfS12'};
for iField = 1:numel(required)
    if ~isfield(optics, required{iField})
        error('MSMCSimGPU:InvalidPreparedOptics', ...
            'preparedOptics is missing field "%s".', required{iField});
    end
end
if isfield(optics, 'lambda') && optics.lambda ~= lidarCfg.lambda
    error('MSMCSimGPU:PreparedOpticsWavelength', ...
        'Prepared wavelength does not match lidarCfg.lambda.');
end
if numel(optics.muExt) ~= cloudCfg.nClLayers
    error('MSMCSimGPU:PreparedOpticsLayers', ...
        'Prepared layer count does not match cloudCfg.nClLayers.');
end
end
