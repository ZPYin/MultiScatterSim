function [optics, info] = prepareMSMCOpticsCached( ...
    lidarCfg, cloudCfg, varargin)
%PREPAREMSMCOPTICSCACHED Cache shared Nd-scalable MSMC optical properties.
%
% [optics, info] = prepareMSMCOpticsCached(lidarCfg, cloudCfg, ...)
%
% For a fixed wavelength, Reff/kappa profile, scattering database, and
% normalized Nd(z) shape, changing the absolute Nd scale changes muExt and
% muSca linearly while leaving the normalized Mueller matrices and their
% CDFs unchanged.  This function stores one reference calculation and
% reconstructs other Nd amplitudes exactly by linear scaling.
%
% Options:
%   scatMatrixFile   Mie database used by prepareMSMCOptics
%   cacheDir         directory for versioned MAT cache entries
%   cacheNamespace   filename namespace (default 'msmc-optics-v1')
%   useDiskCache     reuse entries across MATLAB processes (default true)
%   useMemoryCache   retain the most recent entry in-process (default true)
%   cloud_distribution distribution passed to prepareMSMCOptics

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'lidarCfg', @isstruct);
addRequired(p, 'cloudCfg', @isstruct);
addParameter(p, 'scatMatrixFile', '', @ischar);
addParameter(p, 'cacheDir', '', @ischar);
addParameter(p, 'cacheNamespace', 'msmc-optics-v1', @ischar);
addParameter(p, 'useDiskCache', true, @islogical);
addParameter(p, 'useMemoryCache', true, @islogical);
addParameter(p, 'cloud_distribution', 'modified-gamma', @ischar);
parse(p, lidarCfg, cloudCfg, varargin{:});

tTotal = tic;
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
scatMatrixFile = p.Results.scatMatrixFile;
if isempty(scatMatrixFile)
    scatMatrixFile = fullfile(projectDir, 'data', ...
        'water_droplet_scattering_532nm_finer_validation_r_step_0p01um.mat');
end
assert(exist(scatMatrixFile, 'file') == 2, ...
    'Scattering database does not exist: %s', scatMatrixFile);

cacheDir = p.Results.cacheDir;
if isempty(cacheDir)
    cacheDir = fullfile(projectDir, 'tests', 'data', 'gpu_optics_cache');
end
cacheNamespace = regexprep(p.Results.cacheNamespace, ...
    '[^A-Za-z0-9_-]', '_');

extraArgs = reshape( ...
    [fieldnames(p.Unmatched).'; struct2cell(p.Unmatched).'], 1, []);

info = struct();
info.source = 'computed-uncached';
info.cacheKey = '';
info.cacheFile = '';
info.lookupSeconds = 0;
info.computeSeconds = 0;
info.saveSeconds = 0;
info.scale = 1;

if ~isfield(cloudCfg, 'N_d') || isempty(cloudCfg.N_d)
    tCompute = tic;
    optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
        'scatMatrixFile', scatMatrixFile, ...
        'cloud_distribution', p.Results.cloud_distribution, ...
        extraArgs{:});
    info.computeSeconds = toc(tCompute);
    info.totalSeconds = toc(tTotal);
    return;
end

nLayers = cloudCfg.nClLayers;
ndVector = expandLayerVector(cloudCfg.N_d, nLayers, 'N_d');
assert(all(isfinite(ndVector)) && all(ndVector > 0), ...
    'N_d must be finite and positive for scalable optical caching.');
ndShape = ndVector ./ ndVector(1);

fileInfo = dir(scatMatrixFile);
keyData = struct();
keyData.formatVersion = 1;
keyData.lambda = double(lidarCfg.lambda);
keyData.scatMatrixFile = char(scatMatrixFile);
keyData.scatMatrixBytes = double(fileInfo.bytes);
keyData.scatMatrixDatenum = double(fileInfo.datenum);
keyData.cloudDistribution = p.Results.cloud_distribution;
keyData.nLayers = double(nLayers);
keyData.kappa = double(cloudCfg.kappa(:).');
keyData.reff = double(cloudCfg.clReff(:).');
keyData.ndShape = double(ndShape(:).');
keyData.extraParameters = p.Unmatched;

cacheKey = sha256Text(jsonencode(keyData));
cacheFile = fullfile(cacheDir, sprintf('%s_%s.mat', ...
    cacheNamespace, cacheKey(1:24)));
info.cacheKey = cacheKey;
info.cacheFile = cacheFile;

persistent memoryCacheKey memoryCacheEntry
cacheEntry = [];
tLookup = tic;
if p.Results.useMemoryCache && ...
        ~isempty(memoryCacheKey) && strcmp(memoryCacheKey, cacheKey)
    cacheEntry = memoryCacheEntry;
    info.source = 'memory';
elseif p.Results.useDiskCache && exist(cacheFile, 'file') == 2
    try
        loaded = load(cacheFile, 'cacheEntry');
        if isfield(loaded, 'cacheEntry') && ...
                isValidCacheEntry(loaded.cacheEntry, keyData, nLayers)
            cacheEntry = loaded.cacheEntry;
            info.source = 'disk';
        else
            warning('prepareMSMCOpticsCached:InvalidEntry', ...
                'Ignoring incompatible optical cache entry: %s', cacheFile);
        end
    catch ME
        warning('prepareMSMCOpticsCached:UnreadableEntry', '%s', ...
            sprintf('Ignoring unreadable optical cache entry %s: %s', ...
            cacheFile, ME.message));
    end
end
info.lookupSeconds = toc(tLookup);

if isempty(cacheEntry)
    tCompute = tic;
    referenceOptics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
        'scatMatrixFile', scatMatrixFile, ...
        'cloud_distribution', p.Results.cloud_distribution, ...
        extraArgs{:});
    info.computeSeconds = toc(tCompute);

    referenceNd = expandLayerVector( ...
        cloudCfg.N_d, numel(referenceOptics.muExt), 'N_d');
    cacheEntry = struct();
    cacheEntry.formatVersion = 1;
    cacheEntry.keyData = keyData;
    cacheEntry.referenceNd = referenceNd;
    cacheEntry.optics = referenceOptics;
    cacheEntry.createdAt = char(datetime( ...
        'now', 'Format', 'yyyy-MM-dd HH:mm:ss Z'));
    info.source = 'computed';

    if p.Results.useDiskCache
        if ~exist(cacheDir, 'dir')
            mkdir(cacheDir);
        end
        tSave = tic;
        temporaryFile = [tempname(cacheDir) '.mat'];
        save(temporaryFile, 'cacheEntry', '-v7.3');
        [moveOK, moveMessage] = movefile(temporaryFile, cacheFile, 'f');
        if ~moveOK
            if exist(temporaryFile, 'file') == 2
                delete(temporaryFile);
            end
            error('prepareMSMCOpticsCached:CacheMoveFailed', ...
                'Could not commit optical cache entry: %s', moveMessage);
        end
        info.saveSeconds = toc(tSave);
    end
end

if p.Results.useMemoryCache
    memoryCacheKey = cacheKey;
    memoryCacheEntry = cacheEntry;
end

referenceNd = expandLayerVector( ...
    cacheEntry.referenceNd, numel(cacheEntry.optics.muExt), 'referenceNd');
targetNd = expandLayerVector( ...
    cloudCfg.N_d, numel(cacheEntry.optics.muExt), 'N_d');
scaleByLayer = targetNd ./ referenceNd;
tolerance = 32 * eps(max(1, max(abs(scaleByLayer))));
if max(abs(scaleByLayer - scaleByLayer(1))) > tolerance
    error('prepareMSMCOpticsCached:NonuniformScaling', ...
        ['Cached and requested Nd profiles do not differ by one scalar. ' ...
         'The cache key or profile normalization is inconsistent.']);
end

info.scale = scaleByLayer(1);
optics = scaleOptics(cacheEntry.optics, info.scale);
info.totalSeconds = toc(tTotal);
end

function values = expandLayerVector(values, nLayers, parameterName)
values = double(values(:));
if isscalar(values)
    values = repmat(values, nLayers, 1);
elseif numel(values) ~= nLayers
    error('prepareMSMCOpticsCached:LayerCount', ...
        '%s has %d values; expected 1 or %d.', ...
        parameterName, numel(values), nLayers);
end
end

function valid = isValidCacheEntry(entry, keyData, nLayers)
required = {'formatVersion', 'keyData', 'referenceNd', 'optics'};
valid = isstruct(entry) && all(isfield(entry, required)) && ...
    entry.formatVersion == 1 && isequaln(entry.keyData, keyData) && ...
    numel(entry.referenceNd) == nLayers && ...
    isfield(entry.optics, 'muExt') && ...
    numel(entry.optics.muExt) == nLayers;
end

function optics = scaleOptics(reference, scale)
optics = reference;
if scale == 1
    return;
end
optics.muExt = reference.muExt .* scale;
optics.muSca = reference.muSca .* scale;
optics.albedo = optics.muSca ./ optics.muExt;
optics.cloudScatMatrix.muExt = ...
    reference.cloudScatMatrix.muExt .* scale;
optics.cloudScatMatrix.muSca = ...
    reference.cloudScatMatrix.muSca .* scale;
end

function hash = sha256Text(value)
digest = java.security.MessageDigest.getInstance('SHA-256');
digest.update(uint8(unicode2native(value, 'UTF-8')));
bytes = typecast(digest.digest(), 'uint8');
hash = lower(reshape(dec2hex(bytes, 2).', 1, []));
end
