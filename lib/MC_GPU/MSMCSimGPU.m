function [SVReturn, diag] = MSMCSimGPU(lidarCfg, cloudCfg, varargin)
%MSMCSIMGPU Public entry point for the GPU Monte Carlo backend.
%
% [SVReturn, diag] = MSMCSimGPU(lidarCfg, cloudCfg, ...)
%
% The GPU implementation lives in lib/MC_GPU.  This facade preserves the
% existing public function name.  Use setupMSToolbox, or recursively add
% the complete lib directory, so the shared CPU/optics code is also visible.

persistent gpuBackendReady
if isempty(gpuBackendReady)
    backendDir = fileparts(mfilename('fullpath'));
    backendFile = fullfile(backendDir, 'MSMCSimGPUBackend.m');
    if exist(backendFile, 'file') ~= 2
        error('MSMCSimGPU:MissingBackend', ...
            'GPU backend implementation was not found: %s', backendFile);
    end
    if isempty(which('MSMCSimGPUBackend'))
        addpath(backendDir);
    end
    gpuBackendReady = true;
end

[SVReturn, diag] = MSMCSimGPUBackend(lidarCfg, cloudCfg, varargin{:});
end
