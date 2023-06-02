function [samples] = rejectSampling(x, y, nSamples, varargin)
% REJECTSAMPLING rejection sampling.
%
% USAGE:
%    [samples] = rejectSampling(x, y, nSamples)
%
% INPUTS:
%    x: numeric
%    y: numeric
%    nSamples: numeric
%
% OUTPUTS:
%    samples: numeric
%
% HISTORY:
%    2023-05-24: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'x', @isnumeric);
addRequired(p, 'y', @isnumeric);
addRequired(p, 'nSamples', @isnumeric);

parse(p, x, y, nSamples, varargin{:});

sampleCounter = 0;
samples = [];
while (sampleCounter < nSamples)
    xSamples = rand(1, 1e3) * (x(end) - x(1)) + x(1);
    maxY = max(y);
    accepted = (rand(1, 1e3)) < (interp1(x, y, xSamples) / maxY);

    sampleCounter = sampleCounter + sum(accepted);
    samples = cat(2, samples, xSamples(accepted));
end

samples = samples(1:nSamples);

end