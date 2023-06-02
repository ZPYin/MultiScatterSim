function [samples] = inverseSampling(x, y, nSamples, varargin)
% INVERSESAMPLING inverse sampling.
%
% USAGE:
%    [samples] = inverseSampling(x, y, nSamples)
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

x = reshape(x, 1, length(x));
y = reshape(y, 1, length(y));
cumsumY = cumsum([0, diff(x)] .* y);
xRand = rand(1, nSamples);
samples = NaN(1, nSamples);
for iSample = 1:nSamples
    samples(iSample) = x(find(cumsumY > (xRand(iSample) * cumsumY(end)), 1));
end

end