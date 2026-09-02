function idx = sampleScatteringAngleGPU( ...
    cdf11, cdf12, intensity, qRot, layer, nLayers, proto)
%SAMPLESCATTERINGANGLEGPU Sample a scattering-grid index in MC_GPU.
u = rand(size(intensity), 'like', proto);
nTheta = size(cdf11, 2);
lastIndex = layer + (nTheta - 1) * nLayers;
lastCDF11 = cdf11(lastIndex);
lastCDF12 = cdf12(lastIndex);
lastCDF11 = lastCDF11(:);
lastCDF12 = lastCDF12(:);
target = u .* (intensity .* lastCDF11 + qRot .* lastCDF12);
lo = ones(size(intensity), 'like', proto);
hi = nTheta * ones(size(intensity), 'like', proto);
for k = 1:ceil(log2(double(nTheta)))
    mid = floor((lo + hi) / 2);
    opticalIndex = layer + (mid - 1) * nLayers;
    midCDF11 = cdf11(opticalIndex);
    midCDF12 = cdf12(opticalIndex);
    midCDF11 = midCDF11(:);
    midCDF12 = midCDF12(:);
    cdfMid = intensity .* midCDF11 + qRot .* midCDF12;
    moveLeft = cdfMid > target;
    hi(moveLeft) = mid(moveLeft);
    lo(~moveLeft) = mid(~moveLeft) + 1;
end
idx = min(lo, nTheta);
end

