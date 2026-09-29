function [pOut, inside, ds, outLayer] = freeRunningGPU( ...
    pIn, direction, G, stateLayout)
%FREERUNNINGGPU Vectorized MC_GPU free-path propagation in a layered cloud.
%
% This is the gpuArray counterpart of freeRunning.  The packed path avoids
% scalar GPU-to-host synchronization while retaining the validated legacy
% path for compatibility.

if strcmp(stateLayout, 'packed')
    [pOut, inside, ds, outLayer] = ...
        freeRunningPackedImpl(pIn, direction, G);
else
    [pOut, inside, ds, outLayer] = ...
        freeRunningLegacyImpl(pIn, direction, G);
end
end

function [pOut, inside, ds, outLayer] = freeRunningLegacyImpl(pIn, direction, G)
tau = -log(rand(size(pIn, 1), 1, 'like', pIn));
ds = zeros(size(tau), 'like', tau);
pOut = nan(size(pIn), 'like', pIn);
outLayer = zeros(size(tau), 'like', tau);

inLayer = ceil((pIn(:, 3) - G.cloudBase) ./ G.layerThickness + 1e-8);
inLayer = max(1, min(G.nLayers, inLayer));

horizontal = abs(direction(:, 3)) < 1e-6;
if any(horizontal)
    horizontalLayer = inLayer(horizontal);
    ds(horizontal) = tau(horizontal) ./ G.muExt(horizontalLayer);
    pOut(horizontal, :) = pIn(horizontal, :) + ...
        direction(horizontal, :) .* ds(horizontal);
    outLayer(horizontal) = horizontalLayer;
end

oblique = ~horizontal;
if any(oblique)
    obliqueLayer = inLayer(oblique);
    layerBottom = G.cloudBase + ...
        (obliqueLayer - 1) .* G.layerThickness;
    tau1 = G.lowerTau(obliqueLayer) + ...
        (pIn(oblique, 3) - layerBottom) .* G.muExt(obliqueLayer);
    tau2 = tau1 + tau(oblique) .* direction(oblique, 3);

    destinationLayer = 1 + sum(tau2 >= G.cumulativeTau.', 2);
    remainsResolvable = destinationLayer <= G.nLayers;
    obliqueRows = find(oblique);
    resolvedRows = obliqueRows(remainsResolvable);
    if ~isempty(resolvedRows)
        resolvedLayer = destinationLayer(remainsResolvable);
        newZ = resolvedLayer .* G.layerThickness + G.cloudBase + ...
            (tau2(remainsResolvable) - G.cumulativeTau(resolvedLayer)) ./ ...
            G.muExt(resolvedLayer);
        ds(resolvedRows) = (newZ - pIn(resolvedRows, 3)) ./ ...
            direction(resolvedRows, 3);
        pOut(resolvedRows, :) = pIn(resolvedRows, :) + ...
            direction(resolvedRows, :) .* ds(resolvedRows);
        outLayer(resolvedRows) = resolvedLayer;
    end
end

base = G.cloudBase;
top = G.cloudTop;
candidate = outLayer >= 1 & outLayer <= G.nLayers & isfinite(pOut(:, 3));
inside = false(size(candidate), 'like', candidate);
if any(candidate)
    rows = find(candidate);
    bound = 4 ./ G.muExt(outLayer(rows));
    inside(rows) = abs(pOut(rows, 1)) <= bound & ...
        abs(pOut(rows, 2)) <= bound & ...
        pOut(rows, 3) >= base & pOut(rows, 3) <= top;
    if isfield(G,'horizontallyHomogeneous') && G.horizontallyHomogeneous
        inside(rows)=pOut(rows,3)>=base & pOut(rows,3)<=top;
    end
end

end

function [pOut, inside, ds, outLayer] = freeRunningPackedImpl(pIn, direction, G)
% The packed path permits empty GPU index sets, avoiding scalar any()
% results that would otherwise synchronize the GPU with the host.
tau = -log(rand(size(pIn, 1), 1, 'like', pIn));
ds = zeros(size(tau), 'like', tau);
pOut = nan(size(pIn), 'like', pIn);
outLayer = zeros(size(tau), 'like', tau);

inLayer = ceil((pIn(:, 3) - G.cloudBase) ./ G.layerThickness + 1e-8);
inLayer = max(1, min(G.nLayers, inLayer));

horizontal = abs(direction(:, 3)) < 1e-6;
horizontalRows = reshape(find(horizontal), [], 1);
horizontalLayer = inLayer(horizontalRows);
horizontalLayer = horizontalLayer(:);
ds(horizontalRows) = tau(horizontalRows) ./ G.muExt(horizontalLayer);
pOut(horizontalRows, :) = pIn(horizontalRows, :) + ...
    direction(horizontalRows, :) .* ...
    repmat(reshape(ds(horizontalRows), [], 1), 1, 3);
outLayer(horizontalRows) = horizontalLayer;

oblique = ~horizontal;
obliqueRows = reshape(find(oblique), [], 1);
obliqueLayer = inLayer(obliqueRows);
obliqueLayer = obliqueLayer(:);
layerBottom = G.cloudBase + ...
    (obliqueLayer - 1) .* G.layerThickness;
tau1 = G.lowerTau(obliqueLayer) + ...
    (pIn(obliqueRows, 3) - layerBottom) .* G.muExt(obliqueLayer);
tau2 = tau1 + tau(obliqueRows) .* direction(obliqueRows, 3);

destinationLayer = 1 + sum(tau2 >= G.cumulativeTau.', 2);
remainsResolvable = destinationLayer <= G.nLayers;
resolvedLocal = reshape(find(remainsResolvable), [], 1);
resolvedRows = obliqueRows(resolvedLocal);
resolvedLayer = destinationLayer(resolvedLocal);
resolvedLayer = resolvedLayer(:);
newZ = resolvedLayer .* G.layerThickness + G.cloudBase + ...
    (tau2(resolvedLocal) - G.cumulativeTau(resolvedLayer)) ./ ...
    G.muExt(resolvedLayer);
ds(resolvedRows) = (newZ - pIn(resolvedRows, 3)) ./ ...
    direction(resolvedRows, 3);
pOut(resolvedRows, :) = pIn(resolvedRows, :) + ...
    direction(resolvedRows, :) .* ...
    repmat(reshape(ds(resolvedRows), [], 1), 1, 3);
outLayer(resolvedRows) = resolvedLayer;

candidate = outLayer >= 1 & outLayer <= G.nLayers & isfinite(pOut(:, 3));
inside = false(size(candidate), 'like', candidate);
rows = reshape(find(candidate), [], 1);
bound = 4 ./ G.muExt(outLayer(rows));
inside(rows) = abs(pOut(rows, 1)) <= bound & ...
    abs(pOut(rows, 2)) <= bound & ...
    pOut(rows, 3) >= G.cloudBase & pOut(rows, 3) <= G.cloudTop;
if isfield(G,'horizontallyHomogeneous') && G.horizontallyHomogeneous
    inside(rows)=pOut(rows,3)>=G.cloudBase & pOut(rows,3)<=G.cloudTop;
end
end

