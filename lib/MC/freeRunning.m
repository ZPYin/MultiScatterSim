function [pOutPos, isInMedium, ds, outClLayer] = freeRunning(pInPos, pInDir, cloudExt, dClThickness, cloudBaseHeight, cloudBoundX, cloudBoundY)
% FREERUNNING calculate the next move of photon
%
% USAGE:
%    [pOutPos, isInMedium, ds] = freeRunning(pInPos, pInDir, cloudExt, dClThickness, cloudBaseHeight, cloudBoundX, cloudBoundY)
%
% INPUTS:
%    pInPos: numeric
%        position of incident photon.
%    pInDir: numeric
%        photon move direction.
%    cloudExt: numeric
%        cloud layer extinction. (m-1)
%    dClThickness: numeric
%        cloud layer thickness in vertical. (m)
%    cloudBaseHeight: numeric
%        cloud base height. (m)
%    cloudBoundX: numeric
%        cloud boundary in x axis. (m)
%    cloudBoundY: numeric
%        cloud boundary in y axis. (m)
%
% OUTPUTS:
%    pOutPos: numeric
%        photon position after the move.
%    isInMedium: logical
%        logical variable to determine whether photon is still in the medium.
%    ds: numeric
%        moving length.
%    outClLayer: nmumeric
%        index of cloud layer where the photon is located.
%
% HISTORY:
%    2023-06-02: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

pOutPos = NaN(1, 3);
isInMedium = true;
ds = 0;

tau0 = - log(rand(1, 1));

inClLayer = ceil((pInPos(3) - cloudBaseHeight) / dClThickness + 1e-8);
if abs(pInDir(3)) < 1e-6

    % horizontal running
    ds = tau0 / cloudExt(inClLayer);
    pOutPos = pInPos + pInDir * ds;
    outClLayer = inClLayer;

else

    cumsumTau = cumsum(cloudExt .* dClThickness);

    % slant running
    if inClLayer > 1
        tau1 = (sum(cloudExt(1:(inClLayer - 1)) * dClThickness)) + (pInPos(3) - cloudBaseHeight - (inClLayer - 1) * dClThickness) * cloudExt(inClLayer);
    else
        tau1 = (pInPos(3) - cloudBaseHeight - (inClLayer - 1) * dClThickness) * cloudExt(inClLayer);
    end
    tau2 = tau1 + tau0 * pInDir(3);
    idx = find(tau2 < cumsumTau, 1);
    if isempty(idx)
        outClLayer = NaN;
    else
        outClLayer = idx;
        pOutPos(3) = idx * dClThickness + cloudBaseHeight +(tau2 - cumsumTau(idx)) / cloudExt(idx);
        ds = (pOutPos(3) - pInPos(3)) / pInDir(3);

        pOutPos = pInPos + pInDir * ds;
    end

end

% determine whether photon is still in medium
if isnan(outClLayer)
    isInMedium = false;
elseif (abs(pOutPos(1)) > cloudBoundX(outClLayer)) || (abs(pOutPos(2)) > cloudBoundY(outClLayer))
    isInMedium = false;
elseif (pOutPos(3) > (cloudBaseHeight + dClThickness * length(cloudExt))) || (pOutPos(3) < cloudBaseHeight)
    isInMedium = false;
end

end