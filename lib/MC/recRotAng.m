function [rotationIn, rotationOut] = ...
    recRotAng(directionIn, directionOut, proto)
%RECROTANG Build the two Stokes rotations from vectors.
%
% The incident packet Stokes vector is expressed in its meridian-plane
% basis.  The returned signal is expressed in the receiver's fixed x-axis
% basis.  Constructing both rotations from projected unit vectors avoids
% ambiguity from supplementary-angle formulae near lidar backscatter.

n = cross(directionIn, directionOut, 2);
nNorm = sqrt(sum(n.^2, 2));
regular = nNorm > cast(1e-10, 'like', proto);
n = n ./ max(nNorm, cast(1e-10, 'like', proto));

eMeridian = meridianParallelAxis(directionIn, proto);
eDetector = fixedXAxisProjection(directionOut, proto);

pIn = cross(n, directionIn, 2);
pOut = cross(n, directionOut, 2);
pIn = normalizeRows(pIn, proto);
pOut = normalizeRows(pOut, proto);

% In the collinear limit the scattering plane is arbitrary.  Use the
% incident meridian axis as that arbitrary plane; this preserves the
% expected zero depolarization for single 180-degree scattering by spheres.
if any(~regular)
    pIn(~regular, :) = eMeridian(~regular, :);
    pOut(~regular, :) = eMeridian(~regular, :);
end

rotationIn = signedAxisAngle(eMeridian, pIn, directionIn);
rotationOut = signedAxisAngle(pOut, eDetector, directionOut);
end

function [rotationIn, rotationOut] = ...
    propagationRotationsVectorBasis(directionIn, directionOut, proto)
%PROPAGATIONROTATIONSVECTORBASIS Rotate between meridian/scattering planes.
n = cross(directionIn, directionOut, 2);
nNorm = sqrt(sum(n.^2, 2));
regular = nNorm > cast(1e-10, 'like', proto);
n = n ./ max(nNorm, cast(1e-10, 'like', proto));

eMeridianIn = meridianParallelAxis(directionIn, proto);
eMeridianOut = meridianParallelAxis(directionOut, proto);
pIn = normalizeRows(cross(n, directionIn, 2), proto);
pOut = normalizeRows(cross(n, directionOut, 2), proto);
if any(~regular)
    pIn(~regular, :) = eMeridianIn(~regular, :);
    pOut(~regular, :) = eMeridianOut(~regular, :);
end
rotationIn = signedAxisAngle(eMeridianIn, pIn, directionIn);
rotationOut = signedAxisAngle(pOut, eMeridianOut, directionOut);
end

function axis = meridianParallelAxis(direction, proto)
horizontal = sqrt(direction(:, 1).^2 + direction(:, 2).^2);
axis = zeros(size(direction), 'like', proto);
regular = horizontal > cast(1e-10, 'like', proto);
axis(regular, 1) = direction(regular, 3) .* ...
    direction(regular, 1) ./ horizontal(regular);
axis(regular, 2) = direction(regular, 3) .* ...
    direction(regular, 2) ./ horizontal(regular);
axis(regular, 3) = -horizontal(regular);
axis(~regular, 1) = 1;
end

function axis = fixedXAxisProjection(direction, proto)
axis = zeros(size(direction), 'like', proto);
axis(:, 1) = 1;
axis = axis - direction .* direction(:, 1);
axis = normalizeRows(axis, proto);
end

function rows = normalizeRows(rows, proto)
rowNorm = sqrt(sum(rows.^2, 2));
rows = rows ./ max(rowNorm, cast(1e-12, 'like', proto));
end

function angle = signedAxisAngle(axisFrom, axisTo, propagationDirection)
crossTerm = cross(axisFrom, axisTo, 2);
sineTerm = sum(propagationDirection .* crossTerm, 2);
cosineTerm = sum(axisFrom .* axisTo, 2);
angle = atan2(sineTerm, cosineTerm);
end

