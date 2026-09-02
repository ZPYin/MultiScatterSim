function [scaAng, rotAng1, rotAng2] = ...
    scaRotAngGPU(in, out, stateLayout)
%SCAROTANGGPU Vectorized MC_GPU counterpart of scaRotAng.
%
% It returns the scattering angle and the two Stokes reference-plane
% rotation angles for every photon direction pair.

if strcmp(stateLayout, 'packed')
    [scaAng, rotAng1, rotAng2] = ...
        scatteringRotationAnglesPackedImpl(in, out);
else
    [scaAng, rotAng1, rotAng2] = ...
        scatteringRotationAnglesLegacyImpl(in, out);
end
end

function [scaAng, rotAng1, rotAng2] = scatteringRotationAnglesLegacyImpl(in, out)
dotProduct = max(-1, min(1, sum(in .* out, 2)));
scaAng = acos(dotProduct);
rotAng1 = zeros(size(scaAng), 'like', scaAng);
rotAng2 = zeros(size(scaAng), 'like', scaAng);

parallel = abs(sin(scaAng)) < 1e-12;
inVertical = ~parallel & (1 - abs(in(:, 3))) <= 1e-12;
outVertical = ~parallel & ~inVertical & (1 - abs(out(:, 3))) <= 1e-12;
general = ~parallel & ~inVertical & ~outVertical;

inPhi = atan2(in(:, 2), in(:, 1));
outPhi = atan2(out(:, 2), out(:, 1));
if any(inVertical)
    rotAng1(inVertical) = acos(max(-1, min(1, ...
        in(inVertical, 3) .* cos(outPhi(inVertical) - inPhi(inVertical)))));
    rotAng2(inVertical) = acos(max(-1, min(1, in(inVertical, 3))));
end
if any(outVertical)
    rotAng1(outVertical) = acos(max(-1, min(1, ...
        out(outVertical, 3) .* cos(outPhi(outVertical) - inPhi(outVertical)))));
    rotAng2(outVertical) = acos(max(-1, min(1, out(outVertical, 3))));
end
if any(general)
    inTheta = acos(max(-1, min(1, in(general, 3))));
    outTheta = acos(max(-1, min(1, out(general, 3))));
    sa = scaAng(general);
    c1 = (out(general, 3) - in(general, 3) .* cos(sa)) ./ ...
        (sin(inTheta) .* sin(sa));
    c2 = (in(general, 3) - out(general, 3) .* cos(sa)) ./ ...
        (sin(outTheta) .* sin(sa));
    rotAng1(general) = acos(max(-1, min(1, c1)));
    rotAng2(general) = acos(max(-1, min(1, c2)));
end

rotationSign = in(:, 1) .* out(:, 2) - in(:, 2) .* out(:, 1);
% scaRotAng returns before quadrant correction for its parallel and
% vertical special cases.  Apply the sign only to the general branch to
% preserve the CPU reference behavior exactly.
negative = general & rotationSign < 0;
rotAng1(negative) = -rotAng1(negative);
rotAng2(negative) = -rotAng2(negative);
end

function [scaAng, rotAng1, rotAng2] = ...
    scatteringRotationAnglesPackedImpl(in, out)
dotProduct = max(-1, min(1, sum(in .* out, 2)));
scaAng = acos(dotProduct);
rotAng1 = zeros(size(scaAng), 'like', scaAng);
rotAng2 = zeros(size(scaAng), 'like', scaAng);

parallel = abs(sin(scaAng)) < 1e-12;
inVertical = ~parallel & (1 - abs(in(:, 3))) <= 1e-12;
outVertical = ~parallel & ~inVertical & (1 - abs(out(:, 3))) <= 1e-12;
general = ~parallel & ~inVertical & ~outVertical;

inPhi = atan2(in(:, 2), in(:, 1));
outPhi = atan2(out(:, 2), out(:, 1));
inVerticalRows = reshape(find(inVertical), [], 1);
outVerticalRows = reshape(find(outVertical), [], 1);
generalRows = reshape(find(general), [], 1);
rotAng1(inVerticalRows) = acos(max(-1, min(1, ...
    in(inVerticalRows, 3) .* ...
    cos(outPhi(inVerticalRows) - inPhi(inVerticalRows)))));
rotAng2(inVerticalRows) = acos(max(-1, min(1, in(inVerticalRows, 3))));
rotAng1(outVerticalRows) = acos(max(-1, min(1, ...
    out(outVerticalRows, 3) .* ...
    cos(outPhi(outVerticalRows) - inPhi(outVerticalRows)))));
rotAng2(outVerticalRows) = acos(max(-1, min(1, out(outVerticalRows, 3))));

inTheta = acos(max(-1, min(1, in(generalRows, 3))));
outTheta = acos(max(-1, min(1, out(generalRows, 3))));
sa = scaAng(generalRows);
c1 = (out(generalRows, 3) - in(generalRows, 3) .* cos(sa)) ./ ...
    (sin(inTheta) .* sin(sa));
c2 = (in(generalRows, 3) - out(generalRows, 3) .* cos(sa)) ./ ...
    (sin(outTheta) .* sin(sa));
rotAng1(generalRows) = acos(max(-1, min(1, c1)));
rotAng2(generalRows) = acos(max(-1, min(1, c2)));

rotationSign = in(:, 1) .* out(:, 2) - in(:, 2) .* out(:, 1);
negative = general & rotationSign < 0;
negativeRows = reshape(find(negative), [], 1);
rotAng1(negativeRows) = -rotAng1(negativeRows);
rotAng2(negativeRows) = -rotAng2(negativeRows);
end

