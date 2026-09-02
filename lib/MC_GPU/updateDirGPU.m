function out = updateDirGPU(in, phi, theta, stateLayout)
%UPDATEDIRGPU Vectorized MC_GPU counterpart of updateDir.
%
% The two state-layout implementations are kept local so the public
% physical-process interface remains aligned with the CPU model.

if strcmp(stateLayout, 'packed')
    out = updateDirectionPackedImpl(in, phi, theta);
else
    out = updateDirectionLegacyImpl(in, phi, theta);
end
end

function out = updateDirectionLegacyImpl(in, phi, theta)
out = zeros(size(in), 'like', in);
vertical = (1 - abs(in(:, 3))) < 1e-12;
if any(vertical)
    inPhi = atan2(in(vertical, 2), in(vertical, 1));
    signZ = 2 * cast(in(vertical, 3) >= 0, 'like', in) - 1;
    out(vertical, 1) = sin(theta(vertical)) .* cos(inPhi + phi(vertical));
    out(vertical, 2) = sin(theta(vertical)) .* sin(inPhi + phi(vertical));
    out(vertical, 3) = cos(theta(vertical)) .* signZ;
end
regular = ~vertical;
if any(regular)
    temp = sqrt(1 - in(regular, 3).^2);
    out(regular, 1) = sin(theta(regular)) .* ...
        (in(regular, 1) .* in(regular, 3) .* cos(phi(regular)) - ...
        in(regular, 2) .* sin(phi(regular))) ./ temp + ...
        in(regular, 1) .* cos(theta(regular));
    out(regular, 2) = sin(theta(regular)) .* ...
        (in(regular, 2) .* in(regular, 3) .* cos(phi(regular)) + ...
        in(regular, 1) .* sin(phi(regular))) ./ temp + ...
        in(regular, 2) .* cos(theta(regular));
    out(regular, 3) = -sin(theta(regular)) .* cos(phi(regular)) .* temp + ...
        in(regular, 3) .* cos(theta(regular));
end
end

function out = updateDirectionPackedImpl(in, phi, theta)
out = zeros(size(in), 'like', in);
vertical = (1 - abs(in(:, 3))) < 1e-12;
verticalRows = reshape(find(vertical), [], 1);
inPhi = atan2(in(verticalRows, 2), in(verticalRows, 1));
signZ = 2 * cast(in(verticalRows, 3) >= 0, 'like', in) - 1;
out(verticalRows, 1) = sin(theta(verticalRows)) .* ...
    cos(inPhi + phi(verticalRows));
out(verticalRows, 2) = sin(theta(verticalRows)) .* ...
    sin(inPhi + phi(verticalRows));
out(verticalRows, 3) = cos(theta(verticalRows)) .* signZ;

regular = ~vertical;
regularRows = reshape(find(regular), [], 1);
temp = sqrt(1 - in(regularRows, 3).^2);
out(regularRows, 1) = sin(theta(regularRows)) .* ...
    (in(regularRows, 1) .* in(regularRows, 3) .* cos(phi(regularRows)) - ...
    in(regularRows, 2) .* sin(phi(regularRows))) ./ temp + ...
    in(regularRows, 1) .* cos(theta(regularRows));
out(regularRows, 2) = sin(theta(regularRows)) .* ...
    (in(regularRows, 2) .* in(regularRows, 3) .* cos(phi(regularRows)) + ...
    in(regularRows, 1) .* sin(phi(regularRows))) ./ temp + ...
    in(regularRows, 2) .* cos(theta(regularRows));
out(regularRows, 3) = -sin(theta(regularRows)) .* ...
    cos(phi(regularRows)) .* temp + ...
    in(regularRows, 3) .* cos(theta(regularRows));
end

