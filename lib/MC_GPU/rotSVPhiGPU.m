function out = rotSVPhiGPU(in, phi)
%ROTSVPHIGPU Vectorized MC_GPU counterpart of rotSVPhi.
%
c = cos(2 * phi);
s = sin(2 * phi);
out = in;
out(:, 2) = in(:, 2) .* c + in(:, 3) .* s;
out(:, 3) = -in(:, 2) .* s + in(:, 3) .* c;
end

