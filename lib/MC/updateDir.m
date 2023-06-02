function [outDir] = updateDir(inDir, phi, theta)
% UPDATEDIR description
%
% USAGE:
%    [outSV] = updateDir(inSV, phi, theta)
%
% INPUTS:
%    inSV, phi, theta
%
% OUTPUTS:
%    outSV
%
% EXAMPLE:
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

if (1 - abs(inDir(3))) < 1e-12
    outDir = [sin(theta) * cos(phi), ...
              sin(theta) * sin(phi), ...
              cos(theta) * ((inDir(3) >= 0) - 0.5) * 2];
else 
    temp = sqrt(1 - inDir(3)^2);
    outDir = [sin(theta) * (inDir(1) * inDir(3) * cos(phi) - inDir(2) * sin(phi)) / temp + inDir(1) * cos(theta), ...
              sin(theta) * (inDir(2) * inDir(3) * cos(phi) + inDir(1) * sin(phi)) / temp + inDir(2) * cos(theta), ...
              -sin(theta) * cos(phi) * temp + inDir(3) * cos(theta)];
end

end