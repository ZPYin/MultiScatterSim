function [outDir] = updateDir(inDir, phi, theta)
% UPDATEDIR Update the transmit direction.
%
% USAGE:
%    [outSV] = updateDir(inSV, phi, theta)
%
% INPUTS:
%    inSV: vector
%        Stokes Vector of the incident photon.
%    phi: numeric
%        relative azimuth angle. (radians)
%    theta: numeric
%        scattering angle. (radians)
%
% OUTPUTS:
%    outSV: vector
%        Stokes Vector of the scattered photon.
%
% REFERENCES:
%    [1] Hu Y-X, Winker D, Yang P, Baum B, Poole L, Vann L. Identification of cloud phase from PICASSO-CENA lidar depolarization: a multiple scattering sensitivity study. Journal of Quantitative Spectroscopy and Radiative Transfer. 2001;70(4):569-79.
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

if (1 - abs(inDir(3))) < 1e-12
    inPhi = atan2(inDir(2), inDir(1));
    outDir = [sin(theta) * cos(phi+inphi); ...
              sin(theta) * sin(phi+inphi); ...
              cos(theta) * ((inDir(3) >= 0) - 0.5) * 2];
else 
    temp = sqrt(1 - inDir(3)^2);
    outDir = [sin(theta) * (inDir(1) * inDir(3) * cos(phi) - inDir(2) * sin(phi)) / temp + inDir(1) * cos(theta); ...
              sin(theta) * (inDir(2) * inDir(3) * cos(phi) + inDir(1) * sin(phi)) / temp + inDir(2) * cos(theta); ...
              -sin(theta) * cos(phi) * temp + inDir(3) * cos(theta)];
end

end
