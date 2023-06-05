function [SVOut] = rotSVPhi(SVIn, phi)
% ROTSVPHI Rotate the Stokes vection in degree of phi.
%
% USAGE:
%    [SVOut] = rotSVPhi(SVIn, phi)
%
% INPUTS:
%    SVIn: vector
%        input stokes vector.
%    phi: numeric
%        rotation angle. (radians)
%
% OUTPUTS:
%    SVOut: vector
%        output stokes vector.
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

SVOut = [SVIn(1); ...
         SVIn(2) * cos(2 * phi) + SVIn(3) * sin(2 * phi); ...
         -SVIn(2) * sin(2 * phi) + SVIn(3) * cos(2 * phi); ...
         SVIn(4)];

end