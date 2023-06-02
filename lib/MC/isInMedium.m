function [status] = isInMedium(pos, clBound)
% ISINMEDIUM determine whether the current photon is in the medium.
%
% USAGE:
%    [status] = isInMedium(pos, clBound)
%
% INPUTS:
%    pos: numeric
%        photon position [x, y, z]
%    clBound: numeric
%        cloud boundaries. [cloud-base, cloud top, cloudBoundX, cloudBoundY]
%
% OUTPUTS:
%    status: logical
%        true or false.
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

status = ((pos(1) - clBound(3)) * (pos(1) + clBound(3)) < 0) && ...
         ((pos(2) - clBound(4)) * (pos(1) + clBound(4)) < 0) && ...
         ((pos(3) - clBound(1)) * (pos(3) - clBound(2)) < 0);

end