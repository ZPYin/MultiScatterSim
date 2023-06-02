function [scaAng, rotAng1, rotAng2] = scaRotAng(inDir, outDir)
% SCAROTANG calculate scattering angle and rotation angle
%
% USAGE:
%    [scaAng, rotAng1, rotAng2] = scaRotAng(inDir, outDir)
%
% INPUTS:
%    params
%
% OUTPUTS:
%    output
%
% EXAMPLE:
%
% HISTORY:
%    2023-05-31: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

inPhi = atan(inDir(1) / inDir(2));
outPhi = atan(outDir(1) / outDir(2));
inTheta = acos(inDir(3));
outTheta = acos(outDir(3));

 if (sum(inDir .* outDir) >= 1)
    scaAng = 0.0;
 elseif (sum(inDir .* outDir) <= -1)
    scaAng = pi;
 else
    scaAng = acos(sum(inDir .* outDir));
 end

 if (abs(sin(scaAng)) < 1e-12)
    rotAng1 = 0;
    rotAng2 = 0;
    return;
 end

 if (1 - abs(inDir(3))) <= 1e-12
    rotAng1 = acos(inDir(3) * cos(outPhi - inPhi));
    rotAng2 = acos(inDir(3));
    return;
 end

 if (1 - abs(outDir(3))) <= 1e-12
    rotAng1 = acos(outDir(3) * cos(outPhi - inPhi));
    rotAng2 = acos(outDir(3));
    return;
 end

 rotAng1Cos = (outDir(3) - inDir(3) * cos(scaAng)) / (sin(inTheta) * sin(scaAng));
 rotAng2Cos = (inDir(3) - outDir(3) * cos(scaAng)) / (sin(outTheta) * sin(scaAng));

 if rotAng1Cos >= 1
    rotAng1Cos = 1;
 elseif rotAng1Cos <= -1
    rotAng1Cos = -1;
 end
 if rotAng2Cos >= 1
    rotAng2Cos = 1;
 elseif rotAng2Cos <= -1
    rotAng2Cos = -1;
 end

 rotAng1 = acos(rotAng1Cos);
 rotAng2 = acos(rotAng2Cos);

end