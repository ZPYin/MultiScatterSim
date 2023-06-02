function cloud_scattering_gen(r, theta, inWL, varargin)
% CLOUD_SCATTERING_GEN generate cloud droplet scattering function.
%
% USAGE:
%    cloud_scattering_gen(r, theta, inWL)
%
% INPUTS:
%    r: numeric
%        radius. [um]
%    theta: numeric
%        scattering angle. [degree]
%    inWL: numeric
%        incident wavelength. [nm]
%
% KEYWORDS:
%    fileType: char
%        'mat' or 'netcdf'
%    saveFile: char
%        absolute path for saving scattering phase function.
%
% HISTORY:
%    2023-05-24: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'r', @isnumeric);
addRequired(p, 'theta', @isnumeric);
addRequired(p, 'inWL', @isnumeric);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'fileType', 'netcdf', @ischar);
addParameter(p, 'saveFile', '', @ischar);

parse(p, r, theta, inWL, varargin{:});

%% Mie Calculation
S1 = NaN(length(inWL), length(r), length(theta));
S2 = NaN(length(inWL), length(r), length(theta));
Qext = NaN(length(inWL), length(r));
Qsca = NaN(length(inWL), length(r));
for iWL = 1:length(inWL)
    [m1, m2] = refIdxWater(inWL(iWL));

    for iR = 1:length(r)
        x = 2 * pi * r(iR) / inWL(iWL) * 1e3;

        for iTheta = 1:length(theta)
            [S1(iWL, iR, iTheta), S2(iWL, iR, iTheta), Qext(iWL, iR), Qsca(iWL, iR)] = MieCalc(m1 + 1i * m2, x, theta(iTheta));
        end
    end
end

S1_Real = real(S1);
S1_Imag = imag(S1);
S2_Real = real(S2);
S2_Imag = imag(S2);

if ~ isempty(p.Results.saveFile)
    if strcmpi(p.Results.fileType, 'mat')

        save(p.Results.saveFile, 'r', 'theta', 'inWL', 'S1_Real', 'S1_Imag', 'S2_Real', 'S2_Imag', 'Qext', 'Qsca');

    else if strcmpi(p.Results.fileType, 'netcdf')
        mode = netcdf.getConstant('NETCDF4');
        mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
        mode = bitor(mode, netcdf.getConstant('CLOBBER'));
        ncID = netcdf.create(p.Results.saveFile, mode);

        % define dimensions
        dimID_r = netcdf.defDim(ncID, 'radius', length(r));
        dimID_t = netcdf.defDim(ncID, 'scattering_angle', length(theta));
        dimID_w = netcdf.defDim(ncID, 'wavelength', length(inWL));

        % define variables
        varID_r = netcdf.defVar(ncID, 'radius', 'NC_DOUBLE', [dimID_r]);
        varID_t = netcdf.defVar(ncID, 'scattering_angle', 'NC_DOUBLE', [dimID_t]);
        varID_WL = netcdf.defVar(ncID, 'wavelength', 'NC_DOUBLE', [dimID_w]);
        varID_S1_Real = netcdf.defVar(ncID, 'S1_Real', 'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S1_Imag = netcdf.defVar(ncID, 'S1_Imaginary', 'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S2_Real = netcdf.defVar(ncID, 'S2_Real', 'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S2_Imag = netcdf.defVar(ncID, 'S2_Imaginary', 'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_Qext = netcdf.defVar(ncID, 'Qext', 'NC_DOUBLE', [dimID_w, dimID_r]);
        varID_Qsca = netcdf.defVar(ncID, 'Qsca', 'NC_DOUBLE', [dimID_w, dimID_r]);

        % define data compression
        netcdf.defVarDeflate(ncID, varID_S1_Real, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S1_Imag, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S2_Real, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S2_Imag, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_Qext, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_Qsca, true, true, 5);

        % leave define mode
        netcdf.endDef(ncID);

        % write data to NC file
        netcdf.putVar(ncID, varID_r, r);
        netcdf.putVar(ncID, varID_t, theta);
        netcdf.putVar(ncID, varID_c, inWL);
        netcdf.putVar(ncID, varID_S1_Real, S1_Real);
        netcdf.putVar(ncID, varID_S1_Imag, S1_Imag);
        netcdf.putVar(ncID, varID_S2_Real, S1_Real);
        netcdf.putVar(ncID, varID_S2_Imag, S1_Imag);
        netcdf.putVar(ncID, varID_Qext, Qext);
        netcdf.putVar(ncID, varID_Qsca, Qsca);

        % re-enter define mode
        netcdf.reDef(ncID);

        netcdf.putAtt(ncID, varID_r, 'unit', 'micron');
        netcdf.putAtt(ncID, varID_r, 'name', 'r');
        netcdf.putAtt(ncID, varID_r, 'definition', 'radius of particle');

        netcdf.putAtt(ncID, varID_t, 'unit', 'degree');
        netcdf.putAtt(ncID, varID_t, 'name', '\theta');
        netcdf.putAtt(ncID, varID_t, 'definition', 'scattering angle');

        netcdf.putAtt(ncID, varID_c, 'unit', 'nm');
        netcdf.putAtt(ncID, varID_c, 'name', 'lambda');
        netcdf.putAtt(ncID, varID_c, 'definition', 'incident wavelength');

        netcdf.putAtt(ncID, varID_S1_Real, 'unit', '');
        netcdf.putAtt(ncID, varID_S1_Real, 'name', 'S1_Real');
        netcdf.putAtt(ncID, varID_S1_Real, 'definition', 'Real part of S1');

        netcdf.putAtt(ncID, varID_S1_Imag, 'unit', '');
        netcdf.putAtt(ncID, varID_S1_Imag, 'name', 'S1_Imag');
        netcdf.putAtt(ncID, varID_S1_Imag, 'definition', 'Imaginary part of S1');

        netcdf.putAtt(ncID, varID_S2_Real, 'unit', '');
        netcdf.putAtt(ncID, varID_S2_Real, 'name', 'S2_Real');
        netcdf.putAtt(ncID, varID_S2_Real, 'definition', 'Real part of S2');

        netcdf.putAtt(ncID, varID_S2_Imag, 'unit', '');
        netcdf.putAtt(ncID, varID_S2_Imag, 'name', 'S2_Imag');
        netcdf.putAtt(ncID, varID_S2_Imag, 'definition', 'Imaginary part of S2');

        netcdf.putAtt(ncID, varID_S2_Imag, 'unit', '');
        netcdf.putAtt(ncID, varID_S2_Imag, 'name', 'Qext');
        netcdf.putAtt(ncID, varID_S2_Imag, 'definition', 'Extinction efficiency');

        netcdf.putAtt(ncID, varID_S2_Imag, 'unit', '');
        netcdf.putAtt(ncID, varID_S2_Imag, 'name', 'Qsca');
        netcdf.putAtt(ncID, varID_S2_Imag, 'definition', 'Scattering efficiency');

        % write global attributes
        varID_global = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncID, varID_Z, 'License', 'CF-1.0');

        % close NC file
        netcdf.close(ncID);
    else
        error('Unknown file type: %s', p.Results.fileType);
    end
end

end