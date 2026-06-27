function water_droplet_scattering_gen(r, theta, inWL, varargin)
% WATER_DROPLET_SCATTERING_GEN 使用米氏散射（Mie theory）生成水滴散射数据。
%
% 作用
%   在给定的半径 r[微米]、散射角度 theta[度]、入射波长 inWL[纳米] 的网格上，
%   计算米氏散射振幅函数 S1、S2（复数），以及消光效率 Qext 和散射效率 Qsca。
%   可选择将结果保存为 MAT 或 NetCDF 文件。
%
% 用法示例
%   water_droplet_scattering_gen(r, theta, inWL)
%   water_droplet_scattering_gen(r, theta, inWL, 'fileType', 'mat', 'saveFile', 'out.mat')
%   water_droplet_scattering_gen(r, theta, inWL, 'fileType', 'netcdf', 'saveFile', 'out.nc')
%
% 输入
%   r     : 数值向量，水滴半径，单位微米 [um]
%   theta : 数值向量，散射角度，单位度 [deg]
%   inWL  : 数值向量，入射波长，单位纳米 [nm]
%
% 名值参数（可选）
%   'debug'    : 逻辑值，默认 false。为 true 时打印计算进度信息。
%   'fileType' : 字符串，'mat' 或 'netcdf'，默认 'netcdf'，决定输出格式。
%   'saveFile' : 字符串，输出文件路径。不为空则保存；为空则只计算不保存。
%
% 输出（写入文件中的变量）
%   S1_Real, S1_Imag, S2_Real, S2_Imag : 三维数组 [nWL × nR × nTheta]
%   Qext, Qsca                         : 二维数组 [nWL × nR]
%   r, theta, inWL                     : 与结果一并保存的坐标网格
%
% 说明
%   - S1/S2 是复数，保存时将其实部、虚部分开，避免在文件中写复数类型。
%   - Qext/Qsca 与散射角无关；本实现在角度循环中重复赋值，便于代码结构，
%     也可以优化为每个 (WL, R) 只计算一次。
%   - 传入 MieCalc 的角度单位为“度”，与库函数保持一致。
%   - 水的复折射率由 refIdxWater(inWL) 提供。
%
% 历史
%   2023-05-24: Zhenping 首次编写
%   2025-10-24: 补充中文注释，并整理 NetCDF 保存段中的变量命名与属性写入（罗澄炜）
% .. Authors: - zp.yin@whu.edu.cn

%-----------------------------
% 参数解析
%-----------------------------
p = inputParser;
p.KeepUnmatched = true;

addRequired(p, 'r', @isnumeric);
addRequired(p, 'theta', @isnumeric);
addRequired(p, 'inWL', @isnumeric);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'fileType', 'netcdf', @ischar);
addParameter(p, 'saveFile', '', @ischar);

parse(p, r, theta, inWL, varargin{:});

%-----------------------------
% 预分配结果数组（维度顺序：波长 × 半径 × 角度）
%-----------------------------
S1   = NaN(length(inWL), length(r), length(theta));
S2   = NaN(length(inWL), length(r), length(theta));
Qext = NaN(length(inWL), length(r));
Qsca = NaN(length(inWL), length(r));

%-----------------------------
% 米氏散射计算主循环
%-----------------------------
for iWL = 1:length(inWL)
    % 获取当前波长下水的复折射率 m = m1 + i*m2
    [m1, m2] = refIdxWater(inWL(iWL));
    if p.Results.debug
        fprintf('WL=%gnm -> m=%g + %gi\n', inWL(iWL), m1, m2);
    end

    for iR = 1:length(r)
        % 尺寸参数 x = 2*pi*r / lambda
        % 注意单位：r 为微米，需要转换为纳米（×1e3）与 lambda[nm] 匹配
        x = 2 * pi * r(iR) / inWL(iWL) * 1e3;

        for iTheta = 1:length(theta)
            % 计算当前 (WL, R, theta) 的散射振幅与效率
            % Qext/Qsca 与角度无关，重复赋值不会影响结果
            [S1(iWL, iR, iTheta), S2(iWL, iR, iTheta), Qext(iWL, iR), Qsca(iWL, iR)] = ...
                MieCalc(m1 + 1i * m2, x, theta(iTheta));
        end
    end
end

%-----------------------------
% 拆分复数为实部与虚部，便于保存
%-----------------------------
S1_Real = real(S1);
S1_Imag = imag(S1);
S2_Real = real(S2);
S2_Imag = imag(S2);

%-----------------------------
% 根据需要保存结果
%-----------------------------
if ~isempty(p.Results.saveFile)
    if strcmpi(p.Results.fileType, 'mat')
        % 保存为 MAT 文件
        save(p.Results.saveFile, 'r', 'theta', 'inWL', 'S1_Real', 'S1_Imag', 'S2_Real', 'S2_Imag', 'Qext', 'Qsca');

    elseif strcmpi(p.Results.fileType, 'netcdf')
        % 保存为 NetCDF 文件（带压缩与基础属性）
        mode = netcdf.getConstant('NETCDF4');
        mode = bitor(mode, netcdf.getConstant('CLASSIC_MODEL'));
        mode = bitor(mode, netcdf.getConstant('CLOBBER'));
        ncID = netcdf.create(p.Results.saveFile, mode);

        % 定义维度
        dimID_r = netcdf.defDim(ncID, 'radius', length(r));
        dimID_t = netcdf.defDim(ncID, 'scattering_angle', length(theta));
        dimID_w = netcdf.defDim(ncID, 'wavelength', length(inWL));

        % 定义变量（注意维度顺序：[wavelength, radius, angle]）
        varID_r       = netcdf.defVar(ncID, 'radius',           'NC_DOUBLE', [dimID_r]);
        varID_t       = netcdf.defVar(ncID, 'scattering_angle', 'NC_DOUBLE', [dimID_t]);
        varID_WL      = netcdf.defVar(ncID, 'wavelength',       'NC_DOUBLE', [dimID_w]);
        varID_S1_Real = netcdf.defVar(ncID, 'S1_Real',          'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S1_Imag = netcdf.defVar(ncID, 'S1_Imaginary',     'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S2_Real = netcdf.defVar(ncID, 'S2_Real',          'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_S2_Imag = netcdf.defVar(ncID, 'S2_Imaginary',     'NC_DOUBLE', [dimID_w, dimID_r, dimID_t]);
        varID_Qext    = netcdf.defVar(ncID, 'Qext',             'NC_DOUBLE', [dimID_w, dimID_r]);
        varID_Qsca    = netcdf.defVar(ncID, 'Qsca',             'NC_DOUBLE', [dimID_w, dimID_r]);

        % 开启压缩（提高存储效率）
        netcdf.defVarDeflate(ncID, varID_S1_Real, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S1_Imag, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S2_Real, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_S2_Imag, true, true, 5);
        netcdf.defVarDeflate(ncID, varID_Qext,    true, true, 5);
        netcdf.defVarDeflate(ncID, varID_Qsca,    true, true, 5);

        % 写入数据（退出定义模式）
        netcdf.endDef(ncID);
        netcdf.putVar(ncID, varID_r,  r);
        netcdf.putVar(ncID, varID_t,  theta);
        netcdf.putVar(ncID, varID_WL, inWL);
        netcdf.putVar(ncID, varID_S1_Real, S1_Real);
        netcdf.putVar(ncID, varID_S1_Imag, S1_Imag);
        netcdf.putVar(ncID, varID_S2_Real, S2_Real);
        netcdf.putVar(ncID, varID_S2_Imag, S2_Imag);
        netcdf.putVar(ncID, varID_Qext,    Qext);
        netcdf.putVar(ncID, varID_Qsca,    Qsca);

        % 写属性（重新进入定义模式）
        netcdf.reDef(ncID);
        netcdf.putAtt(ncID, varID_r,  'unit', 'micron');
        netcdf.putAtt(ncID, varID_r,  'name', 'r');
        netcdf.putAtt(ncID, varID_r,  'definition', 'radius of particle');

        netcdf.putAtt(ncID, varID_t,  'unit', 'degree');
        netcdf.putAtt(ncID, varID_t,  'name', '\theta');
        netcdf.putAtt(ncID, varID_t,  'definition', 'scattering angle');

        netcdf.putAtt(ncID, varID_WL, 'unit', 'nm');
        netcdf.putAtt(ncID, varID_WL, 'name', 'lambda');
        netcdf.putAtt(ncID, varID_WL, 'definition', 'incident wavelength');

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

        netcdf.putAtt(ncID, varID_Qext, 'unit', '');
        netcdf.putAtt(ncID, varID_Qext, 'name', 'Qext');
        netcdf.putAtt(ncID, varID_Qext, 'definition', 'Extinction efficiency');

        netcdf.putAtt(ncID, varID_Qsca, 'unit', '');
        netcdf.putAtt(ncID, varID_Qsca, 'name', 'Qsca');
        netcdf.putAtt(ncID, varID_Qsca, 'definition', 'Scattering efficiency');

        % 全局属性与文件关闭
        varID_global = netcdf.getConstant('GLOBAL');
        netcdf.putAtt(ncID, varID_global, 'License', 'CF-1.0');
        netcdf.close(ncID);
    else
        error('Unknown file type: %s', p.Results.fileType);
    end
end

end