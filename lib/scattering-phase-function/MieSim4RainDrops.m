function [S1, S2] = MieSim4RainDrops(r, theta, inWL, varargin)
% MIESIM4RAINDROPS Mie Scattering simulation for certain shperes.
%
% 中文说明：
% 该函数用于计算给定半径数组 r、散射角数组 theta 和入射波长数组 inWL 下的米氏散射结果。
% 返回的 S1、S2 为米氏散射的两条复散射幅度分量（通常用于构造穆勒矩阵/相位函数）。
% 同时在内部计算消光效率 Qext 与散射效率 Qsca（不作为函数返回，但可选择保存到文件）。
%
% 使用示例：
%   r = linspace(0.5, 10, 50);         % 半径，单位：um
%   theta = linspace(0, 180, 361);     % 散射角，单位：度
%   inWL = [532];                      % 波长，单位：nm
%   [S1, S2] = MieSim4RainDrops(r, theta, inWL, 'saveFile', 'mie_outputs.mat');
%
% USAGE:
%    [S1, S2] = MieSim4RainDrops(r, theta, inWL)
%
% INPUTS:
%    r: array
%        radius of rain drops. (um)
%        中文：水滴半径数组，单位为微米 μm。
%    theta: array
%        scattering angles. (degree)
%        中文：散射角数组，单位为度（0~180°）。
%    inWL: array
%        incident wavelengths. (nm)
%        中文：入射光波长数组，单位为纳米 nm。
%
% KEYWORDS:
%    saveFile
%    debug
%    中文：
%      - saveFile：若提供非空字符串路径，则把计算得到的 S1/S2/Qext/Qsca 以及输入参数保存为 MAT 文件；
%      - debug：逻辑值，true 时在计算过程中打印当前进度信息（会显著降低计算速度）。
%
% OUTPUTS:
%    S1: 3-D matrix (wavelength x radius x theta)
%        S1 component from Mie Scattering.
%        中文：S1 分量，维度顺序为 [波长 x 半径 x 角度]。
%    S2: 3-D matrix (wavelength x radius x theta)
%        S2 component from Mie Scattering.
%        中文：S2 分量，维度顺序为 [波长 x 半径 x 角度]。
%    Qext: 3-D matrix (wavelength x radius x theta)
%        extinction efficiency.
%        中文：消光效率，维度顺序与 S1/S2 相同（函数内部计算；默认不返回）。
%    Qsca: 3-D matrix (wavelength x radius x theta)
%        scattering efficiency.
%        中文：散射效率，维度顺序与 S1/S2 相同（函数内部计算；默认不返回）。
%
% HISTORY:
%    2023-10-19: first edition by Zhenping
% .. Authors: - zp.yin@whu.edu.cn

% 参数解析器：支持必选参数 r/theta/inWL，以及可选参数 saveFile/debug
p = inputParser;
p.KeepUnmatched = true;           % 允许保留未匹配的可选参数（扩展兼容）

addRequired(p, 'r', @isnumeric);         % 半径数组（um）
addRequired(p, 'theta', @isnumeric);     % 散射角数组（度）
addRequired(p, 'inWL', @isnumeric);      % 入射波长数组（nm）
addParameter(p, 'saveFile', '', @ischar);% 保存文件路径（字符串），默认不保存
addParameter(p, 'debug', false, @islogical); % 调试开关，默认关闭

parse(p, r, theta, inWL, varargin{:});   % 解析输入参数

%% Parameter Definition
% r = [0.062,0.187,0.312,0.437,0.562,0.687,0.812,0.937,1.062,1.187,1.375,1.625,1.875,2.125,2.375,2.750,3.250,3.750,4.250,4.750,5.500] / 2 * 1e3;   % radius of water droplet. [um]
% r = (4000:0.1:4001) / 2;   % radius of water droplet. [um]
% theta = 180 - (0:0.001:0.1);
% inWL = 355;

%% Mie Scattering 主计算
% 预分配 3D 数组以提升速度，维度顺序为 [波长 x 半径 x 角度]
S1   = NaN(length(inWL), length(r), length(theta));
S2   = NaN(length(inWL), length(r), length(theta));
Qext = NaN(length(inWL), length(r), length(theta));
Qsca = NaN(length(inWL), length(r), length(theta));

for iWL = 1:length(inWL)
    % 计算该波长下的水的复折射率 m = n + i*k
    [m1, m2] = refIdxWater(inWL(iWL));

    for iR = 1:length(r)
        % 尺度参数 x = 2π r / λ，其中 r 为 μm、λ 为 nm，因此乘 1e3 做单位换算（μm/nm = 1e3）
        x = 2 * pi * r(iR) / inWL(iWL) * 1e3;

        for iTheta = 1:length(theta)
            % 可选的调试输出：当前波长/半径/角度（开启 debug 会极大降低速度）
            if (p.Results.debug)
                fprintf('Current sim: wavelength=%f nm; r=%f um; theta=%f\\deg\n', inWL(iWL), r(iR), theta(iTheta));
            end

            % 调用米氏散射计算核心：返回复散射幅度 S1/S2 以及效率系数 Qext/Qsca
            [S1(iWL, iR, iTheta), S2(iWL, iR, iTheta), Qext(iWL, iR, iTheta), Qsca(iWL, iR, iTheta), ] = MieCalc(m1 + 1i * m2, x, theta(iTheta));
        end
    end
end

% 以下为基于 S1/S2 构造相位矩阵（穆勒矩阵部分）的示例推导，默认注释不启用：
% P11 = 1 / 2 * (abs(S1).^2 + abs(S2).^2);                 % 总强度分量
% factor = trapz(theta / 180 * pi, P11, 2);                 % 按角度积分做归一化因子
% P12 = 1 / 2 * (abs(S2).^2 - abs(S1).^2);                 % 线偏振差分分量
% P33 = 1 / 2 * (conj(S2) .* S1 + S2 .* conj(S1));         % 交叉项（与圆偏振相关）
% P11 = P11 ./ repmat(factor, 1, length(theta));
% P12 = P12 ./ repmat(factor, 1, length(theta));
% P33 = P33 ./ repmat(factor, 1, length(theta));

% L1 = (P11 + P12);
% L2 = (P11 - P12);
% L3 = P33;                                           % 用于推导某些简化相位函数形式

% u = cos(repmat(theta, length(r), 1) / 180 * pi);
% P_p = 1 / 8 * (3 * L1 + 3 * L2 .* u.^2 + 2 * L3 .* u);   % 线偏振相函数示例
% P_c = 1 / 8 * (L1 + L2 .* u.^2 - 2 * L3 .* u);            % 圆偏振相函数示例

% 若需要保存相函数，也可以按需扩展保存内容：
% save(fullfile(MS_ENVS.RootPath, 'data', saveFile), 'theta', 'r', 'inWL', 'P_p', 'P_c', 'P11', 'P12', 'P33');
% 可选：保存输出到 MAT 文件（仅当提供了 'saveFile' 字符串路径时）
if ~ isempty(p.Results.saveFile)
    save(p.Results.saveFile, 'theta', 'r', 'inWL', 'S1', 'S2', 'Qext', 'Qsca');
end

end