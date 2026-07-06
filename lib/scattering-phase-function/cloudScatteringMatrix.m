function [cloudScatMatrix] = cloudScatteringMatrix(lambda, wdScatFile, varargin)
% CLOUDSCATTERINGMATRIX 计算云层的体平均散射矩阵及光学系数
%
% 用途：
%   根据水滴散射查找表（Mie 计算输出），结合云滴粒径分布（目前支持修正伽马分布），
%   计算每一云层的体平均偏振散射矩阵元素（S11, S12, S33, S34）以及消光系数 muExt 与散射系数 muSca。
%
% 调用方法：
%   cloudScatMatrix = cloudScatteringMatrix(lambda, wdScatFile, 'cloud_distribution', 'modified-gamma', ...)
%
% 输入参数：
%   lambda: 数值型，入射波长（单位：nm）
%   wdScatFile: 字符串，水滴散射查找表的绝对路径（由 ./lib/scattering-phase-function/water_droplet_scattering_gen.m 生成或等价格式）
%
% 关键字参数（Name-Value）：
%   cloud_distribution: 云滴粒径分布类型。当前支持 'modified-gamma'（修正伽马分布）。
%   kappa: 修正伽马分布形状参数（可为标量或与层数相同的向量）。
%   sigma_e: 云层的消光系数（m^-1），可为标量或与层数相同的向量。
%   R_eff: 云层的有效半径（um），可为标量或与层数相同的向量。
%
% 输出：
%   cloudScatMatrix: 结构体，包含以下字段：
%       muExt, muSca: 每层的消光与散射系数（单位：m^-1）
%       S11, S12, S33, S34: 体平均的偏振散射矩阵元素（每行为一层，列为散射角度）
%       theta: 散射角数组（单位：deg）
%
% 数据库要求与维度说明（典型）：
%   - mieScaInfo.theta: [1 x nTheta] 或 [nTheta x 1]，单位为度（deg）
%   - mieScaInfo.r: [1 x nR] 或 [nR x 1]，半径（um）
%   - mieScaInfo.S1_Real/S1_Imag/S2_Real/S2_Imag: 若为单波长，通常为 [nR x nTheta]；
%     若包含多个波长，则可能为 [nWL x nR x nTheta]。本实现自动判断并处理这两种情况。
%   - mieScaInfo.Qext/Qsca: 可选字段；若不存在，则用输入 sigma_e 作为 muExt，并假设单次散射反照率≈1（muSca≈muExt）。
%
% 计算流程概述：
%   1) 读取散射查找表，检查波长是否匹配（若 inWL 缺失则跳过检查）。
%   2) 构造 S1/S2 复数散射振幅，计算 Mueller 矩阵中用到的 S11, S12, S33, S34。
%   3) 对角度进行数值积分归一化（使 ∫ S11 dθ ≈ 1，θ单位为度），保证后续体平均不受角度采样密度影响。
%   4) 针对每一云层，根据形参 kappa、sigma_e、R_eff 计算粒径分布 n(r)。
%   5) 按半径维进行加权积分：得到 muExt/muSca（若可用）以及体平均的 S11/S12/S33/S34。
%
% 历史：
%   2023-06-02: 初版 by Zhenping
%   2025-11-04: 增强对不同数据库格式的兼容性与数值稳定性
% .. Authors: - zp.yin@whu.edu.cn

% 添加必要的路径
projectDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectDir, 'lib', 'cloud-microphysics'));  % 引入 smmgDist 等云微物理函数

p = inputParser;            % 解析 Name-Value 输入
p.KeepUnmatched = true;     % 保留未匹配参数，便于扩展

addRequired(p, 'lambda', @isnumeric);         % 波长（nm）
addRequired(p, 'wdScatFile', @ischar);        % 水滴散射查找表路径（绝对路径）
addParameter(p, 'cloud_distribution', 'modified-gamma', @ischar);
addParameter(p, 'kappa', 7, @isnumeric);      % 修正伽马分布形状参数（原版）
addParameter(p, 'b', 1/9, @isnumeric);        % 论文版修正伽马分布形状参数（smmgDist_paper）
addParameter(p, 'sigma_e', 2e-2, @isnumeric); % 云层消光系数（m^-1）—— 没有 N_d 时使用
addParameter(p, 'R_eff', 5, @isnumeric);      % 有效半径（um）
addParameter(p, 'N_d', [], @isnumeric);        % 数浓度 [m^-3]，优先级高于 sigma_e（多波长同云模式）

parse(p, lambda, wdScatFile, varargin{:});

cloudScatMatrix = struct();                   % 预分配输出结构体

mieScaInfo = load(wdScatFile);                % 读入 Mie 散射查找表（MAT 文件）

% 小工具：支持标量参数自动扩展到各云层。
pickLayerValue = @(x, iLayer) x(min(iLayer, numel(x)));

% 根据云微物理分布类型，选择不同的体平均策略（当前仅实现 modified-gamma）
switch lower(p.Results.cloud_distribution)
case 'modified-gamma'
    % 云层数：由最长的层参数决定，标量参数自动扩展。
    nClLayer = max([numel(p.Results.sigma_e), numel(p.Results.R_eff), ...
        numel(p.Results.kappa), max(1, numel(p.Results.N_d))]);

    % 角度与半径矢量，统一形状
    theta = mieScaInfo.theta;
    theta = theta(:)';        % 行向量 [1 x nTheta]，单位：deg
    nTheta = length(theta);

    r = mieScaInfo.r;
    r = r(:);                 % 列向量 [nR x 1]，单位：um
    nR = length(r);

    % 波长兼容性检查（若数据库包含 inWL 则进行检查，否则跳过）
    if isfield(mieScaInfo, 'inWL')
        try
            inWL = mieScaInfo.inWL(:);
        catch
            inWL = mieScaInfo.inWL; % 兼容非向量存储
        end
        if ~any(abs(inWL - lambda) < 1)
            error('Scattering properties at given wavelength (%f nm) are not available in the database.', lambda);
        end
        wlIdx = find(abs(inWL - lambda) < 1, 1, 'first');
    else
        wlIdx = 1; % 兼容仅含单波长数据的文件
    end

    % 预分配各层输出变量：
    muExt = NaN(nClLayer, 1);
    muSca = NaN(nClLayer, 1);
    S11_PD = NaN(nClLayer, nTheta);
    S12_PD = NaN(nClLayer, nTheta);
    S33_PD = NaN(nClLayer, nTheta);
    S34_PD = NaN(nClLayer, nTheta);

    % 构造复数散射振幅 S1/S2：自动兼容 2D/3D 数据格式
    S1r = mieScaInfo.S1_Real; S1i = mieScaInfo.S1_Imag;
    S2r = mieScaInfo.S2_Real; S2i = mieScaInfo.S2_Imag;
    if ndims(S1r) == 3
        S1 = squeeze(S1r(wlIdx, :, :)) + 1i * squeeze(S1i(wlIdx, :, :));
        S2 = squeeze(S2r(wlIdx, :, :)) + 1i * squeeze(S2i(wlIdx, :, :));
    else
        S1 = S1r + 1i * S1i;
        S2 = S2r + 1i * S2i;
    end
    % 保证 S1/S2 维度为 [nR x nTheta]
    [s1_m, s1_n] = size(S1);
    if s1_m == nTheta && s1_n == nR
        S1 = S1'; S2 = S2';
    elseif s1_m ~= nR || s1_n ~= nTheta
        error('S1/S2 dimensions do not match r/theta: expected [%d x %d], got [%d x %d].', nR, nTheta, s1_m, s1_n);
    end

    % 由 S1/S2 计算 Mueller 矩阵中涉及的 4 个元素（标准公式）：
    S11 = 0.5 * (abs(S1).^2 + abs(S2).^2);
    S12 = 0.5 * (abs(S2).^2 - abs(S1).^2);
    S33 = real(0.5 * (conj(S2) .* S1 + S2 .* conj(S1))); % 取实部以抑制数值噪声
    S34 = real(imag(conj(S1) .* (S2)));                    % 取实部以抑制数值噪声

    % 角度归一化（单位：deg）：使 ∫ S11(θ) sin(θ) dθ ≈ 1
    % 正确的三维空间相函数归一化应该包含立体角权重 sin(θ)
    thetaDeg = theta;
    thetaRad = thetaDeg * pi / 180;
    
    % 标准公式：1/2 * 积分( P(θ)*sinθ dθ ) = 1 => 积分值应为 2
    % 注意：积分操作必须在弧度(thetaRad)下进行！
    integral_val = trapz(thetaRad, S11 .* repmat(sin(thetaRad), [size(S11, 1), 1]), 2);    % [nR x 1]
    norm_ang = integral_val / 2;
    
    norm_ang(norm_ang == 0) = eps;
    S11 = S11 ./ norm_ang;
    S12 = S12 ./ norm_ang;
    S33 = S33 ./ norm_ang;
    S34 = S34 ./ norm_ang;

    % 可选的 Qext/Qsca（若不存在则回退到输入 sigma_e）
    hasQext = isfield(mieScaInfo, 'Qext');
    hasQsca = isfield(mieScaInfo, 'Qsca');

    for iLayer = 1:nClLayer
        % —— 针对第 iLayer 层进行体平均 ——
        % 计算云滴粒径分布 n(r)
        kappa_layer = pickLayerValue(p.Results.kappa, iLayer);
        reff_layer = pickLayerValue(p.Results.R_eff, iLayer);
        % 构造 smmgDist 参数：N_d 优先（多波长同云），否则用 sigma_e
        smmgArgs = {'kappa', kappa_layer, 'R_eff', reff_layer};
        if ~isempty(p.Results.N_d)
            smmgArgs = [smmgArgs, {'N_d', pickLayerValue(p.Results.N_d, iLayer)}];
        else
            smmgArgs = [smmgArgs, {'sigma_e', pickLayerValue(p.Results.sigma_e, iLayer)}];
        end
        n_r = smmgDist(r, smmgArgs{:});

        % 计算体平均的消光/散射系数（稳健地处理 Qext/Qsca 的维度差异）
        if hasQext && hasQsca
            Qext_raw = mieScaInfo.Qext; Qsca_raw = mieScaInfo.Qsca;
            Qext_v = []; Qsca_v = [];

            % —— 提取 Qext 与 r 对应的一维向量 ——
            if isvector(Qext_raw)
                tmp = Qext_raw(:);
                if length(tmp) == nR
                    Qext_v = tmp;
                end
            else
                sz = size(Qext_raw);
                if numel(sz) == 2
                    m = sz(1); n = sz(2);
                    if n == nR
                        row = min(wlIdx, m);
                        Qext_v = Qext_raw(row, :).';
                    elseif m == nR
                        col = min(wlIdx, n);
                        Qext_v = Qext_raw(:, col);
                    end
                elseif numel(sz) == 3
                    if sz(2) == nR
                        Qext_v = squeeze(Qext_raw(min(wlIdx, sz(1)), :, 1));
                    elseif sz(1) == nR
                        Qext_v = squeeze(Qext_raw(:, min(wlIdx, sz(2)), 1));
                    end
                end
            end

            % —— 提取 Qsca 与 r 对应的一维向量 ——
            if isvector(Qsca_raw)
                tmp = Qsca_raw(:);
                if length(tmp) == nR
                    Qsca_v = tmp;
                end
            else
                sz = size(Qsca_raw);
                if numel(sz) == 2
                    m = sz(1); n = sz(2);
                    if n == nR
                        row = min(wlIdx, m);
                        Qsca_v = Qsca_raw(row, :).';
                    elseif m == nR
                        col = min(wlIdx, n);
                        Qsca_v = Qsca_raw(:, col);
                    end
                elseif numel(sz) == 3
                    if sz(2) == nR
                        Qsca_v = squeeze(Qsca_raw(min(wlIdx, sz(1)), :, 1));
                    elseif sz(1) == nR
                        Qsca_v = squeeze(Qsca_raw(:, min(wlIdx, sz(2)), 1));
                    end
                end
            end

            % 若无法稳定匹配 nR，则回退到 sigma_e
            if isempty(Qext_v) || isempty(Qsca_v) || length(Qext_v) ~= nR || length(Qsca_v) ~= nR
                muExt(iLayer) = pickLayerValue(p.Results.sigma_e, iLayer);
                muSca(iLayer) = pickLayerValue(p.Results.sigma_e, iLayer);
            else
                muExt(iLayer) = pi * (1e-6)^2 * trapz(r, n_r .* (r.^2) .* Qext_v);
                muSca(iLayer) = pi * (1e-6)^2 * trapz(r, n_r .* (r.^2) .* Qsca_v);
            end
        else
            % 若数据库不含 Qext/Qsca，则使用输入的 sigma_e 作为 muExt，假设 muSca≈muExt
            muExt(iLayer) = pickLayerValue(p.Results.sigma_e, iLayer);
            muSca(iLayer) = pickLayerValue(p.Results.sigma_e, iLayer);
        end

        % 计算体平均的偏振散射矩阵元素：使用散射截面加权
        % weight = n(r) × r² × Qsca(r)，大粒子散射主导，抑制小粒子Mie振荡
        if hasQext && hasQsca
            weight = n_r .* (r.^2) .* max(Qsca_v, 1e-30);
        else
            weight = n_r .* (r.^2);   % 退化为几何截面加权
        end
        denom = trapz(r, weight);
        if denom == 0, denom = eps; end
        w = weight / denom;                     % 归一化权重 [nR x 1]
        W = repmat(w, 1, nTheta);
        S11_PD(iLayer, :) = trapz(r, S11 .* W, 1);
        S12_PD(iLayer, :) = trapz(r, S12 .* W, 1);
        S33_PD(iLayer, :) = trapz(r, S33 .* W, 1);
        S34_PD(iLayer, :) = trapz(r, S34 .* W, 1);
    end

    % 汇总输出：
    cloudScatMatrix.muExt = muExt;
    cloudScatMatrix.muSca = muSca;
    cloudScatMatrix.S11 = S11_PD;
    cloudScatMatrix.S12 = S12_PD;
    cloudScatMatrix.S33 = S33_PD;
    cloudScatMatrix.S34 = S34_PD;
    cloudScatMatrix.theta = theta; % deg, 行向量

case 'modified-gamma-paper'
    % 使用论文版修正伽马分布（概率密度）进行体平均散射矩阵计算
    nClLayer = length(p.Results.sigma_e);

    theta = mieScaInfo.theta; theta = theta(:)'; nTheta = length(theta);
    r = mieScaInfo.r; r = r(:); nR = length(r);

    if isfield(mieScaInfo, 'inWL')
        try, inWL = mieScaInfo.inWL(:); catch, inWL = mieScaInfo.inWL; end
        if ~any(abs(inWL - lambda) < 1)
            error('Scattering properties at given wavelength (%f nm) are not available in the database.', lambda);
        end
        wlIdx = find(abs(inWL - lambda) < 1, 1, 'first');
    else
        wlIdx = 1;
    end

    muExt = NaN(nClLayer, 1);
    muSca = NaN(nClLayer, 1);
    S11_PD = NaN(nClLayer, nTheta);
    S12_PD = NaN(nClLayer, nTheta);
    S33_PD = NaN(nClLayer, nTheta);
    S34_PD = NaN(nClLayer, nTheta);

    S1r = mieScaInfo.S1_Real; S1i = mieScaInfo.S1_Imag;
    S2r = mieScaInfo.S2_Real; S2i = mieScaInfo.S2_Imag;
    if ndims(S1r) == 3
        S1 = squeeze(S1r(wlIdx, :, :)) + 1i * squeeze(S1i(wlIdx, :, :));
        S2 = squeeze(S2r(wlIdx, :, :)) + 1i * squeeze(S2i(wlIdx, :, :));
    else
        S1 = S1r + 1i * S1i; S2 = S2r + 1i * S2i;
    end
    [s1_m, s1_n] = size(S1);
    if s1_m == nTheta && s1_n == nR
        S1 = S1'; S2 = S2';
    elseif s1_m ~= nR || s1_n ~= nTheta
        error('S1/S2 dimensions do not match r/theta: expected [%d x %d], got [%d x %d].', nR, nTheta, s1_m, s1_n);
    end

    S11 = 0.5 * (abs(S1).^2 + abs(S2).^2);
    S12 = 0.5 * (abs(S2).^2 - abs(S1).^2);
    S33 = real(0.5 * (conj(S2) .* S1 + S2 .* conj(S1)));
    S34 = real(-imag(S1 .* conj(S2)));

    thetaDeg = theta;
    thetaRad = thetaDeg * pi / 180;
    integral_val = trapz(thetaRad, S11 .* repmat(sin(thetaRad), [size(S11, 1), 1]), 2);
    norm_ang = integral_val / 2;
    norm_ang(norm_ang == 0) = eps;
    S11 = S11 ./ norm_ang; S12 = S12 ./ norm_ang; S33 = S33 ./ norm_ang; S34 = S34 ./ norm_ang;

    % 提取 Qsca（若数据库包含），用于散射截面加权
    hasQsca_paper = isfield(mieScaInfo, 'Qsca');
    if hasQsca_paper
        Qsca_raw = mieScaInfo.Qsca;
        Qsca_paper = [];
        if isvector(Qsca_raw)
            tmp = Qsca_raw(:);
            if length(tmp) == nR, Qsca_paper = tmp; end
        elseif ndims(Qsca_raw) == 2
            sz = size(Qsca_raw);
            if sz(2) == nR, Qsca_paper = Qsca_raw(min(wlIdx, sz(1)), :).'; end
            if sz(1) == nR, Qsca_paper = Qsca_raw(:, min(wlIdx, sz(2))); end
        elseif ndims(Qsca_raw) == 3
            Qsca_paper = squeeze(Qsca_raw(min(wlIdx, size(Qsca_raw, 1)), :, 1));
            Qsca_paper = Qsca_paper(:);
        end
        if isempty(Qsca_paper)
            hasQsca_paper = false;
        end
    end

    for iLayer = 1:nClLayer
        % 取论文版形参 b
        if length(p.Results.b) == 1
            b_layer = p.Results.b;
        else
            b_layer = p.Results.b(iLayer);
        end
        % 计算概率密度 f(r)，在给定 r 网格上做归一化（论文版参数：a=R_eff）
        f_r = smmgDist_paper(r, 'b', b_layer, 'a', p.Results.R_eff(iLayer), 'normalize', 'grid');

        % 散射截面权重：w(r) = f(r) × r² × Qsca(r)，大粒子散射主导
        if hasQsca_paper
            weight = f_r .* (r.^2) .* max(Qsca_paper, 1e-30);
        else
            weight = f_r .* (r.^2);   % 退化：几何截面加权
        end
        denom = trapz(r, weight);
        if denom == 0, denom = eps; end
        w = weight / denom; W = repmat(w, 1, nTheta);

        S11_PD(iLayer, :) = trapz(r, S11 .* W, 1);
        S12_PD(iLayer, :) = trapz(r, S12 .* W, 1);
        S33_PD(iLayer, :) = trapz(r, S33 .* W, 1);
        S34_PD(iLayer, :) = trapz(r, S34 .* W, 1);

        % 消光/散射系数: N_d 模式由 Mie Qext 积分得到, sigma_e 模式直接使用输入
        if ~isempty(p.Results.N_d) && isfield(mieScaInfo, 'Qext')
            Qe = mieScaInfo.Qext;
            if ndims(Qe) == 3
                Qe = squeeze(Qe(wlIdx, :, :));
            end
            Qe = Qe(:);
            if length(Qe) == nR
                muExt(iLayer) = pi * (1e-6)^2 * trapz(r, f_r .* (r.^2) .* Qe);
                muSca(iLayer) = muExt(iLayer);   % 近似
            else
                muExt(iLayer) = p.Results.sigma_e(iLayer);
                muSca(iLayer) = p.Results.sigma_e(iLayer);
            end
        else
            muExt(iLayer) = p.Results.sigma_e(iLayer);
            muSca(iLayer) = p.Results.sigma_e(iLayer);
        end
    end

    cloudScatMatrix.muExt = muExt;
    cloudScatMatrix.muSca = muSca;
    cloudScatMatrix.S11 = S11_PD;
    cloudScatMatrix.S12 = S12_PD;
    cloudScatMatrix.S33 = S33_PD;
    cloudScatMatrix.S34 = S34_PD;
    cloudScatMatrix.theta = theta;

case 'log-normal'
otherwise
    error('Unknown cloud microphysical distribution: %s', p.Results.cloud_distribution);
end

end
