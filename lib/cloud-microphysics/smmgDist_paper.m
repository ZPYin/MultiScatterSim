function [f_r, const] = smmgDist_paper(r, varargin)
% SMMGDIST_PAPER 修正伽马分布（对方论文公式）
%
% 功能：
%   根据论文中给出的修正伽马分布公式
%       f(r) = const * r^((1-3b)/b) * exp(- r / (a*b))
%   计算云滴半径 r 的概率密度函数（归一化常数 const 自动计算）。
%
%   其中：
%     - const 为归一化常数（使得 ∫_0^∞ f(r) dr = 1，连续归一化）
%     - a 等于有效半径 r_eff（单位：μm）
%     - b 为形状参数（论文给定 b = 1/9）
%
% 输入：
%   r : 数值数组，云滴半径（单位：μm）。可为标量或向量。
%
% 可选参数（Name-Value）：
%   'a'        : 有效半径 r_eff，默认 9（μm）。
%   'b'        : 形状参数 b，默认 1/9。
%   'normalize': 归一化方式，'continuous' 或 'grid'，默认 'continuous'。
%                - 'continuous' 使用解析式计算 const，使 ∫_0^∞ f(r) dr = 1。
%                - 'grid'       在给定 r 网格上数值归一化，使 trapz(r, f(r)) = 1。
%
% 输出：
%   f_r   : 在 r 处的概率密度（单位：μm^{-1}）。
%   const : 使用的归一化常数。
%
% 用途：
%   该分布通常用作对相函数/散射矩阵随半径的加权（体积平均或数浓度平均）。
%   若需要得到“数浓度谱”而非概率密度，可将 f_r 乘以总数浓度 N_d。
%
% 示例：
%   r = linspace(0.01, 50, 2000);
%   [f, c] = smmgDist_paper(r, 'a', 10, 'b', 1/9);
%   fprintf('连续归一化常数 const = %g\n', c);
%   fprintf('∫ f(r) dr ≈ %g\n', trapz(r, f));
%
% 注意：
%   - 输入与参数单位需一致（本函数以 μm 为单位）。
%   - 当 r 网格范围较窄或不包含接近 0 的小半径时，建议使用 'grid' 归一化。

% 参数解析
p = inputParser;
addRequired(p, 'r', @isnumeric);
addParameter(p, 'a', 9, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'b', 1/9, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'normalize', 'continuous', @(s) ismember(lower(s), {'continuous','grid'}));
parse(p, r, varargin{:});

% 取参数
r = p.Results.r(:);   % 列向量，便于向量化
a = p.Results.a;
b = p.Results.b;
normMode = lower(p.Results.normalize);

% 幂指数 pwr = (1 - 3b)/b
pwr = (1 - 3*b) / b;

% 先计算未归一化的核
kernel = r.^pwr .* exp(- r ./ (a*b));

switch normMode
    case 'continuous'
        % 连续解析归一化：
        % ∫_0^∞ r^pwr * exp(-r/(a*b)) dr = Γ(pwr+1) * (a*b)^(pwr+1)
        const = 1 / ( gamma(pwr + 1) * (a*b)^(pwr + 1) );
        f_r = const * kernel;

    case 'grid'
        % 在给定 r 网格上数值归一化
        % 避免 r 范围有限导致连续归一化偏差
        I = trapz(r, kernel);
        if I <= 0 || ~isfinite(I)
            error('smmgDist_paper:badGrid', '数值归一化失败：请检查 r 网格与参数 (a,b)。');
        end
        const = 1 / I;
        f_r = const * kernel;
end

% 将输出形状与输入一致
f_r = reshape(f_r, size(p.Results.r));

end