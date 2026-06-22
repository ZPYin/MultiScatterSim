function [refIdxReal, refIdxImg] = refIdxWater(inWL, varargin)
% REFIDXWATER 水的复折射率计算函数
%
% 功能描述：
%    该函数根据给定的波长计算水的复折射率，包括实部和虚部。
%    通过插值预存的水折射率数据库，为Mie散射计算和激光雷达仿真
%    提供准确的水滴光学参数。这是云光学特性计算的基础函数。
%
% 物理意义：
%    水的复折射率 m = n + ik，其中：
%    - 实部 n：决定光在水中的传播速度和折射特性
%    - 虚部 k：决定水对光的吸收特性
%    复折射率是Mie散射理论计算的核心输入参数，直接影响散射和消光特性。
%
% 使用方法：
%    [refIdxReal, refIdxImg] = refIdxWater(532)              % 532nm波长
%    [refIdxReal, refIdxImg] = refIdxWater(355, 'method', 'spline')  % 使用样条插值
%
% 输入参数：
%    inWL: 数值型标量或数组
%        输入波长，单位：纳米 (nm)
%        支持单个波长或波长数组
%
% 可选参数：
%    method: 字符串
%        插值方法 (默认值: 'linear')
%        可选：'linear', 'nearest', 'spline', 'cubic' 等
%        不同方法影响插值精度和计算速度
%
% 输出参数：
%    refIdxReal: 数值型
%        水的折射率实部 (无量纲)
%        典型值范围：1.33-1.34 (可见光波段)
%    refIdxImg: 数值型
%        水的折射率虚部 (无量纲)
%        典型值范围：10^-9 到 10^-6 (可见光波段)
%
% 应用场景：
%    1. Mie散射理论计算中的复折射率输入
%    2. 云滴和雨滴的光学特性计算
%    3. 激光雷达信号仿真中的散射参数
%    4. 水云的消光和散射系数计算
%    5. 偏振激光雷达的散射矩阵计算
%
% 数据来源：
%    函数依赖于 'refractive_index_water.mat' 数据文件，
%    该文件包含不同波长下水的折射率实验测量数据。
%
% 注意事项：
%    - 输入波长必须在数据库覆盖范围内，否则会报错
%    - 波长单位为纳米，函数内部会转换为微米进行插值
%    - 建议在可见光和近红外波段使用，精度最高
%
% 使用示例：
%    % 计算532nm激光的水折射率
%    [n, k] = refIdxWater(532);
%    m_water = n + 1i*k;  % 复折射率
%
% 历史记录：
%    2022-06-29: 首次版本 by Zhenping
% 作者联系方式：zp.yin@whu.edu.cn

% 参数解析器设置
p = inputParser;
p.KeepUnmatched = true;  % 保留未匹配的参数

% 定义必需参数和可选参数
addRequired(p, 'inWL', @isnumeric);                    % 必需参数：输入波长(nm)
addParameter(p, 'method', 'linear', @ischar);          % 可选参数：插值方法，默认线性插值

% 解析输入参数
parse(p, inWL, varargin{:});

% 波长单位转换：纳米 -> 微米
% 数据库中的波长以微米为单位存储
inWL = inWL / 1e3;   % (micron)

% 加载水的折射率数据库
% 该.mat文件包含wavelength、refractive_index_real、refractive_index_imaginary三个变量
res = load('refractive_index_water.mat');

% 检查输入波长是否在数据库覆盖范围内
% 超出范围会导致插值失败，因此需要提前检查
if (inWL > res.wavelength(end)) || (inWL < res.wavelength(1))
    error('Wavelength is out of range.');
end

% 使用一维插值计算折射率实部
% interp1函数根据已知的波长-折射率数据点进行插值
refIdxReal = interp1(res.wavelength, res.refractive_index_real, inWL, p.Results.method);

% 使用一维插值计算折射率虚部
% 虚部决定了水对光的吸收特性
refIdxImg = interp1(res.wavelength, res.refractive_index_imaginary, inWL, p.Results.method);

end