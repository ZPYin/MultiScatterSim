function [samples] = rejectSampling(x, y, nSamples, varargin)
% REJECTSAMPLING 拒绝采样算法 - 基于概率密度函数的随机采样
%
% 功能描述：
%    该函数实现拒绝采样（Rejection Sampling）算法，也称为接受-拒绝方法。
%    根据给定的概率密度函数（由x和y定义）生成符合该分布的随机样本。
%    这是蒙特卡罗仿真中重要的采样技术，用于从复杂分布中生成随机数。
%
% 算法原理：
%    拒绝采样是一种通用的随机采样方法：
%    1. 在定义域内均匀生成候选样本
%    2. 计算候选样本处的概率密度值
%    3. 以该密度值与最大密度值的比值作为接受概率
%    4. 通过随机数判断是否接受该样本
%    5. 重复直到获得足够数量的有效样本
%
% 使用方法：
%    [samples] = rejectSampling(x, y, 1000)              % 生成1000个样本
%    [samples] = rejectSampling(theta, phase_func, 500)  % 从相函数中采样散射角
%
% 输入参数：
%    x: 数值型数组
%        自变量数组，定义概率密度函数的定义域
%        例如：散射角、粒径、波长等物理量
%    y: 数值型数组
%        因变量数组，表示概率密度函数值
%        必须与x数组长度相同，且为非负值
%    nSamples: 数值型标量
%        需要生成的样本数量
%        正整数，表示期望的采样点个数
%
% 可选参数：
%    (当前版本暂无可选参数，保留扩展接口)
%
% 输出参数：
%    samples: 数值型数组
%        生成的随机样本数组，长度为nSamples
%        样本服从输入的概率密度函数分布
%
% 应用场景：
%    1. 蒙特卡罗光子传输仿真中的散射角采样
%    2. 云滴粒径谱的随机采样
%    3. 散射相函数的方向采样
%    4. 复杂概率分布的随机数生成
%    5. 激光雷达信号仿真中的随机过程建模
%
% 算法特点：
%    - 适用于任意形状的概率密度函数
%    - 实现简单，概念清晰
%    - 采样效率取决于分布的形状
%    - 对于峰值较尖的分布效率较低
%
% 注意事项：
%    - 输入的y数组应为归一化的概率密度函数
%    - 算法效率与分布形状相关，平坦分布效率更高
%    - 对于高维分布，建议考虑其他采样方法
%
% 历史记录：
%    2023-05-24: 首次版本 by Zhenping
% 作者联系方式：zp.yin@whu.edu.cn

% 参数解析器设置
p = inputParser;
p.KeepUnmatched = true;  % 保留未匹配的参数

% 定义必需参数
addRequired(p, 'x', @isnumeric);         % 必需参数：自变量数组
addRequired(p, 'y', @isnumeric);         % 必需参数：概率密度函数值
addRequired(p, 'nSamples', @isnumeric);  % 必需参数：样本数量

% 解析输入参数
parse(p, x, y, nSamples, varargin{:});

% 初始化采样计数器和样本数组
sampleCounter = 0;    % 已接受的样本计数
samples = [];         % 存储接受的样本

% 拒绝采样主循环：持续采样直到获得足够数量的样本
while (sampleCounter < nSamples)
    % 步骤1：在定义域内均匀生成1000个候选样本
    % 使用均匀分布在[x(1), x(end)]区间内生成随机数
    xSamples = rand(1, 1e3) * (x(end) - x(1)) + x(1);
    
    % 步骤2：找到概率密度函数的最大值，用于归一化
    maxY = max(y);
    
    % 步骤3：计算每个候选样本的接受概率
    % 通过插值计算候选样本处的概率密度值，然后除以最大值进行归一化
    % 生成[0,1]区间的随机数，与归一化概率密度比较决定是否接受
    accepted = (rand(1, 1e3)) < (interp1(x, y, xSamples) / maxY);
    
    % 步骤4：更新计数器并保存被接受的样本
    sampleCounter = sampleCounter + sum(accepted);
    samples = cat(2, samples, xSamples(accepted));
end

% 步骤5：截取所需数量的样本（可能会超出nSamples）
samples = samples(1:nSamples);

end