function [figPos] = subfigPos(pos, nRow, nColumn, xPad, yPad)
% SUBFIGPOS 计算多子图布局中每个子图的标准化位置坐标
%
% 功能描述：
%   这是一个用于MATLAB图形显示的实用工具函数，专门用于计算多子图布局中
%   每个子图的精确位置。该函数可以自动计算在给定的总体区域内，按照指定
%   的行列数和间距要求，每个子图应该占据的标准化坐标位置。
%
%   主要特点：
%   - 支持任意行列数的子图网格布局
%   - 可自定义水平和垂直间距
%   - 返回MATLAB标准的标准化坐标系统 [0,1]
%   - 子图编号顺序：从左到右，从上到下
%
% 使用方法：
%    [figPos] = subfigPos(pos, nRow, nColumn, xPad, yPad)
%
% 输入参数：
%    pos: 4元素数组
%        总体布局区域的标准化坐标 [left, bottom, width, height]
%        其中：left-左边距，bottom-下边距，width-总宽度，height-总高度
%    nRow: 整数
%        子图网格的总行数 (默认值: 1)
%    nColumn: 整数  
%        子图网格的总列数 (默认值: 1)
%    xPad: 数值
%        子图间的水平间距（标准化坐标）
%    yPad: 数值
%        子图间的垂直间距（标准化坐标）
%
% 输出参数：
%    figPos: 矩阵 [nRow*nColumn × 4]
%        每个子图的位置坐标矩阵，每行包含一个子图的位置信息
%        格式：[xPos, yPos, width, height]
%        子图编号顺序：第1个子图位于左上角，按行优先顺序排列
%        (从左到右，从上到下)
%
% 应用场景：
%   - 激光雷达数据的多参数同时显示
%   - 散射特性的多角度对比分析  
%   - 仿真结果的网格化展示
%   - 科学论文图表的标准化布局
%
% 历史记录：
%    - 2018-11-09: 首版 by Zhenping
%    - 2019-12-15: 增加x-y间距设置功能
%
% 作者: zhenping@tropos.de

% 检查可选参数是否存在，如不存在则设置默认值
if ~ exist('xPad', 'var')
    xPad = 0;    % 默认水平间距为0
end

if ~ exist('yPad', 'var')
    yPad = 0;    % 默认垂直间距为0
end

% 初始化输出矩阵：存储所有子图的位置信息
% 矩阵大小：[总子图数 × 4]，每行存储一个子图的[x, y, width, height]
figPos = zeros(nRow * nColumn, 4);

% ========== 计算边距参数 ==========
% 基于输入的总体区域pos，计算各方向的边距
leftMargin = pos(1);                    % 左边距
rightMargin = 1 - pos(1) - pos(3);     % 右边距 = 1 - 左边距 - 总宽度
topMargin = 1 - pos(2) - pos(4);       % 上边距 = 1 - 下边距 - 总高度  
baseMargin = pos(2);                    % 下边距

% ========== 计算单个子图的尺寸 ==========
% 考虑边距和子图间距，计算每个子图的实际宽度和高度
% 可用宽度 = 总宽度 - 左右边距 - (列数-1)×水平间距，然后除以列数
widthCol = (1 - leftMargin - rightMargin - (nColumn - 1) * xPad) / nColumn;
% 可用高度 = 总高度 - 上下边距 - (行数-1)×垂直间距，然后除以行数
widthRow = (1 - topMargin - baseMargin - (nRow - 1) * yPad) / nRow;

% ========== 计算每个子图的具体位置 ==========
% 双重循环遍历所有子图位置（行优先顺序）
for iRow = 1:nRow
    for iCol = 1:nColumn
        % 计算当前子图的x坐标（左下角）
        % x位置 = 左边距 + (列索引-1) × (子图宽度 + 水平间距)
        xPos = leftMargin + (iCol - 1) * (widthCol + xPad);
        
        % 计算当前子图的y坐标（左下角）
        % y位置 = 下边距 + (总行数-当前行数) × (子图高度 + 垂直间距)
        % 注意：(nRow - iRow)实现了从上到下的排列顺序
        yPos = baseMargin + (nRow - iRow) * (widthRow + yPad);
        
        % 将计算得到的位置信息存储到输出矩阵中
        % 子图索引 = (行索引-1) × 列数 + 列索引
        figPos((iRow - 1)*nColumn + iCol, :) = [xPos, yPos, widthCol, widthRow];
    end 
end

end
