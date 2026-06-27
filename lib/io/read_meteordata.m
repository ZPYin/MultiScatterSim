function [temp, pres, relh, meteor_time] = read_meteordata(measTime, altitude, varargin)
% READ_METEORDATA 根据输入的气象数据类型读取气象数据
%
% 功能描述：
%   这是一个通用的气象数据读取函数，支持多种气象数据源，包括：
%   - GDAS1 (Global Data Assimilation System) 全球数据同化系统
%   - Radiosonde 无线电探空仪数据
%   - ERA-5 欧洲中期天气预报中心再分析数据
%   - Standard Atmosphere 标准大气模型
%   
%   该函数通过插值方法将气象数据映射到指定的高度层，为激光雷达仿真
%   提供必要的大气环境参数（温度、压力、相对湿度）。
%
% 使用示例：
%    % 用例1: 读取GDAS1数据
%    [temp, pres, relh] = read_meteordata(measTime, altitude, 'meteor_data', 'GDAS1', 'station', 'wuhan', 'GDAS1Folder', '/GDAS1');
%    % 用例2: 读取无线电探空数据
%    [temp, pres, relh] = read_meteordata(measTime, altitude, 'meteor_data', 'Radiosonde', 'RadiosondeFolder', '/Radiosonde', 'station', 'wuhan');
%    % 用例3: 读取ERA-5再分析数据
%    [temp, pres, relh] = read_meteordata(measTime, altitude, 'meteor_data', 'ERA-5', 'ERA5Folder', '/ERA-5', 'station', 'wuhan');
%
% 输入参数：
%    measTime: datenum格式
%        测量时间 (UTC时间)
%    altitude: 数组
%        海拔高度 [m]
%
% 可选参数：
%    meteor_data: 字符串
%        气象数据源类型
%        'GDAS1' (默认), 'Radiosonde', 'ERA-5', 'standard_atmosphere'
%    station: 字符串
%        观测站标识 (默认: 'wuhan')
%    GDAS1Folder: 字符串
%        GDAS1数据文件夹路径
%    RadiosondeFolder: 字符串
%        无线电探空数据文件夹路径
%    RadiosondeType: 整数
%        无线电探空文件类型
%        - 1: MOSAiC项目探空文件格式
%        - 2: MUA项目探空文件格式 (默认)
%    ERA5Folder: 字符串
%        ERA-5数据文件夹路径
%
% 输出参数：
%    temp: 数组
%        各高度层的温度 [°C]
%    pres: 数组
%        各高度层的压力 [hPa]
%    relh: 数组
%        各高度层的相对湿度 [%]
%    meteor_time: 数值
%        气象数据的时间戳 (datenum格式)
%
% 历史记录：
%    2020-05-28. 首版 by Zhenping
% 联系方式：
%    zhenping@tropos.de

% 创建输入参数解析器
p = inputParser;
p.KeepUnmatched = true;  % 保留未匹配的参数

% 添加必需参数
addRequired(p, 'measTime', @isnumeric);    % 测量时间（datenum格式）
addRequired(p, 'altitude', @isnumeric);    % 高度数组

% 添加可选参数及其默认值
addParameter(p, 'meteor_data', 'GDAS1', @ischar);        % 气象数据源类型
addParameter(p, 'station', 'wuhan', @ischar);            % 观测站标识
addParameter(p, 'GDAS1Folder', '', @ischar);             % GDAS1数据文件夹
addParameter(p, 'RadiosondeFolder', '', @ischar);        % 探空数据文件夹
addParameter(p, 'RadiosondeType', 2, @isnumeric);        % 探空文件类型
addParameter(p, 'ERA5Folder', '', @ischar);              % ERA5数据文件夹

% 解析输入参数
parse(p, measTime, altitude, varargin{:});

% 初始化输出变量为NaN（与高度数组同尺寸）
temp = NaN(size(altitude));        % 温度数组
pres = NaN(size(altitude));        % 压力数组
relh = NaN(size(altitude));        % 相对湿度数组
meteor_time = NaN(1);              % 气象数据时间戳

% 根据气象数据源类型进行分支处理
switch lower(p.Results.meteor_data)

case 'gdas1'
    % ========== GDAS1数据处理 ==========
    % GDAS1 (Global Data Assimilation System) 是NCEP的全球数据同化系统
    % 提供全球范围的气象再分析数据，时间分辨率为6小时
    
    % 读取GDAS1气象廓线数据
    [gdas1_alt, gdas1_temp, gdas1_pres, gdas1_relh, gdas1_file] = read_gdas1(measTime, ...
        p.Results.station, p.Results.GDAS1Folder);

    % 检查数据有效性
    if isempty(gdas1_alt) || (~ (sum(~ isnan(gdas1_temp)) >= 2))
        % 如果没有数据或有效数据点少于2个，则返回NaN
        return;
    end

    % 提取气象数据时间戳并进行插值
    meteor_time = gdas1FileTimestamp(basename(gdas1_file));
    temp = interp1(gdas1_alt, gdas1_temp, altitude);    % 温度插值
    pres = interp1(gdas1_alt, gdas1_pres, altitude);    % 压力插值
    relh = interp1(gdas1_alt, gdas1_relh, altitude);    % 相对湿度插值

case 'radiosonde'
    % ========== 无线电探空数据处理 ==========
    % Radiosonde是通过气球携带的探测仪器进行的大气垂直探测
    % 提供高精度的温度、压力、湿度垂直廓线
    
    % 搜索对应时间的探空文件
    sondeFile = radiosonde_search(fullfile(p.Results.RadiosondeFolder, p.Results.station), measTime, p.Results.RadiosondeType);
    
    % 读取探空数据
    [rs_alt, rs_temp, rs_pres, rs_relh, rs_time] = ...
        read_radiosonde(sondeFile, p.Results.RadiosondeType);

    % 检查数据有效性
    if isempty(rs_alt) || (~ (sum(~ isnan(rs_temp)) >= 2))
        % 如果没有数据或有效数据点少于2个，则返回NaN
        return;
    end

    % 按高度升序排列测量数据（探空数据可能不是按高度排序的）
    [rs_alt, sortIndxAlt] = sort(rs_alt);
    rs_temp = rs_temp(sortIndxAlt);
    rs_pres = rs_pres(sortIndxAlt);
    rs_relh = rs_relh(sortIndxAlt);

    % 去除同一高度的重复测量值（保证插值函数的单调性）
    [rs_alt, iUniq, ~] = unique(rs_alt);
    rs_temp = rs_temp(iUniq);
    rs_pres = rs_pres(iUniq);
    rs_relh = rs_relh(iUniq);

    % 设置时间戳并进行插值
    meteor_time = rs_time;
    temp = interp1(rs_alt, rs_temp, altitude);    % 温度插值
    pres = interp1(rs_alt, rs_pres, altitude);    % 压力插值
    relh = interp1(rs_alt, rs_relh, altitude);    % 相对湿度插值

case 'era-5'
    % ========== ERA-5再分析数据处理 ==========
    % ERA-5是欧洲中期天气预报中心(ECMWF)的第五代大气再分析数据
    % 提供1979年至今的全球高分辨率气象数据，时间分辨率为1小时
    
    % 读取ERA-5气象廓线数据
    [ERA5_alt, ERA5_temp, ERA5_pres, ERA5_relh, ERA5_file] = read_ERA5(measTime, ...
        p.Results.station, p.Results.ERA5Folder);

    % 检查数据有效性
    if isempty(ERA5_alt) || (~ (sum(~ isnan(ERA5_temp)) >= 2))
        % 如果没有数据或有效数据点少于2个，则返回NaN
        return;
    end

    % 提取气象数据时间戳并进行插值
    meteor_time = ERA5FileTimestamp(basename(ERA5_file));
    temp = interp1(ERA5_alt, ERA5_temp, altitude);    % 温度插值
    pres = interp1(ERA5_alt, ERA5_pres, altitude);    % 压力插值
    relh = interp1(ERA5_alt, ERA5_relh, altitude);    % 相对湿度插值

case 'standard_atmosphere'
    % ========== 标准大气模型处理 ==========
    % 使用1976年美国标准大气模型作为默认的大气参数
    % 适用于没有实际气象观测数据时的理论计算
    
    % 调用标准大气模型函数（高度单位转换为km）
    [sa_alt, ~, ~, sa_temp, sa_pres] = atmo(max(altitude/1000)+1, 0.03, 1);

    % 插值计算各高度层的气象参数
    pres = interp1(sa_alt, sa_pres / 1e2, altitude / 1000);        % 压力：Pa转hPa
    temp = interp1(sa_alt, sa_temp - 273.17, altitude / 1000);     % 温度：K转°C
    relh = NaN(size(temp));                                        % 标准大气模型不提供湿度

otherwise
    % 未知的气象数据类型，抛出错误
    error('Unknown meteorological data: %s', p.Results.meteor_data);

end

end
