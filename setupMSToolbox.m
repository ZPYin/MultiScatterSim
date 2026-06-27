global MS_ENVS;

MS_ENVS.RootPath = fileparts(mfilename('fullpath'));
MS_ENVS.Author = 'Zhenping Yin; Chengwei Luo';
MS_ENVS.Email = 'zp.yin@whu.edu.cn';
MS_ENVS.UpdateTime = '2023-05-22';
MS_ENVS.Version = '0.0.1';

addpath(genpath(fullfile(MS_ENVS.RootPath, 'lib')));
addpath(genpath(fullfile(MS_ENVS.RootPath, 'include')));
addpath(genpath(fullfile(MS_ENVS.RootPath, 'tests')));