function run_model(modelRoot,nPhoton,seeds)
% Run the formal model in the delivered Hu2001 satellite case.
here=fileparts(mfilename('fullpath'));
% Run directly from the MATLAB editor, or pass [] to locate this project.
if nargin<1 || isempty(modelRoot)
    modelRoot=fileparts(fileparts(fileparts(here)));
end
if nargin<2 || isempty(nPhoton),nPhoton=2000000;end
if nargin<3 || isempty(seeds),seeds=[20260929,20260930,20261001];end
oldPath=path;
modelRoot=strrep(char(modelRoot),'/',filesep);
cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(genpath(fullfile(modelRoot,'lib')));
addpath(fullfile(modelRoot,'lib','MC_GPU'),'-begin');
clear MSMCSimGPUBackend simulateLayeredBatchGPU freeRunningGPU;
assert(strcmp(which('MSMCSimGPUBackend'),fullfile(modelRoot,'lib','MC_GPU','MSMCSimGPUBackend.m')));
s=load(fullfile(here,'optics.mat'),'optics');
cloud=struct('clBase',704700,'clDh',300,'nClLayers',1,'kappa',9,'clReff',6,'N_d',1,'atmExt',0);
lidar=struct('lambda',532,'rBeam',0,'divBeam',0,'rTel',1,'binWidth',20,'nBins',100, ...
    'distBeamTel',0,'FOV',.001,'SBeam',[1;1;0;0]);
rows=cell(numel(seeds)*3,1);entry=0;depth=((1:100)'-.5)*3;
for seed=seeds(:)'
    for fov=[.001,.13,13]
        lidar.FOV=fov;
        fprintf('START satellite / %g mrad / seed %d / N %d\n',fov,seed,nPhoton);
        caseTimer=tic;
        [sv,d]=MSMCSimGPUBackend(lidar,cloud,'preparedOptics',s.optics,'nPhoton',nPhoton, ...
            'seed',seed,'precision','double','batchSize',min(nPhoton,131072), ...
            'maxScatterOrder',256,'diagnosticsLevel','full','useRussianRoulette',false, ...
            'fallbackToCPU',false,'geometryMode','mscart','proposalTailDegrees',90);
        assert(d.truncatedPhotons==0);
        fprintf('DONE satellite / %g mrad / seed %d: %.3f s\n',fov,seed,toc(caseTimer));
        conversion=1/(nPhoton*pi*(lidar.rTel*1e-3)^2*2*pi);
        raw=sv.*(cloud.clBase+depth)'.^2*conversion;
        nth=d.SVByOrder.*reshape((cloud.clBase+depth).^2,1,100,1).*conversion;
        T=table(repmat("hu2001",100,1),repmat(fov,100,1),repmat(seed,100,1),depth, ...
            (raw(1,:)'+raw(2,:)')/2,(raw(1,:)'-raw(2,:)')/2, ...
            'VariableNames',{'extinction','fov_full_mrad','repeat','depth_m','parallel','perpendicular'});
        for k=1:4
            if k<4,w=nth(:,:,k);else,w=raw-sum(nth(:,:,1:3),3);end
            T.(sprintf('order%d_parallel',k))=(w(1,:)'+w(2,:)')/2;
            T.(sprintf('order%d_perpendicular',k))=(w(1,:)'-w(2,:)')/2;
        end
        entry=entry+1;rows{entry}=T;
    end
end
writetable(vertcat(rows{:}),fullfile(here,'model_results.csv'));
fprintf('Saved satellite/model_results.csv\n');
end
