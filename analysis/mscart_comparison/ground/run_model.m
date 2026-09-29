function run_model(modelRoot,nPhoton,seeds)
% Run the formal model in the delivered low/high ground cases.
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
s=load(fullfile(here,'optics.mat'),'optics','cloud','lidar');
levels={'low','high'};allRows=cell(2*numel(seeds)*2,1);entry=0;
for scale=1:2
    for seed=seeds(:)'
        for fov=[.5,2]
            optics=s.optics;cloud=s.cloud;lidar=s.lidar;
            optics.muExt=optics.muExt*scale;optics.muSca=optics.muSca*scale;
            optics.cloudScatMatrix.muExt=optics.muExt;optics.cloudScatMatrix.muSca=optics.muSca;
            cloud.N_d=cloud.N_d*scale;lidar.FOV=fov;lidar.rBeam=0;
            fprintf('START ground %s / %g mrad / seed %d / N %d\n',levels{scale},fov,seed,nPhoton);
            caseTimer=tic;
            [sv,d]=MSMCSimGPUBackend(lidar,cloud,'preparedOptics',optics,'nPhoton',nPhoton, ...
                'seed',seed,'precision','double','batchSize',min(nPhoton,262144), ...
                'maxScatterOrder',1024,'diagnosticsLevel','full','useRussianRoulette',false, ...
                'fallbackToCPU',false,'geometryMode','mscart','exactTimeWindowStop',true);
            assert(d.truncatedPhotons==0);
            fprintf('DONE ground %s / %g mrad / seed %d: %.3f s\n',levels{scale},fov,seed,toc(caseTimer));
            conversion=1/(nPhoton*pi*(lidar.rTel*1e-3)^2*2*pi);
            raw=sv*conversion;nth=d.SVByOrder.*conversion;
            T=table(repmat(string(levels{scale}),60,1),repmat(fov,60,1),repmat(seed,60,1), ...
                ((1:60)'-.5)*5,(raw(1,:)'+raw(2,:)')/2,(raw(1,:)'-raw(2,:)')/2, ...
                'VariableNames',{'extinction','fov_full_mrad','repeat','depth_m','parallel','perpendicular'});
            for k=1:4
                if k<4,w=nth(:,:,k);else,w=raw-sum(nth(:,:,1:3),3);end
                T.(sprintf('order%d_parallel',k))=(w(1,:)'+w(2,:)')/2;
                T.(sprintf('order%d_perpendicular',k))=(w(1,:)'-w(2,:)')/2;
            end
            entry=entry+1;allRows{entry}=T;
        end
    end
end
writetable(vertcat(allRows{:}),fullfile(here,'model_results.csv'));
fprintf('Saved ground/model_results.csv\n');
end
