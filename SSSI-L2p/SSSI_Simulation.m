clc
clear
tic
%brainstorm %norgui
%%
algorithms = {'SSSI-L2p'};
ResultsLoading = [0];
WGN = 1; % Using white Gaussian Nosie (WGN = 1) or Human EEG noise (WGN  = 0);
LFNormlization = 0; % Whether normalize the LeadField Matrix
Uniform = 1; % Uniform/NonUniform Sources
VariousExtents = 0;
VariousSNRs = 1;
VariousSNIRs = 0;
VariousPatches = 0;
VariousCorrelation = 0;
VariousChannels = 0;
Test = 0;
if VariousExtents+VariousSNRs+VariousSNIRs+VariousPatches+VariousCorrelation+VariousChannels+Test ~= 1
    error('There will be one and only one scenario.');
end
% tic
%% Export the Channel,Cortex and HeadModel to workspace
if WGN
    channelselect=[1:32,34:42,44:64]; %Select EEG data
%     channelselect=[1:32 34:42 44:59 61:63];
else
    channelselect=[1:32 34:42 44:59 61:63]; % for Real noise simulations
end
% [sStudy, iStudy] = bst_get('StudyWithCondition', 'LiDaoli/Stim_2');
[sStudy, iStudy] = bst_get('StudyWithCondition','Subject01/Test');%GaussianSources');%Extents');%SNRs'); 
% [sStudy, iStudy] =bst_get('StudyWithCondition','Subject01/mind004_050924_median01_raw_clean_notch');
index = 1;
bst_call(@export_matlab, {char(sStudy.Data(1,index).FileName)},'data');
%=======================Import the LFM and Cortex==================================%
[sSurface, iSurface] = bst_get('SurfaceFileByType',[],'Cortex');
% [sSurface, iSurface] = bst_get('SurfaceFileByType',2,'Cortex');
bst_call(@export_matlab, {char(sSurface.FileName)},'Cortex');
Atlas = Cortex.Atlas(2);

[sHeadModel] = bst_get('HeadModelForStudy', iStudy);
bst_call(@export_matlab, {char(sHeadModel.FileName)},'model');
Gain=model.Gain(channelselect,:);
GridLoc=model.GridLoc;
GridOrient=model.GridOrient;

Gain = bst_gain_orient(Gain,GridOrient);
% load('MEGLeadfield.mat')
clear GridOrient
% load 'MEGLeadfield.mat';
% load 'CortexMEG.mat';
[nSensor,nSource] = size(Gain);
%% Reducing the leadfield matrix
%  L = Gain;
%  u = spm_svd(Gain*Gain');
%  L = u'*Gain;
% B = u'*B;
% clear u
%% Make output directory structure for imaging results (if doesn't already exist)
if VariousSNRs
   scenario = 'various SNRs';%'various SNRs';
   SNR1 = [0,3,5,10];
   SNIR1 = zeros(4,1)+5;
   condition = SNR1';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousSNIRs
   scenario = 'various SNIRs';
   SNR1 = zeros(4,1)+5;
   SNIR1 = [0,3,5,10];
   condition = SNIR1';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousExtents
   scenario = 'various extents';
   SNR1 = 5*ones(5,1);
   SNIR1 = 5*ones(5,1);
   condition = [1:5]';
   K = ones(5,1);
   DefinedArea = [2 5 10 18 32]'*1e-4*ones(1,max(K));%[0.5 4 8 14 22 32]'*1e-4*ones(1,2);% 38 48]'*1e-4;
elseif VariousChannels
   scenario = 'various channels';
   SNR1 = 5*ones(4,1);
   SNIR1 = 5*ones(4,1);
   condition = [62, 46, 32, 16]';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif Test
    algorithms = {'SSSI-L2p'}%'VSSI-sL2p','VSSI-L2p'};%'VB-SCCD','SSSI-L2p','VSSI-GGD'
    ResultsLoading = [0 0 0];
    scenario = 'test';
    SNR1 = 0;
    SNIR1 = 0;
    condition = [1];
    K = 1;
    DefinedArea = 5*1e-4;
end

outpath = 'E:\result\SSSI-L2p\';
for i = 1 : size(condition,1)
    path{i} = fullfile(outpath,scenario,'\',num2str(condition(i)));
    if ~exist(path{i})
        mkdir(path{i});
    end
     if ~exist([path{i} '\' 'metrics.mat'], 'file')
         metrics = [];
         save([path{i} '\' 'metrics.mat'], 'metrics')
     end
end

%% Iteration
    dim = 0;
    Miter = 50 - dim;   
    Eccentricity = sqrt(sum(GridLoc.^2,2));
    Ec = find(Eccentricity > 70*1e-3);
    EstimatedArea_Mean = zeros(Miter,5);
for iter = 1:Miter    
%     ind = randperm(size(Gain,2)); 
    ind = randperm(numel(Ec));
    Thr = zeros(size(condition,1),numel(algorithms));
for iteration = 1:size(condition,1)
    fprintf('iter = %g, iteration = %g\n', iter,iteration)
    savepath = path{iteration};
    load ([path{iteration} '\' 'metrics'])
    SD = []; DLE = []; RMSE = []; AUC = []; PRE = []; REC = [];
    if any(ResultsLoading)
        SD = metrics.SD(iter,:); DLE = metrics.DLE(iter,:); RMSE = metrics.RMSE(iter,:); AUC = metrics.AUC(iter,:); PRE = metrics.PRE(iter,:); REC = metrics.REC(iter,:); nRMSE = metrics.nRMSE(iter,:); SE = metrics.SE(iter,:);
    end
%     if VariousSNIRs && (iteration==2 || iteration==3 || iteration==4)
%         continue;
%     end
%% Generate Simulated EEG Data
SNR = SNR1(iteration);
SNIR = SNIR1(iteration);
seedvox = Ec(ind(1:K(iteration)));
tau = [0.1 0.35 0.5 0.6];omega = [0.1 0.15 0.15 0.15];%[0.07 0.05 0.10 0.10];%[0.035 0.035 0.035 0.035];%*max(Time(Activetime));
%tau = [0.1 0.2 0.5 0.6];omega = [0.1 0.15 0.15 0.15];%[0.07 0.05 0.10 0.10];%[0.035 0.035 0.035 0.035];%*max(Time(Activetime));
f = [10 11 8 9];%10*ones(1,4);%5 (Hz);
Amp = 1e-8;
StimTime = find(abs(data.Time) == min(abs(data.Time)));
TimeLen = 300;
Time = data.Time(StimTime-0.5*TimeLen:StimTime+0.5*TimeLen-1);

OPTIONS.DefinedArea    = DefinedArea(iteration,:);
OPTIONS.seedvox        = seedvox;
OPTIONS.frequency      = f;
OPTIONS.tau            = tau;
OPTIONS.omega          = omega;
OPTIONS.Amp            = Amp;
OPTIONS.GridLoc        = GridLoc;
if VariousPatches
OPTIONS.MixMatrix = [1    0    0     0;
                     0    1    0     0;
                     0    0    1     0;
                     0    0    0     1;
                     0    0    .5     .5];
elseif VariousCorrelation
    xi = condition(iteration);
%     OPTIONS.MixMatrix = [1-xi  xi   0    0;
%                          0     1    0    0;
%                          0  0   0    1;
%                          0  0   1    0];
OPTIONS.MixMatrix = [1            0            0          0;
                     xi     sqrt(1-xi^2)       0          0;
                     xi           0        sqrt(1-xi^2)   0;
                     0            0            1          0];
end
OPTIONS.uniform       = Uniform;
OPTIONS.WGN           = WGN;
OPTIONS.SNR           = SNR;
OPTIONS.SNIR          = SNIR;
OPTIONS.ar            = 0;
OPTIONS.params(:,:,1) = [ 0.8    0    0 ;
                            0  0.9  0.5 ;
                          0.4    0  0.5];

OPTIONS.params(:,:,2) = [-0.5    0    0 ;
                            0 -0.8    0 ;
                            0    0 -0.2];

OPTIONS.noisecov      = [ 0.3    0    0 ;
                            0    1    0 ;
                            0    0  0.2];

if ~any(ResultsLoading)
    [Data,s_real,Result] = Simulation_Data_Generate (Gain,Cortex,Time,OPTIONS);
    ActiveVoxSeed = Result.ActiveVoxSeed;
else
    % %================= Recording from previous Results ====================%
    load ([savepath '\' 'result' num2str(iter+dim)]);
    Data = Result.B;
    s_real = Result.real;
    seedvox = Result.seedvox;
    Result.B = Data;
    [~, VertArea] = tess_area(Cortex.Vertices, Cortex.Faces);
    AreaDef = DefinedArea(iteration,:);
    ActiveVox = [];
    for k = 1:numel(seedvox)
        ActiveVoxSeed{k} = PatchGenerate(seedvox(k),Cortex.VertConn,VertArea,AreaDef(k));
        ActiveVox = union(ActiveVoxSeed{k},ActiveVox);
    end
    StimTime = size(Data,2)/2 + 1;
end
% %=======================================================================%
fprintf('Actual SNR is %g\n',20*log10(norm(Gain*s_real,'fro')/norm(Data-Gain*s_real,'fro')));
corr(s_real(seedvox,StimTime+1:end)','type','Pearson')
%% Data scale
Scale = 1;
if Scale
    ScaleType = 0;  % ScaleType = 1, Simultaneously scale the MEG data and leadfiled matrix; ScaleType = 0, only scale the MEG data
    ratio = 1e-6%1e-6;%Ratio(iiter);
    if ScaleType == 0
        B = Data./ratio;
        Gain_scale = Gain;
    else
        B  = Data./ratio;
        Gain_scale = Gain./ratio;
        ratio = 1;
    end
else
    ratio = 1;
    B = Data;
    Gain_scale = Gain;
end
%% 
%======================================%
%         Leadfield Matrix normalization
%=====================================%
if LFNormlization
    LfvW = sqrt(sum(Gain_scale.^2,1).^0.3);
    Gain_scale = Gain.*kron(ones(nSensor,1),1./LfvW);
end

%% Whiten measurements and lead field matrix
Bpre = B(:,1:StimTime);
Cov_n = NoiseEstimate(B,StimTime);
NoiseMethod = 'median';%'reg';%'none';%{'shrink', 'reg', 'diag', 'none', 'median'};
FourthMoment = Bpre*Bpre'./StimTime;nSamples = StimTime;
NoiseReg = 0.1;
[Cov,W] = truncate_and_regularize_covariance(Cov_n,NoiseMethod,'EEG',NoiseReg,FourthMoment,nSamples);
L = W*Gain_scale;
B = W*B;
if ~any(ResultsLoading)
    Result.Whiter = W;
end
clear W;
%%  SVD for TBFs
[Dic] = TBFSelection(B,0);%,'threshold','Permutation');'Kaiser'

MetricsInterval = [];

%% Channel select
if VariousChannels
    if ~any(ResultsLoading)
        cnumber = condition(iteration);
        load(['C:\Users\guest0\Desktop\GGD-fig\channel\',num2str(cnumber),'c.mat']);
        channels(channels>33 & channels<43) = channels(channels>33 & channels<43) - 1;
        channels(channels>43) = channels(channels>43) - 2;
        B = B(channels,:);
        L = L(channels,:);
        Result.channels = channels;
    else
        channels = Result.channels;
        B = B(channels,:);
        L = L(channels,:);

    end
end

%% SSSI-L2p
Weight = logspace(-4,0,20);
% Sourcenorm = zeros(10,1);variationnorm = zeros(10,1);
variation = 'Variation';%'Sparse+Variation';%'Sparse';% 'Laplacian';%'Laplacian+Variation';%'Sparse+Laplacian';%'Laplacian';%
opts.sparse = 0.5;%Weight(14);%Weight(8);%0.05;%%0.15;
% opts.laplacian = 0.8;
if any(strcmpi('SSSI-L2p', algorithms))
    MethodIndex = find(strcmpi('SSSI-L2p', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    if ~ResultsLoading(MethodIndex)
        [S_sssi] = SSSI_L2p_ADMM(B*Dic',L,Cortex.VertConn,opts.sparse,'transform',variation,'p',0.8,'tol',1e-6);
        S_sssil2p = S_sssi{end}*Dic;
        if LFNormlization
            S_sssil2p = repmat(1./LfvW',1,size(S_sssil2p,2)).*S_sssil2p;
        end
        S_sssil2p = S_sssil2p*ratio;
    else
        S_sssil2p = Result.SSSIL2p;
    end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_sssil2p(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_sssil2p(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_sssil2p,'fro')^2/norm(B*ratio,'fro')^2;
    Result.SSSIL2p = S_sssil2p;
end 

% Source{1} = s_real;
% Comment{1} = ['GroundTruth-',num2str(seedvox)];
% for j = 1:size(S_sssi,2)
%     Source{j+1} = repmat(sqrt(sum(S_sssi{j}*Dic.^2,2)),1,300);
%     Comment{j+1} = ['SSSI-L2p-',num2str(ThresholdSelect(S_sssi{j}*Dic)*100),'%'];
% end
% % for jj = 1:size(S_L2p,2)
% %     Source{j+jj} = repmat(sqrt(sum(S_L2p{jj}*Dic.^2,2)),1,300);
% %     Comment{j+jj} = ['VSSI-L2p-',num2str(ThresholdSelect(S_L2p{jj}*Dic)*100),'%'];
% % end
% temp = 1;
% %load Node.mat;
% %bst_call(@tree_callbacks, bstNodes, 'doubleclick');
% for j = 1:size(Source,2)
%     bst_call(@export_matlab, char(sStudy.Result(size(Source,2)*(temp-1)+j+1).FileName),'result');
%     result.ImageGridAmp=Source{j};
%     result.Comment = Comment{j};
%     Node_Import_PS(char(sStudy.Result(size(Source,2)*(temp-1)+j+1).FileName),iStudy,'result');
% end
% for j = 1:size(Source,2)
%     view_surface_data([],char(sStudy.Result(size(Source,2)*(temp-1)+j+1).FileName));
% end
% bst_memory('UnloadAll', 'Forced');
% bst_progress('stop'); 

continue;
%% Save Results on Disk
method = 1:numel(algorithms);
metrics.AUC(dim+iter,method) = AUC;%(method);
metrics.SD(dim+iter,method) = SD;%(method); 
metrics.DLE(dim+iter,method) = DLE;%(method);
metrics.RMSE(dim+iter,method) = RMSE;%(method);
metrics.nRMSE(dim+iter,method) = nRMSE;%(method);
metrics.SE(dim+iter,method) = SE;%(method);
metrics.PRE(dim+iter,method) = PRE;
metrics.REC(dim+iter,method) = REC;

save_file_name=[savepath '\' 'result' num2str(dim+iter)];
save (save_file_name,'Result');
save([savepath '\' 'metrics.mat'], 'metrics')

end
end
toc
runningtime = toc;