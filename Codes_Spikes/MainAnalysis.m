%% User Input
DataDir = '..\Subjects'
MiceOpt = {'EB002'}%,'EB001','EB003'}%{'EB001','EB002','EB003','CB007'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'
SaveDir = '..\Data\Preproccesed\'
tmpdatafolder = '..\tmpdata\';
LocalDir = '..\Data\KiloSortOutput' 
IBLDir = '..\Data\IBLPrep' 
AllenCCFPath = 'C:\Users\EnnyB\Documents\MATLAB\allenCCF'
storevideopath='..\Data\Videos'
HistoFolder = '..\Data\Histology\';
nidq_sync_used = zeros(1,length(MiceOpt));
nidq_sync_used(ismember(MiceOpt,{'EB001','CB007','CB008'}))=1;
AREASOfInterest = {'RSP','CA1','DG','LGd','TH','VISp','ZI','LP','MG','SGN','SUB','or','PIL','Pro','CB','MB'}

maxsessnr = 2; %max nr. sessions on a day (doesn't need to be accurate)
MinDist2Include = 50; %Including only trials in which mouse reached at least 70cm for some statistics

% addpath(genpath(cd))
% addpath(genpath('C:\Users\EnnyB\Documents\GitHub\spikes'))
% addpath(genpath('C:\Users\EnnyB\Documents\GitHub\npy-matlab'))
% addpath(genpath('C:\Users\EnnyB\Documents\GitHub\AP_histology'))


%% VR - basic behavioral plots
[AllContr,AllSound,AllPosition,AllTime,AllTraj,AllRewarded,AllRewardPos,SessionOpt] = MainBehaviorVR(DataDir,MiceOpt,maxsessnr,MinDist2Include)

%% Videos to process with DeepLabCut
Listvideofiles(DataDir,MiceOpt,storevideopath,0)

%% Kilosort?
if 0
    addpath(cd,'Neuraldata')
    main_kilosort_Local
    
    disp('You may want to do histology alignment in python now?')
end

%% Timeline
[AllTimeLines,AllInputs] = ExtractTimeline(DataDir,SessionOpt,MiceOpt)

%% Neural Data - Needs to be preprocessed by Kilosort / Phy
LoadSpikes

%% Spatial Coding X Sensory Processing
SpatialCodingAnalysis


%% Microphone - extract data
[AllMicrophone,AllFS,Allnbits] = MicrophoneExtract(DataDir,MiceOpt)