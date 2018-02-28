function kilosort_cluster_quality(data_root)
% for UCLA 64 Chan silocon probe 
% Modified from plotNeuronFigsScript.m
% Lu, June 21 2017 


if nargin<1
    data_root = '.\';
end

INCLUDE_UNIT_TYPE = 2;

%%
ksRoot = [data_root 'Kilosort'];

loadPars.loadPCs = true;
sp = loadKSdir(ksRoot, loadPars);

inclCID = sp.cids(sp.cgs==INCLUDE_UNIT_TYPE);
st = sp.st;
clu = sp.clu;

pcFeat = sp.pcFeat;
pcFeatInd = sp.pcFeatInd;
spikeAmps = sp.tempScalingAmps;

figDir = fullfile(data_root, 'figs'); 
if ~exist(figDir,'dir'); mkdir(figDir); end

params.dataType = sp.dtype;
params.filename = fullfile(ksRoot, sp.dat_path) ;
d = dir(params.filename); 
nSamp = d.bytes/2/sp.n_channels_dat;

params.dataSize = [sp.n_channels_dat nSamp];
params.chanMap = readNPY(fullfile(ksRoot, 'channel_map.npy'));
params.Fs = sp.sample_rate;
params.xcoords = sp.xcoords; 
params.ycoords = sp.ycoords;
params.plotDistance = 100;
params.nWFsToLoad = 1000;
params.nWFsToPlot = 100;
params.gain = 0.6/512/500*1e6; % raw file units to uV
params.nPCsToPlot = 50000;
% params.highlightRegions = [132 1104];

%% extract median WFs (just once)
medwf_file = fullfile(ksRoot, 'medWFs.mat');
if exist(medwf_file, 'file')
    load(medwf_file)
    
else
    inclSP = ismember(clu, sp.cids(sp.cgs==INCLUDE_UNIT_TYPE));
    medWFs = extractMedianWFs(clu(inclSP), st(inclSP), params.Fs, params.filename, ...
        params.dataType, params.dataSize, params.chanMap, params.gain);
    
    save(medwf_file, 'medWFs');
end

%% compute cluster quality stats (just once)
cluster_quality_file = fullfile(ksRoot, 'clusterQualityMetrics.mat');
if exist(cluster_quality_file, 'file')
    load(cluster_quality_file)
    
else
    [cgs, uQ, cR, isiV] = sqKilosort.computeAllMeasures(ksRoot);
    save(cluster_quality_file, 'cgs', 'uQ', 'cR', 'isiV');
    
end
%%

sparsePCfeat = sparsePCs(pcFeat, pcFeatInd, sp.spikeTemplates, 2, 10);

%%

theseISI = isiV(cgs==INCLUDE_UNIT_TYPE);
theseCR = cR(cgs==INCLUDE_UNIT_TYPE);
theseID = uQ(cgs==INCLUDE_UNIT_TYPE);

for q = 1:length(inclCID)
    clusterID = inclCID(q);

    stats.medWF = squeeze(medWFs(inclCID==clusterID,:,:))';
    stats.isiContamination = theseISI(inclCID==clusterID);
    stats.isoDistance = theseID(inclCID==clusterID);
    stats.mahalContamination = theseCR(inclCID==clusterID);
%     stats.isiContamination = 0;
%     stats.isoDistance = 0;
    figHand = neuronFig(clusterID, st, clu, sparsePCfeat, spikeAmps, stats, params);
    set(figHand, 'Position', [-1890         -59        1810        1031]);
    saveas(figHand, fullfile(figDir, sprintf('/cluster%03d.png', clusterID)))
    close(figHand); clear figHand
    
end

%%

% [~, spikeDepths, ~, ~, ~, ~, ~] = ...
%     templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
% cluDepths = clusterAverage(sp.clu(ismember(sp.clu,sp.cids(sp.cgs==INCLUDE_UNIT_TYPE))), spikeDepths(ismember(sp.clu,sp.cids(sp.cgs==INCLUDE_UNIT_TYPE))));
% makeClusterFigsWebsite(figDir, 'sortByDepth.html', cluDepths);
% 
% cluAmps = max(max(medWFs,[],3)-min(medWFs,[],3), [], 2);
% makeClusterFigsWebsite(figDir, 'sortByAmp.html', cluAmps);
% 
% [~, FRs] = countUnique(sp.clu(ismember(sp.clu, sp.cids(sp.cgs==INCLUDE_UNIT_TYPE))));
% makeClusterFigsWebsite(figDir, 'sortByFR.html', FRs);
% 
% makeClusterFigsWebsite(figDir, 'sortByIsoDist.html', uQ(cgs==INCLUDE_UNIT_TYPE));
