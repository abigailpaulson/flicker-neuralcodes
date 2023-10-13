function cf_flicker_timesequences(dirs, params, allindex, metadata, pos, cellTypes)
%cf_flicker_timesequences
%   playing around with calculating time sequences during flicker
%ALP 2/28/23

%%% set up important parameters
frbins = 1; %in s
nblocks = 10;

%%% index, directories, etc
index = [28 200823 2];
dayindex = index(:,1:2);
d = 1;
anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
anspikedir = anprocesseddatadir;

%%% load spiking activity
[spikes, unitInfo, spikevect] = cf_getspikes(anspikedir, dirs, params, index, params.brainReg);

%%% load stim periods
load([anprocesseddatadir, 'stimperiods', num2str(index(1,3)), '.mat'])
stimTime = stimperiods{index(1,1)}{index(1,2)}{index(1,3)}.visualper; 
stimTime = stimTime./stimperiods{index(1,1)}{index(1,2)}{index(1,3)}.samprate; 

%%% get cell types of interest
ct = {unitInfo.cellType};
inclCell = strcmp(cellTypes, ct); 
goodSpikes = spikes(index(1,3)).data(inclCell);

%%% set up gaussian window for FR calculation
% smooth with gaussian window, Wang et al. 2015 use 75ms sd gaussian for
% smoothing. Wang et al also use 1250hz samprate for spike times. maybe I
% could try 2000hz since that is an easy downsample
% winLength = 1000; %in ms
% stdDev = 75; %in ms
% alpha = (winLength-1)/(stdDev*2); 
% gwin = gausswin(winLength, alpha); 
% gwin = gwin(

%%% get firing rate over the stim time
%make spike trains
trainedges = stimTime(1):frbins:stimTime(2); %in seconds
spiketrains = NaN(length(goodSpikes), length(trainedges)-1);
firingrate = NaN(length(goodSpikes), length(trainedges)-1);
for u = 1:length(goodSpikes)
    tmpTimes = goodSpikes(u).spikeTimes; 
    tmpTrain = histcounts(tmpTimes, trainedges); 
%     tmpFR = conv(gwin, tmpTrain, 'same'); 
    
    spiketrains(u,:) = tmpTrain; 
    firingrate(u,:) = tmpTrain./frbins; 
    
    clear tmp
end

%%% break into blocks of time
%%% try 10 blocks? 
blockinds = round((length(trainedges)-1)/nblocks);
blockedges = 1:blockinds:length(trainedges)-1;
plottime = 0:1:blockinds;
plottime = plottime*frbins; 
plottime  = plottime/60;

for b = 1:length(blockedges)-1
    tmpblock = firingrate(:,blockedges(b):blockedges(b+1)-1);
    [peakVal, iMax] = max(tmpblock,[],2, 'omitnan');
    [~, iSort] = sort(iMax);
    tmblock = tmpblock./peakVal;
    
    block_fr(:,:,b) = tmpblock(iSort,:);
end

%%% plot
for b = 1:length(blockedges)-1
    figure
    hold on
    imagesc(plottime(1:end-1), 1:length(goodSpikes), block_fr(:,:,b))
end




end

