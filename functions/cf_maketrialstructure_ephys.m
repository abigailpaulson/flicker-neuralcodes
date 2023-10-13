function cf_maketrialstructure_ephys(dirs, params, allindex, metadata, btype)
%cf_maketrialstructure_ephys
%   add spiking info to the behavior trial structure
%ALP 12/23/2022

dayindex = unique(allindex(:,1:2), 'rows');
positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';

for d = 1:size(dayindex,1)
    index = allindex(allindex(:,2) == dayindex(d,2),:);
    
    %%% load trial data structure
    trialfilename = ['trialInfo_', btype, '_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    trialSpikesfilename = ['trialSpikes_', btype, '_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    load([positiondir, trialfilename])
    
    %%% load spiking for this day 
    anspikedir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
    [spikes, unitInfo, spikevect] = cf_getspikes(anspikedir, dirs, params, index, params.brainReg);

    %%% divide spikes into trials 
    for t = 1:length(trialData)
        iF = trialData(t).file;
        fileSpikes = spikevect(iF).spikeTimes;
        fileIDs = spikevect(iF).spikeIDs;
        trialTime = [trialData(t).starttime trialData(t).endtime];
        
        isTSpike = isExcluded(fileSpikes, trialTime);
        isTSpike = logical(isTSpike); 
        
        %%% get spike position inds
        tTime = trialData(t).time;
        spikePosInds = lookup2(fileSpikes(isTSpike), tTime);
        
        %%% save to output structure
        trialSpikes.data(t).spikeTimes = fileSpikes(isTSpike);
        trialSpikes.data(t).IDs = fileIDs(isTSpike); 
        trialSpikes.data(t).spikePosInds = spikePosInds;
        trialSpikes.data(t).file = iF;
        trialSpikes.info.clusterID = [unitInfo.ID];
        trialSpikes.info.brainReg = {unitInfo.brainReg};
        trialSpikes.info.cellType = {unitInfo.cellType};
        trialSpikes.info.index = index;
        
        clear fileSpikes trialTime iF
    end
    
    info = []; 
    info = addhelpfulinfotostruct(info); 
    save([positiondir, trialSpikesfilename], 'trialSpikes', 'info')
    
    clear trial* file* spikes unitInfo spikevect
end


end

