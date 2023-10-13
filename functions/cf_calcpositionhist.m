function out = cf_calcpositionhist(cellIDs, spiketimes, spikeIDs, vrtime, vrpos, vrspeed, edges, speedthresh, shuffle)
%cf_calcpositionhist
%   for rate maps
%   random shuffle
%ALP 12/19/2022

%%% throw away any early spikes
starttime = vrtime(1);
endtime = vrtime(end);
isVRspike = isExcluded(spiketimes, [starttime endtime]);
spiketimes = spiketimes(logical(isVRspike));
spikeIDs = spikeIDs(logical(isVRspike));
vrduration = endtime-starttime; 

nSpikes = length(spiketimes);
pos_count = []; 

if shuffle
    shiftVect = 20 + (vrduration-20).*rand(nSpikes,1); %shift each spike by a random amount
    spiketimes_shuffled = shiftVect+spiketimes;
    toWrap = spiketimes_shuffled > max(vrtime);
    spiketimes_shuffled(toWrap) = spiketimes_shuffled(toWrap) - vrduration;
    spiketimes = spiketimes_shuffled;
end

isMoving = vrspeed > speedthresh;
movingTimes = vrtime(isMoving);

% spiketimes = spiketimes(spiketimes >= starttime);
% spikeIDs = spikeIDs(spiketimes >= starttime);

%%% get spike positions
spikeVRInd = lookup2(spiketimes, vrtime); %index for all spikes in VR
spikeVRTimes = vrtime(spikeVRInd);
isMovingSpike = ismember(spikeVRTimes, movingTimes);

allspikepositions = vrpos(spikeVRInd);
inclSpikes = isMovingSpike;

%%% only calc position count on moving spikes
use_spikepositions = allspikepositions(inclSpikes);
use_spikeIDs = spikeIDs(inclSpikes);

for u = 1:length(cellIDs)
    unitID = cellIDs(u);
    isCellSpk = use_spikeIDs == unitID;
    tmp_pos = use_spikepositions(isCellSpk);
    pos_count(unitID,:) = histcounts(tmp_pos, edges);
end

out = pos_count;

end

