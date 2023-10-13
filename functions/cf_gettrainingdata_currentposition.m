function out = cf_gettrainingdata_currentposition(trialData, cellIDs, nCells, speedThreshold, edges, posName)
%cf_gettrainingdata_currentposition
%
%   out - training data rate maps [nUnits x position] in Hz
%ALP 12/28/2022

vrSteps = 0.02; %0.02sec in each timebin
vrSamprate = 50; %1/0.02 in hz

%%% initialize output structure
trainingCounts = zeros(nCells, length(edges)-1); 
trainingOcc = zeros(1,length(edges)-1);

for t = 1:length(trialData)
    isMoving = trialData(t).speed > speedThreshold;
    movingInds = find(isMoving == 1); 
    
    %%% occupancy
    tmpPos = trialData(t).(posName); 
    movingPos = tmpPos(isMoving); 
    movingCounts = histcounts(movingPos, edges); 
    movingTime = movingCounts/vrSamprate; %bins/(bins/s) = time in each bin
    trainingOcc = trainingOcc+movingTime; 
    
    %%% spikes over position
    isMovingSpike = ismember(trialData(t).spikePosInds, movingInds);
    tmpSpikes = trialData(t).spikeTimes(isMovingSpike);
    tmpSpikePosInds = trialData(t).spikePosInds(isMovingSpike); 
    tmpSpikePos = tmpPos(tmpSpikePosInds); 
    tmpSpikeIDs = trialData(t).spikeIDs(isMovingSpike);

    for u = 1:nCells
        iC = cellIDs(u);
        isCell = tmpSpikeIDs == iC; 
        tmpCellPos = tmpSpikePos(isCell); 
        tmpPosH = histcounts(tmpCellPos, edges);
        
        trainingCounts(u,:) = trainingCounts(u,:) + tmpPosH; 
        
        clear tmpPosH
    end
    
    clear tmp*
end

%%% smooth the counts
for u = 1:nCells
    smoothCounts(u,:) = gaussSmooth(trainingCounts(u,:), 2);
end
    
%%% smooth the Occ
smoothOcc = gaussSmooth(trainingOcc, 2);

%%% training firing rate maps
trainingMaps = smoothCounts./smoothOcc;

%%% output structure
out.ratemaps = trainingMaps;
out.occ = smoothOcc; 
out.counts = smoothCounts;

%%% plot to check
% figure
% hold on
% title([num2str(length(trialData)), ' trials - ', num2str(sum(trainingOcc)), 's total occ time'])
% for u = 1:5
%     plot(edges(1:end-1), trainingMaps(u,:))
% end
% ylabel( 'firing rate(hz)')
% xlabel('position (Deg)')


end

