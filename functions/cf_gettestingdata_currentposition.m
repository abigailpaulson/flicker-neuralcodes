function testingData = cf_gettestingdata_currentposition(trialData, cellIDs, nCells, speedThreshold, decodingBins, posName)
%cf_gettestingdata_currentposition
%
%ALP 1/4/2023

for t = 1:length(trialData)
    %%% split each trial into time bins
    tmpTime = trialData(t).time;
    tmpTimeEdges = tmpTime(1):decodingBins:tmpTime(end);
    if tmpTimeEdges(end) > tmpTime(end)
        tmpTimeEdges = tmpTimeEdges(1:end-1);
    end
    [~, ~, bins] = histcounts(tmpTime, tmpTimeEdges);
    
    %%% only include moving bins
    tmpSpeed = arrayfun(@(x) nanmean(trialData(t).speed(bins == x)), 1:length(tmpTimeEdges(1:end-1))); 
    inclBins = tmpSpeed > speedThreshold; 
    
    %%% get actual bin position for each time bin
    tmpActualPos = arrayfun(@(x) nanmean(trialData(t).(posName)(bins == x)), 1:length(tmpTimeEdges(1:end-1)));
    tmpActualPos = tmpActualPos(inclBins); %moving bins
    
    %%% get spike histograms over time for each trial
    for u = 1:nCells
        iC = cellIDs(u);
        isCell = trialData(t).spikeIDs == iC; 
        spkTimes = trialData(t).spikeTimes(isCell);
        tmpCellH = histcounts(spkTimes, tmpTimeEdges);
        tmpTestingCounts(u,:) = tmpCellH(inclBins);
        
        clear spkTimes
    end
    
    testingData(t).counts = tmpTestingCounts; % counts for each cell in each time bin [cells x timeBin]
    testingData(t).timeEdges = tmpTimeEdges; %only times during the trial when moving
    testingData(t).realPos = tmpActualPos; %average position in the time bin
    
    clear tmp*
end


end

