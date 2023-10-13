function out = cf_getdecodedposition_ripples(trainingData, testingData)
%cf_getdecodedposition_ripples
%
%ALP 1/11/2023


for i = 1:size(testingData.windows,1)
    eventCells = testingData.activeIDs{i};
    
    eventTrainingCounts = trainingData.ratemaps(eventCells,:).*testingData.binsize;
    for c = 1:size(eventTrainingCounts,1)
        if sum(eventTrainingCounts(c,:)) > 0
            isZero = eventTrainingCounts(c,:) == 0;
            eventTrainingCounts(c,isZero) = 1e-15;
        end
    end
    
    eventTestingCounts = testingData.counts{i}(eventCells,:);
    
    if isempty(eventTrainingCounts) || isempty(eventTestingCounts)
        clear event*
        continue
    end
    
    decodedEvent = cf_getspatialprob(eventTrainingCounts, eventTestingCounts);
    
    isDecodeBin = decodedEvent.totalSpikes > 0; %need at least 1 spike to decode
    
    decodedEvent.spatialProb(:,~isDecodeBin) = NaN(size(decodedEvent.spatialProb,1), sum(~isDecodeBin));
   
    allSpatialProb{i} = decodedEvent.spatialProb;
    allTrainingCounts{i} = eventTrainingCounts;
    allTestingCounts{i} = eventTestingCounts;
    allDecodeBins{i} = isDecodeBin;
    
    clear event* decodedEvent
end

out.spatialProb = allSpatialProb;
out.trainingCounts = allTrainingCounts;
out.allTestingCounts = allTestingCounts;
out.isDecodeBin = allDecodeBins; 


end

