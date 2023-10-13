function out = cf_gettestingdata_ripples(spiketimes, spikeIDs, ripplewindows, binsize, cellIDs, nCells)
%cf_gettestingdata_ripples
%
%ALP 1/11/2023

for i = 1:size(ripplewindows,1)
    rBins(i,:) = ripplewindows(i,1):binsize:ripplewindows(i,2);
    isActiveCell = isExcluded(spiketimes, [ripplewindows(i,1) ripplewindows(i,2)]);
    isActiveCell = logical(isActiveCell);
    activeIDs{i} = unique(spikeIDs(isActiveCell));
    
    for c = 1:nCells
        iC = cellIDs(c); 
        tmpSpikes = spiketimes(spikeIDs == iC);
        testingCounts{i}(iC,:) = histcounts(tmpSpikes, rBins(i,:));
    end
end

out.counts = testingCounts;
out.activeIDs = activeIDs;
out.windows = rBins;
out.binsize = binsize;


end

