function out = cf_getspatialprob(trainingCounts, testingCounts)
%getspatialprob
% Karlsson and Frank 2009
%       trainingCounts [nCells x SpatialBins]
%       testingCounts [nCells x TimeBins]
%
%   out - spatialProb [SpatialBins x TimeBins]
%ALP 6/11/2020
%edited ALP 1/4/2023

isLowOcc = isnan(trainingCounts(1,:));
if sum(isLowOcc) > 1
    A = 1;
end
for iTime = 1:size(testingCounts,2)
    binSpikeCount = repmat(testingCounts(:,iTime), 1, size(trainingCounts,2));
    
    prob = prod(((trainingCounts.^binSpikeCount)./factorial(binSpikeCount)).*exp(-trainingCounts),1)';
    
    if sum(isnan(prob)) > 0
        A = 1;
    end
    
    %if trainingCounts(:,b) is NaN, that means the occupancy was too low
    %in bin b. This should be the same across all cells. make the
    %spatial Prob 0 in those cases
    prob(isLowOcc,:) = 0;
    
    normProb = prob/sum(prob); % normalize across space to make the probabilities add up to 1
    

    spatialProb(:,iTime) = normProb;
end

if sum(isnan(spatialProb)) > 0
    A = 1; 
end

out.spatialProb = spatialProb; %[bins x Time]
out.totalSpikes = sum(testingCounts,1);
end

