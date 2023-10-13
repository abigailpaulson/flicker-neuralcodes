function out = cf_getsignificant_thetaseq_table(decodedData, params, metrics)
%cf_getsignificat_thetaseq_table
%   shuffle to find significance
%
%ALP 3/6/2023

qX = params.qX; 
qY = params.qY; 
% qY = [31 60; 31 60; 1 30; 1 30]; %ALP 3/623 trying this just to see how
% it changes things if included all position

inclTrials = metrics.fullTrial & metrics.rewarded & metrics.engaged;
allDecodedSeq = decodedData.trialSeq; 
allDecodedSeq = cat(3, allDecodedSeq{inclTrials}); 
mnSeq = mean(allDecodedSeq,3, 'omitnan');
day = metrics.day(1);

if isempty(mnSeq)
    dayQRatio = NaN; shuff_thresh = NaN; significantSeq = 0;
    out = table(dayQRatio, shuff_thresh, significantSeq, day); 
    return
end

dayQRatio = cf_calcquadrantratio(qX, qY, mnSeq, 'full');



for s = 1:500
    rVals = [];
    rVals = arrayfun(@(x) randperm(size(allDecodedSeq,2))', 1:size(allDecodedSeq,3), 'UniformOutput', false);
    rVals = cell2mat(rVals)';
    
    shuffSeq = arrayfun(@(x) allDecodedSeq(:, rVals(x,:), x), 1:size(allDecodedSeq,3), 'UniformOutput', false);
    shuffSeq = cat(3, shuffSeq{:});
    shuffMean = nanmean(shuffSeq,3);
    
    shuffQRatio(s) = cf_calcquadrantratio(qX, qY, shuffMean, 'full');
    
    clear shuffSeq shuffMean rVals
end


shuff_thresh = prctile(shuffQRatio, 95);
significantSeq = dayQRatio > prctile(shuffQRatio,95);

out = table(dayQRatio, shuff_thresh, significantSeq, day); 


end

