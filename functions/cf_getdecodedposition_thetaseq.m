function decodedData = cf_getdecodedposition_thetaseq(trainingData, testingData,...
    nCells, inclCells, params)
%cf_getdecodedposition_thetaseq
%
%ALP 1/5/2023

timeAroundTrough = params.timeAroundTrough;
timeEdges = params.timeEdges; 
decodingBinsize = params.decodingBinsize;
posEdges = params.posEdges;
posBinsize = params.posBins; 
adjustedEdges = params.adjustedEdges;

qX = params.qX;
qY = params.qY; 

if max(posEdges) > 180
    nDeg = 360;
    dec_edges = 0:18:360;%edges for calculating prospective quadrant ratio over positio
else 
    nDeg = 180;
    dec_edges = -81:3:99;
end
posEdges = posEdges(1:end-1);

%%% convert trainingData into trainingcounts
tmpTrainingCounts = trainingData.ratemaps.*decodingBinsize;  %expected counts in the bin with size decodingBinsize in s
for u = 1:nCells
    cellI = inclCells(u);
    trainingCounts(u,:) = tmpTrainingCounts(cellI,:); 
    if sum(tmpTrainingCounts(cellI,:) == 0) > 0
        isZero = tmpTrainingCounts(cellI,:) == 0;
        trainingCounts(u,isZero) = 1e-15; 
    end
end

if isempty(testingData)
    decodedData = []; 
    return
end

%%% get theta troughs in each trial
for t = 1:length(testingData)
    for w = 1:size(testingData(t).counts,3)
        tmpTestingCounts = testingData(t).counts(:,:,w);
        decodeOut = cf_getspatialprob(trainingCounts, tmpTestingCounts);
        nSpikes = decodeOut.totalSpikes;
        excludeBins = nSpikes < 2; %only include bins with 2 or more spikes
        
        tmpEstPos = decodeOut.spatialProb;
        %check this nan addition is working OK
        tmpEstPos(:, excludeBins) = NaN(size(tmpEstPos,1), sum(excludeBins));
        
        % transform into relative position
        tmpTroughPos = testingData(t).troughPos(w); 
        posdiff = posEdges - tmpTroughPos; 
        posdiff(posdiff<-(nDeg/2)) = posdiff(posdiff <-(nDeg/2)) + nDeg;
        posdiff(posdiff>(nDeg/2)) = posdiff(posdiff > (nDeg/2)) - nDeg;
        [~, iorder] = sort(posdiff); 
        
        decodedSeq(:,:,w) = tmpEstPos(iorder,:); 
    end
    
    %%% trial average
    mnTrialSeq = mean(decodedSeq,3, 'omitnan');
    trialPRatio = cf_calcquadrantratio(qX, qY, mnTrialSeq, 'mod');
    
    %%% get trial prospective coding ratio across space
    [~,~,binI] = histcounts(testingData(t).troughPos, dec_edges);
    for b = 1:length(dec_edges)-1
        isW = binI == b;
        tmpSeq = mean(decodedSeq(:,:,isW),3, 'omitnan');
        trialPRatio_pos(b) = cf_calcquadrantratio(qX, qY, tmpSeq, 'mod');
    end
%     
%     tmpDecodedTrial = cf_getspatialprob(trainingCounts, testingData(t).fullTrialCounts);
%     % adjust position
%     for b = 1:size(testingData(t).fullTrialCounts,2)
%         tmpTroughPos = testingData(t).fullTrialPos(b);
%         posdiff = posEdges - tmpTroughPos; 
%         posdiff(posdiff < -(nDeg/2)) = posdiff(posdiff < - (nDeg/2)) + nDeg;
%         posdiff(posdiff > (nDeg/2)) = posdiff(posdiff > (nDeg/2)) - nDeg;
%         [~, iorder] = sort(posdiff); 
%         trialDecoding(:,b) = tmpDecodedTrial.spatialProb(iorder,b); 
%     end
    
    % average theta seq across the whole trial
%     figure
%     hold on
%     imagesc(timeEdges(1:end-1), adjustedEdges(1:end-1)+1, mean(decodedSeq,3,'omitnan'))
%     plot(timeEdges(1:end-1), zeros(1,size(decodedSeq,2)), 'w--')
%     plot(zeros(1,size(decodedSeq,1)), adjustedEdges(1:end-1), 'w--')
%     
%     figure
%     hold on
%     imagesc(1:size(testingData(t).fullTrialCounts,2), adjustedEdges, trialDecoding, [0 0.08])
%     
%     figure
%     hold on
%     imagesc(tmpDecodedTrial.spatialProb, [0 0.05])
%     plot(testingData(t).fullTrialPos./2+1, 'r-')

    decodedTrial(t).trialSeq = decodedSeq;
    decodedTrial(t).trialNum = t.*ones(1,size(decodedSeq,3));
    decodedTrial(t).troughPos = testingData(t).troughPos; 
    decodedTrial(t).mnSeq = mnTrialSeq;
    decodedTrial(t).mnSeq_PRatio = trialPRatio; 
    decodedTrial(t).trial_PRatio_pos = trialPRatio_pos;
    clear decodedSeq tmp* 
end

%%% average theta seq across the day
allDecodedSeq = {decodedTrial.trialSeq};
allDecodedSeq = cat(3, allDecodedSeq{:});
mnDecodedSeq = mean(allDecodedSeq,3, 'omitnan');
troughPos = {decodedTrial.troughPos};
troughPos = cell2mat(troughPos');
trialID = [decodedTrial.trialNum];

%%% get overall quadrant value
fullQRatio = cf_calcquadrantratio(qX, qY, mnDecodedSeq, 'full');

%%% shuffle to get significance
for s = 1:500
    rVals = [];
    rVals = arrayfun(@(x) randperm(size(allDecodedSeq,2))', 1:size(allDecodedSeq,3), 'UniformOutput', false);
    rVals = cell2mat(rVals)';
    
    shuffSeq = arrayfun(@(x) allDecodedSeq(:, rVals(x,:), x), 1:size(allDecodedSeq,3), 'UniformOutput', false);
    shuffSeq = cat(3, shuffSeq{:});
    shuffMean = nanmean(shuffSeq,3);
    
    shuffQRatio(s) = cf_calcquadrantratio(qX, qY, shuffMean, 'full');
end

%%% get average across position
[~,~,binI] = histcounts(troughPos, dec_edges);
for b = 1:length(dec_edges)-1
    tmpSeq = [];
    isW = binI == b;
    tmpSeq = mean(allDecodedSeq(:,:,isW),3, 'omitnan');
    pRatio_pos(b) = cf_calcquadrantratio(qX, qY, tmpSeq, 'mod');
end

%%% add stuff to output structure
decodedData.allDecodedSeq = allDecodedSeq;
decodedData.troughPos = troughPos; 
decodedData.trialID = trialID; 
% decodedData.trialData = decodedTrial; 
decodedData.mnDecodedSeq = mnDecodedSeq;
decodedData.mnQuadrantRatio = fullQRatio;
decodedData.shuff_QR_thresh = prctile(shuffQRatio, 95);
decodedData.significantSeq = fullQRatio > prctile(shuffQRatio,95);

figure('Position', [141 547 1227 251])
hold on
subplot(1,4,1)
hold on
imagesc(timeEdges(1:end-1), adjustedEdges(1:end-1), mean(allDecodedSeq,3,'omitnan'))
plot(timeEdges(1:end-1), zeros(1,size(mnDecodedSeq,2)), 'w--')
plot(zeros(1,size(mnDecodedSeq,1)), adjustedEdges(1:end-1), 'w--')
xlim([min(timeEdges) max(timeEdges)])
ylim([min(adjustedEdges), max(adjustedEdges)])
xlabel('time from theta trough')
ylabel('relative decoded position (deg)')
title('total sequence')
colorbar
subplot(1,4,2)
hold on
imagesc(timeEdges(min(min(qX)):max(max(qX))), adjustedEdges(min(min(qY)):max(max(qY))), mnDecodedSeq(min(min(qY)):max(max(qY)), min(min(qX)):max(max(qX))))
plot(timeEdges(1:end-1), zeros(1,size(mnDecodedSeq,2)), 'w--')
plot(zeros(1,size(mnDecodedSeq,1)), adjustedEdges(1:end-1), 'w--')
xlim([timeEdges(min(min(qX)))-decodingBinsize/2, timeEdges(max(max(qX)))+decodingBinsize/2])
ylim([adjustedEdges(min(min(qY))), adjustedEdges(max(max(qY)))])
xlabel('time from theta trough')
ylabel('relative decoded position (deg)')
title('for quadrant')
colorbar
subplot(1,4,3)
hold on
plot(dec_edges(1:end-1), pRatio_pos)
xlim([min(dec_edges) max(dec_edges(1:end-1))])
xlabel('position (deg)')
ylabel('prospective coding ratio')
title('propsective coding ratio over position')
subplot(1,4,4)
hold on
plot(prctile(shuffQRatio, 95), 'k*')
plot(fullQRatio, 'r*')
title('true quad vs. shuffled')


end

