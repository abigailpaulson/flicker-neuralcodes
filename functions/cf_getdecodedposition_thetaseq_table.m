function Out = cf_getdecodedposition_thetaseq_table(trainingData, testingData,...
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
for t = 1:length(testingData) %decode all trials
    decodedSeq = NaN( length(posEdges), length(timeEdges)-1, size(testingData(t).counts,3)); 
    for w = 1:size(testingData(t).counts,3)
        tmpTestingCounts = testingData(t).counts(:,:,w);
        if isempty(tmpTestingCounts)
            continue
        end
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
    trialQRatio = cf_calcquadrantratio(qX, qY, mnTrialSeq, 'full');
    
    %%% get trial prospective coding ratio across space
    [~,~,binI] = histcounts(testingData(t).troughPos, dec_edges);
    for b = 1:length(dec_edges)-1
        isW = binI == b;
        tmpSeq = mean(decodedSeq(:,:,isW),3, 'omitnan');
        trialPRatio_pos(b) = cf_calcquadrantratio(qX, qY, tmpSeq, 'mod');
        trialQRatio_pos(b) = cf_calcquadrantratio(qX, qY, tmpSeq, 'full');
    end
    
    trialtroughpos = testingData(t).troughPos;
    if isempty(trialtroughpos)
        trialtroughpos = NaN;
    end
    
    decodedTrial(t).trialSeq = decodedSeq;
    decodedTrial(t).trialID = t.*ones(size(decodedSeq,3),1);
    decodedTrial(t).troughPos = trialtroughpos;
    decodedTrial(t).mnSeq = mnTrialSeq;
    decodedTrial(t).mnSeq_QRatio = trialQRatio; 
    decodedTrial(t).mnSeq_PRatio = trialPRatio;
    decodedTrial(t).trial_PRatio_pos = trialPRatio_pos;
    decodedTrial(t).trial_QRatio_pos = trialQRatio_pos; 
    clear decodedSeq tmp* trial*
end
% allDecodedSeq = {decodedTrial.trialSeq};
% allDecodedSeq = cat(3, allDecodedSeq{:});
% troughPos = {decodedTrial.troughPos};
% troughPos = cell2mat(troughPos');
% trialID = [decodedTrial.trialNum]';


%%% add stuff to output structure
Out = struct2table(decodedTrial); 



end

