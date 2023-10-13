function out = cf_getdecodedposition_currentposition(trainingData, testingData, nCells, decodingBinsize, posEdges, posName)
%cf_getdecodedposition_currentposition
%
%ALP 1/4/2023

% allBins = length(posEdges)-1;
% halfBin = allBins/2+1;
% quarterBin = allBins/4;
% forwardInds = halfBin+1:1:halfBin+quarterBin;
% behindInds = halfBin-quarterBin:halfBin-1; 

trainingCounts = trainingData.ratemaps.*decodingBinsize;  %expected counts in the bin with size decodingBinsize in s
for u = 1:nCells
    if sum(trainingCounts(u,:) == 0) > 0
        isZero = trainingCounts(u,:) == 0;
        trainingCounts(u,isZero) = 1e-15; 
    end
end

%%% for each trial, decode each time bin
for t = 1:length(testingData)
    % get estimated position
    tmpTestingCounts = testingData(t).counts;
    decodeOut = cf_getspatialprob(trainingCounts, tmpTestingCounts);
    nSpikes = decodeOut.totalSpikes;
    excludeBins = nSpikes < 2; %only include bins with 2 or more spikes
    
    tmpEstPos = decodeOut.spatialProb(:,~excludeBins); 
    
    % get decoding error for each time bin
    [~, iMax] = max(tmpEstPos, [], 1);
    maxEstPos = posEdges(iMax);
    tmpError = maxEstPos - testingData(t).realPos(~excludeBins);
    
    decodedData(t).estPos = tmpEstPos;
    decodedData(t).maxEstPos = maxEstPos;
    decodedData(t).realPos = testingData(t).realPos(~excludeBins);
    decodedData(t).error = tmpError;
%     decodedData(t).ratio = (sum(tmpEstPos(forwardInds,:),1) - sum(tmpEstPos(behindInds,:),1))./(sum(tmpEstPos(forwardInds,:),1) + sum(tmpEstPos(behindInds,:),1));
    
    % plot to check - looks good! ALP 1/4/2023
%         figure
%         hold on
%         subplot(2,1,1)
%         hold on
%         imagesc(tmpEstPos)
%         subplot(2,1,2)
%         hold on
%         plot(1:size(tmpEstPos,2), testingData(t).realPos, 'r-')
%         plot(1:size(tmpEstPos,2), maxEstPos, 'g-')
%         subplot(3,1,3)
%         hold on
%         plot(1:size(tmpEstPos,2), decodedData(t).ratio)
    clear tmp*
end

%%% get decoding estimate average across positions
allEstPos = [decodedData.estPos];
allRealPos = [decodedData.realPos];
allMaxEstPos = [decodedData.maxEstPos];
allErr = [decodedData.error];
if max(posEdges) > 180
    allErr(allErr > 180) = allErr(allErr > 180) - 360;
    allErr(allErr < -180) = 360 + allErr(allErr < -180);
else
    allErr(allErr > 90) = allErr(allErr > 90) - 180; 
    allErr(allErr < -90) = allErr(allErr < -90) + 180; 
end

%bin real positions
[~, ~, bin] = histcounts(allRealPos, posEdges);
mnEstPos = NaN(length(posEdges)-1, length(posEdges)-1); 

for b = 1:length(posEdges)-1
    isBin = bin == b;
    tmpBinEst = mean(allEstPos(:,isBin), 2, 'omitnan'); 
    mnEstPos(:,b) = tmpBinEst; 
    
    binError(b) = mean(allErr(isBin), 'omitnan');
    
    clear tmpBinEst
end

% figure
% hold on
% subplot(2,1,1)
% hold on
% imagesc(posEdges, posEdges, mnEstPos, [0 0.2])
% % xlim([0 360])
% % ylim([0 360])
% xlabel('actual position (deg)')
% ylabel('estimated position (deg)')
% subplot(2,1,2)
% hold on
% plot(posEdges(1:end-1), binError)

%%% create output structure
out.mnEstPos = mnEstPos; 
out.binError = binError; 
out.decodedTrials = decodedData;

end

