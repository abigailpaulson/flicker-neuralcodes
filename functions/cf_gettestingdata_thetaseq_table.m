function [testingData, mnThetaTrace] = cf_gettestingdata_thetaseq_table(index, trialData, theta, ...
    decodingBinsize, timeAroundTrough, speedThreshold, posName, nCells, inclCells)
%cf_gettestingdata_thetaseq
%   adapted from script_thetasequences.m in abby\code\chronicflicker-ephys\
% ALP 1/5/23

testingData = []; allThetaTrace = [];
for i = 1:size(index,1)
    %%% theta properties
    phase{index(i,3)} = theta{index(i,1)}{index(i,2)}{index(i,3)}.data(:,2); %in our preproc pipeline, 0 is the peak, -pi and pi is the trough
    amp{index(i,3)} = theta{index(i,1)}{index(i,2)}{index(i,3)}.data(:,1);
    samprate = theta{index(i,1)}{index(i,2)}{index(i,3)}.samprate;
    
    %%% get theta troughs
    dphase = [0; diff(phase{index(i,3)})];
    iTrough = find(dphase < -6); %itrough should be the index of the first negative value of the cycle
    iTrough = iTrough(iTrough<(length(amp{index(i,3)})-0.2*samprate)); %cut off troughs too close to the end of the recording
    iTrough = iTrough(iTrough>0.2*samprate); %has to be greater than 0.2s
    troughtimes{index(i,3)} = iTrough./samprate; %in sec
    
    
    ampTroughInds = [(iTrough - timeAroundTrough*samprate) (iTrough+timeAroundTrough*samprate)];
    tmptraces = arrayfun(@(x) amp{index(i,3)}(ampTroughInds(x,1):ampTroughInds(x,2)), 1:size(ampTroughInds,1), 'UniformOutput', false);
    allThetaTrace = [allThetaTrace; cell2mat(tmptraces')];
end
mnThetaTrace = mean(allThetaTrace,1); 

for t = 1:length(trialData)
    tmpSpeed = trialData(t).speed; 
    tmpTime = trialData(t).time;
    tmpPos = trialData(t).(posName);
    tmpFile = trialData(t).file; 
    
    %%% get trough times during the trial
    isTrialTrough = isExcluded(troughtimes{tmpFile}, [trialData(t).starttime, trialData(t).endtime]);
    isTrialTrough = logical(isTrialTrough);
    tmpTroughTimes = troughtimes{tmpFile}(isTrialTrough);
    tmpTroughVRInd = lookup2(tmpTroughTimes, tmpTime);

    %%% only troughs during moving periods
    isMoving = tmpSpeed > speedThreshold;
    movingInds = find(isMoving == 1);
%     isMovingTrough = ismember(tmpTroughVRInd, movingInds);
    runs = contiguous(isMoving,1);
    runs = runs{1,2};
    runtimes = tmpTime(runs(:,1));
    runtimes = [runtimes tmpTime(runs(:,2))];
    isMovingTrough = logical(isExcluded(tmpTroughTimes, runtimes));
    movingTroughTimes = tmpTroughTimes(isMovingTrough);
    movingTroughPos = tmpPos(tmpTroughVRInd(isMovingTrough));
    
    %%% get windows around the theta trough 
    timewindows = [movingTroughTimes-timeAroundTrough movingTroughTimes+timeAroundTrough];
    if isempty(timewindows)
        continue
    end
    
    %%% get spiking in the windows
    windowSize = timeAroundTrough*2/decodingBinsize;
    tmpTestingCounts = zeros(nCells, windowSize, size(timewindows,1)); 
    
    for w = 1:size(timewindows,1)
        winEdges = timewindows(w,1):decodingBinsize:timewindows(w,2); 
        for u = 1:nCells
            isCell = trialData(t).spikeIDs == inclCells(u);
            spkTimes = trialData(t).spikeTimes(isCell);
            tmpTestingCounts(u,:,w) = histcounts(spkTimes, winEdges);
        end
    end
    
    %%% bin full trial in small bins
%     tmpTimeBins = tmpTime(1):decodingBinsize:tmpTime(end); 
%     tmpTimeBinCenter = tmpTimeBins+decodingBinsize/2; 
%     tmpCenterVRInds = lookup2(tmpTimeBinCenter, tmpTime);
%     
%     isMovingCenter = ismember(tmpCenterVRInds, movingInds);
%     tmpCenterVRInds = tmpCenterVRInds(isMovingCenter);
%     tmpTimeBinCenterPos = tmpPos(tmpCenterVRInds);
%     
%     fullTrialTestingCounts = [];
%     for u = 1:nCells
%         isCell = trialData(t).spikeIDs == inclCells(u);
%         spkTimes = trialData(t).spikeTimes(isCell);
%         tmpCellH = histcounts(spkTimes, tmpTimeBins);
%         fullTrialTestingCounts(u,:) = tmpCellH(isMovingCenter);
%     end
%     
    %%% add to output structure
    testingData(t).counts = tmpTestingCounts; 
    testingData(t).troughTimes = movingTroughTimes;
    testingData(t).troughPos = movingTroughPos; 
%     testingData(t).fullTrialCounts = fullTrialTestingCounts;
%     testingData(t).fullTrialBinTimes = tmpTimeBinCenter(isMovingCenter);
%     testingData(t).fullTrialPos = tmpTimeBinCenterPos; 
    
    clear tmp* 
end

if isempty(testingData)
    A = 1; 
end



end

