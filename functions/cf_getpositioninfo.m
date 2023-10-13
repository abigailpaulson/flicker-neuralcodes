function cf_getpositioninfo(allindex, virmendir, intandir, params)
%cf_getpositioninfo
% based on loadpositioninfo.m, getpositioninfo.m, and
% getVirmenPositionInfo.m
% for Abby's chronic flicker annular track experiment
%
%   inputs:
%       allindex - list of files, [animal date file#]
%       virmendir - base directory of virmen files
%       intandir - base directory of intan files
%ALP 12/1/2022

rewrite.rawpos = 0;
rewrite.positionInfo = 0;
rewrite.behaviorInfo = 1; 

positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
trackInfo = getTrackInfo_cflicker('chronicflicker_annulartrack'); 

dayindex = unique(allindex(:,1:2), 'rows');

for d = 1:size(dayindex,1)
    files = allindex(allindex(:,2) == dayindex(d,2), 3);
    anintandir = [intandir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
    positionfilename = ['positionInfo_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    behaviorfilename = ['behaviorInfo_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    
    for i = 1:length(files)
        anvirmendir = [virmendir, 'A', num2str(dayindex(d,1)), '\', 'A', ...
            num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '_', num2str(files(i)), '\'];
        index = [dayindex(d,1:2) files(i)];
        
        %%% rawpos
        if ~exist([anintandir, 'rawpos', num2str(files(i)), '.mat'], 'file') || rewrite.rawpos
            rawposALP(index, anvirmendir, anintandir);
        end
        
        if ~exist([anintandir, 'rawpos', num2str(files(i)), '.mat'], 'file')
            rawInfo(i).thetaRaw = [];
            rawInfo(i).timeRaw = [];
            rawInfo(i).thetaBreaks = []; 
            continue
        end
        
        load([anintandir, 'rawpos', num2str(files(i)), '.mat'])

        %%% compile raw positionInfo
        rawInfo(i).thetaRaw = rawpos{dayindex(d,1)}{dayindex(d,2)}{files(i)}.theta;
        rawInfo(i).timeRaw = rawpos{dayindex(d,1)}{dayindex(d,2)}{files(i)}.indices/20000;
        switchTheta = find(diff(rawInfo(i).thetaRaw) < -300);
        thetaBreaks = [0; switchTheta; length(rawInfo(i).thetaRaw)];
        rawInfo(i).thetaBreaks = thetaBreaks;
    end
    
    clear rawPos
    %%% load raw pos for all files
    rawpos = loaddatastruct2(anintandir, dayindex(d,:), 'rawpos', files);
    
    if ~exist([positiondir, positionfilename], 'file') || rewrite.positionInfo
        %%% smooth position
        positionInfo = smoothPosition(rawInfo, files);
        
        %%% find when moving
        positionInfo = findWhenMoving(files, positionInfo, params.speedThreshold);
        
        %%% create other position vectors (distance to reward)
        positionInfo = getdistance2reward(positionInfo, trackInfo, files);
        
        %%% populate index field
        for i = 1:length(files)
            positionInfo(i).index = [dayindex(d,:) files(i)];
        end
    else
        load([positiondir, positionfilename])
    end

    %%% append relevant behavior
    for i = 1:length(files)
        fID = files(i); 
        %%% check files match positionInfo files
        if fID ~= positionInfo(i).index(1,3)
            error('position info file index doesnt match')
        end
        
        if ~isempty(positionInfo(i).timeSmooth)
            posTime = positionInfo(i).timeSmooth;
            speed = positionInfo(i).speedSmooth;
            fulltheta = positionInfo(i).thetaSmooth; 
            tmp = positionInfo(i).distanceToReward;
            tmp = tmp+18;
            tmp(tmp<-81) = tmp(tmp<-81)+180; %this should make it so that there is -81 to 99 with the reward zone at 0 to 18
            tmpBehavior(i).d2r_centered = tmp;
            
            vrTimesRaw = rawpos{dayindex(d,1)}{dayindex(d,2)}{files(i)}.indices./20000;
            vrTimesSmooth = round(positionInfo(i).timeSmooth,6);
            
            isReward = find(diff(rawpos{dayindex(d,1)}{dayindex(d,2)}{files(i)}.reward))+1;
            tmpRewardTime = vrTimesRaw(isReward);  %non resampled
            tmpInds = lookup2(tmpRewardTime, vrTimesSmooth);
            rewardTimes = positionInfo(i).timeSmooth(tmpInds);
            rewardPos = positionInfo(i).thetaSmooth(tmpInds);
            rewardPos_d2r = tmpBehavior(i).d2r_centered(tmpInds);
            
            isLick = find(diff(rawpos{dayindex(d,1)}{dayindex(d,2)}{files(i)}.licks))+1;
            tmpLickTime = vrTimesRaw(isLick);
            tmpIndsLick = lookup2(tmpLickTime, vrTimesSmooth);
            lickTimes = [positionInfo(i).timeSmooth(tmpIndsLick)];
            lickTimes = round(lickTimes,6);
            lickPos = positionInfo(i).thetaSmooth(tmpIndsLick);
            lickPos_d2r = tmpBehavior(i).d2r_centered(tmpIndsLick);
        else
            speed = [];
            posTime = [];
            tmpBehavior(i).d2r_centered = [];
            fulltheta = []; 
            lickTimes = [];
            lickPos = [];
            lickPos_d2r = [];
            rewardTimes = [];
            rewardPos = [];
            rewardPos_d2r = [];
        end
        
        behaviorInfo(fID).time = posTime; 
        behaviorInfo(fID).speed = speed;
        behaviorInfo(fID).theta = fulltheta; 
        behaviorInfo(fID).theta_d2r = tmpBehavior(i).d2r_centered; 
        behaviorInfo(fID).lickTimes = lickTimes;
        behaviorInfo(fID).lickPos = lickPos;
        behaviorInfo(fID).lickPos_d2r = lickPos_d2r; 
        behaviorInfo(fID).rewardTimes = rewardTimes;
        behaviorInfo(fID).rewardPos = rewardPos; 
        behaviorInfo(fID).rewardPos_d2r = rewardPos_d2r; 
        behaviorInfo(fID).index = [dayindex(d,:) fID]; 
    end
    
    %%% save positionInfo
    saveinfo = []; 
    saveinfo = addhelpfulinfotostruct(saveinfo); 
    if ~exist(positiondir, 'dir'); mkdir(positiondir); end
    
    if rewrite.positionInfo || ~exist([positiondir, positionfilename], 'file')
        save([positiondir, positionfilename], 'positionInfo', 'saveinfo');
    end
    
    if rewrite.behaviorInfo || ~exist([positiondir, behaviorfilename],'file')
        save([positiondir, behaviorfilename], 'behaviorInfo', 'saveinfo'); 
    end
    
    clear behaviorInfo positionInfo rawpos
end

end

