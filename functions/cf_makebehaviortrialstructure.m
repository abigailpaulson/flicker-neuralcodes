function cf_makebehaviortrialstructure(dirs, params, allindex, metadata, btype)
%cf_makebehaviortrialstructure
%   make per trial behavior structure!
%   requires behaviorInfo and positionInfo from cf_getpositioninfo.m
%ALP 12/7/2022

rewrite.trialData = 1;

if strcmp(btype, '180')
    nz = 1;
    tracksize = 0:1:360;
    ndeg = 180;
else
    nz = 2;
    tracksize = 0:1:360;
    ndeg = 360;
end

positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
trackInfo = getTrackInfo_cflicker('chronicflicker_annulartrack');
dayindex = unique(allindex(:,1:2), 'rows');

for d = 1:size(dayindex,1)
    trialData = [];
    files = allindex(allindex(:,2) == dayindex(d,2), 3);
%     positionfilename = ['positionInfo_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    behaviorfilename = ['behaviorInfo_A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    trialfilename = ['trialInfo_',btype,'_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    
    if ~exist([positiondir, trialfilename], 'file') || rewrite.trialData
        disp(['getting trial data for ', num2str(dayindex(d,2))])
%         load([positiondir, positionfilename])
        load([positiondir, behaviorfilename])
        
        for i = 1:length(files)
            iP = i;
            iB = files(i);
            if isempty(behaviorInfo(iB).theta_d2r)
                continue
            end
%             time = round(positionInfo(iP).timeSmooth,4);
%             speed = positionInfo(iP).speedSmooth;
            time = round(behaviorInfo(iB).time,4);
            speed = behaviorInfo(iB).speed;
            theta_180 = behaviorInfo(iB).theta_d2r;
            theta = behaviorInfo(iB).theta; 
            licktimes = round(behaviorInfo(iB).lickTimes,4);
            [~, lickind] = ismember(licktimes, time);
            if ~isempty(licktimes) && sum(lickind) == 0
                dbstop
            end
            licks = ind2vect(lickind, length(time), 0);
            
            rewardtimes = round(behaviorInfo(iB).rewardTimes,4);
            [~, rewardind] = ismember(rewardtimes, time);
            rewards = ind2vect(rewardind, length(time), 0);
            
            if strcmp(btype, '180')
                startind = [1; find(diff(theta_180)<-150)+1];
            else
                startind = [1; find(diff(theta)<-150)+1];
            end
            endind = startind-1;
            endind = [endind(2:end); length(time)];
            
            for t = 1:length(startind)
                t_inds = startind(t):endind(t);
                if theta(t_inds(2)) < 180
                    zone = 1;
                else
                    zone = 2; 
                end
                AppendData(t).starttime = time(startind(t));
                AppendData(t).endtime = time(endind(t));
                AppendData(t).duration = time(endind(t)) - time(startind(t));
                AppendData(t).time = time(t_inds); 
                AppendData(t).theta = theta(t_inds);
                AppendData(t).theta_d2r = theta_180(t_inds);
                AppendData(t).speed = speed(t_inds);
                AppendData(t).licks = licks(t_inds);
                AppendData(t).rewards = rewards(t_inds);
                AppendData(t).nLicks = sum(licks(t_inds));
                AppendData(t).nRewards = sum(rewards(t_inds));
                AppendData(t).zone = zone;
                
                engaged = 0;
                if sum(licks(t_inds)) > 1
                    engaged = 1;
                end
                AppendData(t).engaged = engaged;
                
                %%% check for full trial
                roundPos = floor(AppendData(t).theta);
                passThrough = ismember(tracksize, roundPos);
                fulltrial = 0;
                if sum(passThrough) == ndeg
                    fulltrial = 1;
                end
                AppendData(t).fullTrial = fulltrial;
                
                %%% check for manual rewards or other reward issues
                if AppendData(t).nRewards > 4*nz || AppendData(t).nRewards < 3*nz
                    AppendData(t).rewarded = 0;
                else
                    AppendData(t).rewarded = 1;
                end
                
                %%% add file and animal information
                AppendData(t).an = dayindex(d,1);
                AppendData(t).day = dayindex(d,2);
                AppendData(t).file = files(i);
            end
            
            %%% append to trial data structure
            trialData = [trialData AppendData]; 
            AppendData = [];
        end
        
        info = []; 
        info = addhelpfulinfotostruct(info); 
        save([positiondir, trialfilename], 'trialData', 'info')
    end
end




end

