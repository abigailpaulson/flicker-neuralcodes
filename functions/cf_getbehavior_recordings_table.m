function cf_getbehavior_recordings_table(dirs, params, allindex, metadata)
%cf_getbehavior_recordings
%   behavior analysis in the annular track during recording sessions
%   note - recording sessions only! 
%   make sure to run PreProcessing for TrialData structure before this
%   script
%ALP 12/7/22

%%% PARAMS
params.pos_edges = -81:3:99; 
params.RZ = [0 18];
params.AZ = [-18 0]; %ALP 8/16/23 changing this from [-9 0]. Actually changing this now 10/10/23
params.RRZ = [-18 18+18];
params.postRZ = [39 57];
MidBins = [find(params.pos_edges == -18) find(params.pos_edges == 18+18)-1]; %3 zones, AZ, RZ, 1 zone after
RZbins = [find(params.pos_edges == params.RZ(1)) find(params.pos_edges == params.RZ(2))-1];
AZbins = [find(params.pos_edges == params.AZ(1)) find(params.pos_edges == params.AZ(2))-1];
PostBins = [find(params.pos_edges == params.postRZ(1)) find(params.pos_edges == params.postRZ(2))-1];
VRsamprate = 50; %50 Hz sampling rate 

%%% helpful vectors for groups, etc
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 
dayID = 1:height(metadata);
gID.gamma.pre = dayID(isGamma & isPre);
gID.gamma.post = dayID(isGamma & isPost);
gID.random.pre = dayID(isRandom & isPre);
gID.random.post = dayID(isRandom & isPost);
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 

%%% intitialize, directories, etc. 
positiondir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\fig_behavior\';
dayindex = unique(allindex(:,1:2), 'rows');
allTData = []; allTMetrics = table;

%%%%% ----- calculate per trial behavior metrics ----- %%%%%
for d = 1:size(dayindex,1)
    trialfilename = ['trialInfo_A',  num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'];
    load([positiondir, trialfilename])
    nTrials = length(trialData);
    
    %%% inititalize
    lick_h = NaN(nTrials, length(params.pos_edges)-1);
    lick_nh = NaN(nTrials, length(params.pos_edges)-1);
    vel_h = NaN(nTrials, length(params.pos_edges)-1);
    vel_nh = NaN(nTrials, length(params.pos_edges)-1);
    duration = NaN(nTrials,1); 
    lickDI = NaN(nTrials,1);
    nLicks = NaN(nTrials,1); 
    nRewards = NaN(nTrials,1);
    licklatency_s = NaN(nTrials,1);
    licklatency_pos = NaN(nTrials,1);
    trial_speed = NaN(nTrials,1);
    AZ_speed = NaN(nTrials,1);
    RZ_speed = NaN(nTrials,1);
    Ctrl_speed = NaN(nTrials,1);
    postRZ_speed = NaN(nTrials,1);
    AZ_speed_slope = NaN(nTrials,1);
    AZ_lick_slope = NaN(nTrials,1);
    RRZ_time_below_thresh = NaN(nTrials,1); %ALP 5/23/23
    RRZ_time_below_thresh_prop = NaN(nTrials,1); %ALP 5/23/23
    enter_RZ = NaN(nTrials,1);
    leave_RZ = NaN(nTrials,1); 
    low_speed_pos = NaN(nTrials,1);
    high_speed_pos = NaN(nTrials,1);
    day = dayindex(d,2)*ones(nTrials,1);
    animal = dayindex(d,1).*ones(nTrials,1); 
    

    %%% analyze trials
    for t = 1:length(trialData)
        time = trialData(t).time; 
        licks = logical(trialData(t).licks);
        theta = trialData(t).theta_d2r;
        rewards = logical(trialData(t).rewards); 
        speed = trialData(t).speed; 
        lickpos = theta(licks); 
        rewardpos = theta(rewards); 
        duration(t) = time(end)-time(1);
        
        lick_h = histcounts(lickpos, params.pos_edges); 
        [~, ~, theta_bin] = histcounts(theta, params.pos_edges); 
        tmp_vel_h = arrayfun(@(x) mean(speed(theta_bin == x), 'omitnan'), 1:length(params.pos_edges)-1, 'UniformOutput', false); 
        tmp_vel_h = cell2mat(tmp_vel_h); 
        
        inRZ = theta >= params.RZ(1) & theta < params.RZ(2);
        inAZ = theta >= params.AZ(1) & theta < params.AZ(2); 
        inCtrl = ~inRZ & ~inAZ; 
        inRRZ = theta >= params.RRZ(1) & theta < params.RRZ(2);
        inPostRZ = theta >= params.postRZ(1) & theta < params.postRZ(2);
        
        RZtime = time(inRZ);
        if isempty(RZtime)
            RZtime = NaN;
        end
        licktime = time(licks); 
        RZlicktime = licktime(lickpos >= params.RZ(1) & lickpos < params.RZ(2));
        relativelicktime = RZlicktime - RZtime(1);
        lickDI(t) = (sum(licks(inRZ|inAZ))-sum(licks(inCtrl)))/(sum(licks));
        nLicks(t) = length(lickpos);
        nRewards(t) = length(rewards);
        
        if ~isempty(relativelicktime)
            licklatency_s(t) = min(relativelicktime(relativelicktime > 0));
            licklatency_pos(t) = min(lickpos(lickpos > 0 & lickpos < 18));
        end
        
        trial_speed(t) = mean(speed, 'omitnan'); 
        AZ_speed(t) = mean(speed(inAZ), 'omitnan');
        RZ_speed(t) = mean(speed(inRZ), 'omitnan');
        Ctrl_speed(t) = mean(speed(inCtrl), 'omitnan'); 
        postRZ_speed(t) = mean(speed(inPostRZ), 'omitnan');
        
        AZ_speed_slope(t) = (tmp_vel_h(AZbins(2)) - tmp_vel_h(AZbins(1)))/(AZbins(2)-AZbins(1));
        AZ_lick_slope(t) = (lick_h(AZbins(2)) - lick_h(AZbins(1)))/(AZbins(2)-AZbins(1));
        
        lick_nh(t,:) = lick_h'./sum(lick_h);
        vel_h(t,:) = tmp_vel_h';
        vel_nh(t,:) = vel_h(t,:)./max(vel_h(t,:));

        %%% smooth velocity and set a speed threshold (the mean), and get
        %%% time in the RZ below it 
        smooth_vel = smooth(speed, 50); 
        vel_thresh = mean(smooth_vel, 'omitnan');
        RRZ_time_below_thresh(t) = sum(speed(inRRZ) < vel_thresh)/VRsamprate; %in s 
        RRZ_time_below_thresh_prop(t) = RRZ_time_below_thresh(t)/duration(t);
        
        RRZ_speed = smooth_vel(inRRZ); 
        RRZ_pos = theta(inRRZ);
        belowThresh = RRZ_speed < vel_thresh; 
        pos_belowThresh = RRZ_pos(belowThresh); 
        if ~isempty(pos_belowThresh)
        low_speed_pos(t) = pos_belowThresh(1); 
        high_speed_pos(t) = pos_belowThresh(end);
        end
        
        
        %%% get time entering and leaving reward zone
        %then I can get the difference to get time to next reward zone
        enter_RZ(t) = RZtime(1);
        leave_RZ(t) = RZtime(end);  
        
%         figure
%         hold on
%         plot(speed)
%         plot(vel_thresh.*ones(1,length(speed)), '--')
%         
%         pause
%         close  
    end
    
    %get time to next reward zone
    tmp_enter_RZ = [enter_RZ(2:end); NaN];
    time_2_nextRZ = tmp_enter_RZ - leave_RZ;  
    
    trialNum = [1:length(trialData)]';
    engaged = [trialData.engaged]';
    fullTrial = [trialData.fullTrial]';
    rewarded = [trialData.rewarded]';
    group = repmat(metadata.Groups(d), nTrials, 1);
    group = convertCharsToStrings(group); 
    if metadata.FlickerDay(d) < 5
        timepoint = repmat({'pre'}, nTrials, 1);
    else
        timepoint = repmat({'post'}, nTrials, 1); 
    end
    timepoint = convertCharsToStrings(timepoint);
    zone = [trialData.zone]';
    
    metrics = table(lick_nh, vel_h, vel_nh, duration, lickDI, nLicks, ...
        nRewards, licklatency_s, licklatency_pos, trial_speed, AZ_speed, RZ_speed, ...
        Ctrl_speed, AZ_speed_slope, AZ_lick_slope, RRZ_time_below_thresh, RRZ_time_below_thresh_prop, enter_RZ, ...
        leave_RZ, time_2_nextRZ, low_speed_pos, high_speed_pos, postRZ_speed, zone, engaged, ...
        fullTrial, rewarded, trialNum, day, animal, group, timepoint);
    
    %%% append to large structures
    allTMetrics = [allTMetrics; metrics];
    allTData = [allTData trialData];
    
    savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
    savefilename = ['behavioranalysis_prepost_recordings_table_v2_', num2str(dayindex(d,2)), '.mat'];
    save([savedatadir, savefilename], 'metrics')
    
    clear metrics
end

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\'; 
%%% set up data structure for LMM in R
isRandom = strcmp(allTMetrics.group, 'random');
group = isRandom+1; %GROUP 1 is GAMMA and GROUP 2 is RANDOM!!!!!
isPost = strcmp(allTMetrics.timepoint, 'post');
timepoint = isPost+1; 
newInfo = table(group, timepoint);
behaviorData = [allTMetrics(:,3:end-2) newInfo];

inclTrial = allTMetrics.fullTrial & allTMetrics.engaged & allTMetrics.rewarded;
dayData = groupsummary(behaviorData(inclTrial,:), "day", "mean");
writetable(behaviorData, fullfile(statsdir, ['TableData_Behavior_Trial_v2.txt']));   
writetable(dayData, fullfile(statsdir, ['TableData_Behavior_Day_v2.txt']));   

units = {'lick count', 'raw speed (deg/s)', 'fraction of licks', 'norm speed (deg/s)', 'DI', 'latency (s)',...
    'latency (deg)', 'avg speed (deg/s)', 'avg speed (deg/s)', 'avg speed (deg/s)',...
    'avg speed (deg/s)', 'slope', 'slope'}; 

tmpdat1 = allTMetrics.trial_speed(inclTrial & strcmp(allTMetrics.group, 'gamma') & strcmp(allTMetrics.timepoint, 'pre'));
tmpdat2 = allTMetrics.trial_speed(inclTrial & strcmp(allTMetrics.group, 'random') & strcmp(allTMetrics.timepoint, 'pre'));
tmpdat3 = allTMetrics.trial_speed(inclTrial & strcmp(allTMetrics.group, 'gamma') & strcmp(allTMetrics.timepoint, 'post'));
tmpdat4 = allTMetrics.trial_speed(inclTrial & strcmp(allTMetrics.group, 'random') & strcmp(allTMetrics.timepoint, 'post'));
% 
% figure
% hold on
% raincloud_plot(tmpdat1, 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount', 0.15, 'dot_dodge_amount', 0.35, 'color', rgb('slateblue'), 'alpha', 0.5, 'box_col_match', 0, 'linewidth', 1)
% raincloud_plot(tmpdat2, 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount', 0.55, 'dot_dodge_amount', 0.75, 'color', rgb('mediumspringgreen'), 'alpha', 0.5, 'box_col_match', 0, 'linewidth', 1)
% makefigurepretty(gcf)
% 
% figure
% hold on
% subplot(2,1,1)
% raincloud_plot(tmpdat1, 'box_on', 1, 'color', rgb('slightly dark slateblue'), 'alpha', 0.5, 'box_col_match', 0)
% subplot(2,1,2)
% raincloud_plot(tmpdat2, 'box_on', 1, 'color', rgb('slightly dark mediumspringgreen'), 'alpha', 0.5, 'box_col_match', 0)
% makefigurepretty(gcf)


%%%%% ----- SAVE ----- %%%%%
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\']; 
info = [];
info = addhelpfulinfotostruct(info);
info.units = units;
savefilename = 'behavioranalysis_prepost_recordings_table_v2.mat';
save([savedatadir, savefilename], 'allTMetrics', 'allTData', 'info')

end

