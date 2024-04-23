% script_cf_plot_180_PCR_vs_licking
%
%ALP 3/14/24

clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\'];
filename = ['rippledecoding_halfratiozones_360_PC_AllData.mat'];
load([filedir, filename])
figdir = filedir;

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table_v2.mat'];
load([filedir, filename])

%%% set any necessary parameters
params.minRippleDur = 0.015; % in s
params.timeAroundMid = 0.125; % in s, on either side of mid
params.decodingBins = 0.025; %in s
params.posEdges = 0:2:360;
params.posBins = 2;
dType = '360';
params.nDeg = 360;
posName = 'theta';
params.RZ = [54 72; 234 252];
dec_edges = 0:6:360;

%% how often do ripples represent the other goal location?? 

inclRipples = AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos) & strcmp(AllData.timepoint, 'post');
PlotData = AllData(inclRipples,:);
pos = PlotData.rippleposition;

%%% replace the RZ ripple thing because its wrong for 360 degree
isRZ = (pos >= params.RZ(1,1) & pos < params.RZ(1,2)) | (pos >= params.RZ(2,1) & pos < params.RZ(2,2));
PlotData.RZRipple = isRZ; 

days = unique(PlotData.day);
for d = 1:length(days)
    isDay = PlotData.day == days(d);
    DayData = PlotData(isDay,:); 
    group = DayData.group{1};
    
    RippleData(d).group = group;
    RippleData(d).day = days(d);
    RippleData(d).nSigRZ = sum(DayData.RZRipple == 1 & DayData.rippleSigNonLocal == 1);
    RippleData(d).propSigRZ = sum(DayData.RZRipple == 1 & DayData.rippleSigNonLocal == 1)./sum(DayData.RZRipple == 1);
end

%% plotting
gnames = {'random', 'gamma'};
xvals = [1 2]; colororder = [];
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).propSigRZ];
    mn(g) = mean([RippleData(isGroup).propSigRZ], 'omitnan');
    sem(g) = std([RippleData(isGroup).propSigRZ], 'omitnan')./sqrt(sum(~isnan([RippleData(isGroup).propSigRZ])));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title({'RZ ripples', 'non-local zone decoding'})
ylabel('proportion of ripples')
xticks([1 2])
xticklabels(gnames)
makefigurepretty(gcf)
figname = 'prop_RZripples_sig_otherZone';
savefigALP(figdir, figname, 'filetype', 'png')

%% plotting
gnames = {'random', 'gamma'};
xvals = [1 2]; colororder = [];
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).nSigRZ];
    mn(g) = mean([RippleData(isGroup).nSigRZ], 'omitnan');
    sem(g) = std([RippleData(isGroup).nSigRZ], 'omitnan')./sqrt(sum(isGroup));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title({'RZ ripples', 'non-local zone decoding'})
ylabel('number of ripples')
xticks([1 2])
xticklabels(gnames)
makefigurepretty(gcf)
% figname = 'count_wrongtrials';
%savefigALP(figdir, figname, 'filetype', 'png')

%% get behavior data
%%% median speed buckets from theta sequence analysis
buckets = [0.4636 6.7630 13.4049];
speedType = discretize(allTMetrics.Ctrl_speed, buckets); 
AllTrialNonLocalCount = []; AllTrialCount = [];
for d = 1:length(days)
    isDay = PlotData.day == days(d);
    DayData = PlotData(isDay,:); 

    %%% match ripple times to trial subset
    isBehaviorDay = [allTData.day] == days(d);
    DayBehavior = allTData(:,isBehaviorDay);
    DaySpeedType = speedType(isBehaviorDay);
    
    trialEdges = [[DayBehavior.starttime]' [DayBehavior.endtime]'];
    rippleTimes = [DayData.rippleMidTimes];
    isRZ = DayData.RZRipple;
    rippleTimes = rippleTimes(isRZ);
    nonLocalRipple = DayData.rippleSigNonLocal(isRZ);
    rippleTrialNum = getbinindex(rippleTimes, trialEdges);
    

    %%% get info about # trials with ripples, etc. 
    trialEdgesCount = 1:1:length(DayBehavior)+1;
    [trialRippleCountsRaw, ~, bin] = histcounts(rippleTrialNum, trialEdgesCount);
    trialNonLocalCount = arrayfun(@(x) sum(nonLocalRipple(rippleTrialNum == x)), 1:length(DayBehavior));
    
    if length(trialNonLocalCount) ~= length(DayBehavior)
        pause
    end
    
    %%% get good trials
    isGoodTrial = [DayBehavior.engaged] & [DayBehavior.fullTrial] & [DayBehavior.rewarded];
    trialRippleCounts = trialRippleCountsRaw(isGoodTrial);
    goodTrialSpeed = DaySpeedType(isGoodTrial);
    trialNonLocalCountCorr = trialNonLocalCount(isGoodTrial);
    
    for s = 1:2
        ii = goodTrialSpeed == s;
        ripplePropSpeed(s,:) = sum(trialNonLocalCountCorr(ii))./sum(trialNonLocalCountCorr);
    end
    
    %%% get prop trials with > 1 ripple
    tmpprop = sum(trialRippleCounts >= 1)/length(trialRippleCounts);
    tmpn  = sum(trialRippleCounts >= 1); 
    tmpAvg = sum(trialRippleCounts)/length(trialRippleCounts);
    
    %%% add to ripple data
    RippleData(d).prop_trials_ripple = tmpprop;
    RippleData(d).n_trials_ripple = tmpn;
    RippleData(d).avg_ripples_trial = tmpAvg;
    RippleData(d).ripple_prop_nonlocal_speed = ripplePropSpeed;
    
    AllTrialNonLocalCount = [AllTrialNonLocalCount; trialNonLocalCount'];
    AllTrialCount = [AllTrialCount; trialRippleCountsRaw'];
    
    clear tmp* DayBehavior DayData ripplePropSpeed trialNonLocalCount
end

isIncludeDay = ismember(allTMetrics.day, days);
isPostBehavior = strcmp(allTMetrics.timepoint, 'post');
PostBehavior = allTMetrics(isIncludeDay & isPostBehavior,:);
PostBehavior.nonlocalRipples = AllTrialNonLocalCount;
PostBehavior.rippleCount = AllTrialCount;

isGoodTrial = [PostBehavior.engaged] & [PostBehavior.fullTrial] & [PostBehavior.rewarded];
PostBehavior = PostBehavior(isGoodTrial,:); 

%% plot proportion of trials with > 1 ripple
gnames = {'random', 'gamma'};
xvals = [1 2]; colororder = [];
fname = 'prop_trials_ripple';
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).(fname)];
    mn(g) = mean([RippleData(isGroup).(fname)], 'omitnan');
    sem(g) = std([RippleData(isGroup).(fname)], 'omitnan')./sqrt(sum(isGroup));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title({'prop trials > 1 SWR'})
ylabel('proportion of trials')
xticks([1 2])
xticklabels(gnames)
makefigurepretty(gcf)
figname = 'prop_trialcount_withripple_goodtrials';
savefigALP(figdir, figname, 'filetype', 'png')

%% number of trials with > 1 RZ ripple
gnames = {'random', 'gamma'};
xvals = [1 2]; colororder = [];
fname = 'n_trials_ripple';
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).(fname)];
    mn(g) = mean([RippleData(isGroup).(fname)], 'omitnan');
    sem(g) = std([RippleData(isGroup).(fname)], 'omitnan')./sqrt(sum(isGroup));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title({'n trials > 1 SWR'})
ylabel('number of trials')
xticks([1 2])
xticklabels(gnames)
makefigurepretty(gcf)
figname = 'n_trialcount_withripple_goodtrials';
savefigALP(figdir, figname, 'filetype', 'png')

%% avg number of ripples per trial
gnames = {'random', 'gamma'};
xvals = [1 2]; colororder = [];
fname = 'avg_ripples_trial';
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).(fname)];
    mn(g) = mean([RippleData(isGroup).(fname)], 'omitnan');
    sem(g) = std([RippleData(isGroup).(fname)], 'omitnan')./sqrt(sum(isGroup));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title({'avg ripples/trial'})
ylabel('average ripple count')
xticks([1 2])
xticklabels(gnames)
makefigurepretty(gcf)
figname = 'avg_ripplecount_RZripples_goodtrials';
savefigALP(figdir, figname, 'filetype', 'png')

%% proportion of trials with unoccupied reward zone decoding by fast and slow trials
gnames = {'random', 'gamma'};
xvals = [1 2; 4 5]; colororder = [];
fname = 'ripple_prop_nonlocal_speed';

fh = figure('Position', [440 543 560 255]);
hold on
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    tmpdat = [RippleData(isGroup).(fname)];
    scatterdat{1} = tmpdat(1,:);
    scatterdat{2} = tmpdat(2,:);
    mn = mean([RippleData(isGroup).(fname)], 2, 'omitnan');
    sem = std([RippleData(isGroup).(fname)], 0, 2, 'omitnan')./sqrt(sum(~isnan([RippleData(isGroup).(fname)]),2));
    colororder = [params.colors.(gnames{g}).pre; params.colors.(gnames{g}).post];
    
    plotprettypoints(fh, xvals(g,:), scatterdat, 'color', colororder)
    b = bar(xvals(g,:), mn, 'FaceColor', 'flat');
    er = errorbar(xvals(g,:), mn, sem, 'Color','k');
    b.CData = colororder;
    b.FaceAlpha = 0.6;
    er.LineStyle = 'none';
end
title({'prop ripples with unoccupied zone decoding'})
ylabel('proportion of ripples')
xticks([1 2 4 5])
xticklabels({'slow', 'fast', 'slow', 'fast'})
makefigurepretty(gcf)
figname = 'prop_ripples_otherRZ_RZripples_goodtrials_bytrialspeed';
savefigALP(figdir, figname, 'filetype', 'png')

%% plot trial behavior metrics for trials with or without significant nonlocal decoding
dnames = {'pre', 'post'};
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PostBehavior.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples == 0;
        else
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples >= 1;
        end
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PostBehavior.Ctrl_speed(isPlotTrial);
        disp(['unengaged trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
%ylim([-0.75 0.75])
%xlim([0.5 2.5])
xticklabels({'local', 'nonlocal', 'local', 'nonlocal'})
%yticks([-0.75 0 0.75])
ylabel('trial speed (deg/s)')
title({'trials with RZ ripples', 'decoding type'})
makefigurepretty(gcf)
filename = 'trialspeed_RZripples_local_nonlocal_360decoding';
savefigALP(figdir, filename, 'filetype', 'png')

%% plot trial behavior metrics for trials with or without significant nonlocal decoding
dnames = {'pre', 'post'};
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PostBehavior.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples == 0;
        else
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples >= 1;
        end
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PostBehavior.lickDI(isPlotTrial);
        disp(['unengaged trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
%ylim([-0.75 0.75])
%xlim([0.5 2.5])
xticklabels({'local', 'nonlocal', 'local', 'nonlocal'})
%yticks([-0.75 0 0.75])
ylabel('lick discrimination index')
title({'trials with RZ ripples', 'decoding type'})
makefigurepretty(gcf)
filename = 'lickDI_RZripples_local_nonlocal_360decoding';
savefigALP(figdir, filename, 'filetype', 'png')

%% plot trial behavior metrics for trials with or without significant nonlocal decoding
dnames = {'pre', 'post'};
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PostBehavior.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples == 0;
        else
            isPlotTrial = isGroup & PostBehavior.rippleCount >= 1 & PostBehavior.nonlocalRipples >= 1;
        end
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PostBehavior.licklatency_s(isPlotTrial);
        disp(['unengaged trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', true, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
%ylim([-0.75 0.75])
%xlim([0.5 2.5])
xticklabels({'local', 'nonlocal', 'local', 'nonlocal'})
%yticks([-0.75 0 0.75])
ylabel('lick latency (s)')
title({'trials with RZ ripples', 'decoding type'})
makefigurepretty(gcf)
filename = 'lickLatency_RZripples_local_nonlocal_360decoding';
savefigALP(figdir, filename, 'filetype', 'png')

%% save data for stats
clear StatsData
inclData = PostBehavior.rippleCount >= 1;
PostBehavior = PostBehavior(inclData,:);

StatsData.animal = PostBehavior.animal;
StatsData.group = PostBehavior.group;
StatsData.isNonLocal = double(PostBehavior.nonlocalRipples >= 1);
StatsData.licklatency_s = PostBehavior.licklatency_s;
StatsData.lickDI = PostBehavior.lickDI;
StatsData.Ctrl_speed = PostBehavior.Ctrl_speed; 

StatsData = struct2table(StatsData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_TrialData_360Decoding_Ripples_otherzone_pertrial.txt';
writetable(StatsData, fullfile(statsdir, filename))

%% save data for figures 

RippleDecodingUnoccupied = StatsData; 

savedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\';
filename = 'RippleDecoding_360Decoding_otherzone_pertrial_reviewerfiguredata.mat';
save([savedir, filename], 'RippleDecodingUnoccupied', 'RippleData'); 











