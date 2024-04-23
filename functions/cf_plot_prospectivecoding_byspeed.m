function cf_plot_prospectivecoding_byspeed(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%plot prospective coding throughout the trial by speed subsets
%
%ALP 6/29/23

%% set up colors, parameters, etc
datadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'figure5_data_230723.mat';
load([datadir, filename])
postR_group = grpstats(postRPCR, {'group', 'Ctrl_speed_subset'}, {'mean', 'sem'});
behavioredges = -81:3:99;
params.dec_edges = -81:9:99; %could also try 9


gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};
snames = {'speed1', 'speed2'};
groupcolors = params.colors; 


%%% custom colors for these plots
GS1 = hex2rgb('80aec5'); %gamma speed subset 1
GS2 = hex2rgb('196598'); %gamma speed subset 2 
RS1 = hex2rgb('a2cb7d'); %random speed subset 1
RS2 = hex2rgb('2a8924'); %random speed subset 2

c.gamma.speed1 = GS1;
c.gamma.speed2 = GS2;
c.random.speed1 = RS1;
c.random.speed2 = RS2; 

%% figures
%%% average of high and low speeds over trials
iPlot = 1;

figure(fh)
hold on
tmpplot = 1;
for g = 1:2
    axes(ax{iPlot})
    hold on
    patch([-81 -9 -9 -81], [0 0 12 12], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    patch([18 99 99 18], [0 0 12 12], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    for s = 1:2
        mn = []; sem = [];
        nm = [gnames{g}, '_', num2str(s)];
        mn = postR_group(nm,:).mean_velocity;
        sem = postR_group(nm,:).sem_velocity;
        shadedErrorBar(behavioredges(1:end-1), mn, sem, {'Color', c.(gnames{g}).(snames{s}), 'LineWidth', 1},1);
        datforstats{s} = mn(:,1);
        tmpplot = tmpplot+1;
    end
    xlabel('Distance to reward zone (deg)')
    ylabel('Speed (deg/s)')
    ylim([2 12])
    yticks([2 7 12])
    xticks([-81 0 99])
    xlim([-81 99])
    %cf_stats2txt2(datforstats, statfid, panelL{iPanel}, 'trials', 'speed', ['speed', ' ', gnames{g}, snames, tablefilename])

    iPlot = iPlot+1;
end

% speed per trial split by group
iStat = 1;
axes(ax{iPlot})
hold on
vdat = []; cmat = [];
for s = 1:2
    datforstats{g} = [];
    for g = 1:2
        isGroup = strcmp(postRPCR.group, gnames{g});
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.Ctrl_speed(isPlotTrial);
        cmat = [cmat; c.(gnames{g}).(snames{s})];
        datforstats{g} =  postRPCR.Ctrl_speed(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'deg/s', 'speedovertrack', gnames, tablefilename)
    iStat = iStat+1;
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'slow', 'fast', 'slow', 'fast'})
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
ylabel('velocity (deg/s)')
title('control speed')
iPlot  = iPlot + 1;

% prospective coding over position split by velocity - per trial
for g = 1:2
    axes(ax{iPlot})
    hold on
    patch([18 72 72 18], [0 0 0.25 0.25], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
    for s = 1:2
        mn = []; sem = [];
        nm = [gnames{g}, '_', num2str(s)];
        mn = postR_group(nm,:).mean_PCR_pos_trial;
        sem = postR_group(nm,:).sem_PCR_pos_trial;
        shadedErrorBar(params.dec_edges(1:end-1), mn, sem, {'Color', c.(gnames{g}).(snames{s}), 'LineWidth', 1},1);
    end
    xlabel('Distance to reward zone (deg)')
    ylabel('Prospective coding ratio')
    ylim([-0.25 0.25])
    yticks([-0.25 0 0.25])
    xticks([-81 0 99])
    xlim([-81 99])
    %cf_stats2txt2(datforstats, statfid, panelL{iPanel}, 'trials', 'speed', ['speed', ' ', gnames{g}, snames, tablefilename])

    iPlot = iPlot + 1;
end


% make this a violin plot for the main figure
axes(ax{iPlot})
vdat = []; cmat = [];
for g = 1:2
    isGroup = strcmp(postRPCR.group, gnames{g});
    datforstats = [];
    for s = 1:2
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.PostR_PCR_loc_trial(isPlotTrial);
        cmat = [cmat; c.(gnames{g}).(snames{s})];
        datforstats{s} = postRPCR.PostR_PCR_loc_trial(isPlotTrial);
    end
 cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'prospective coding ratio', 'postRewPCR_speed', snames, tablefilename)
iStat = iStat+1;
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'slow', 'fast', 'slow', 'fast'})
ylim([-0.5 0.75])
yticks([-0.5 0 0.5])
ylabel('Prospective coding ratio')
title('post reward PCR')
iPlot = iPlot + 1;

% plot time to next reward zone
% load data, from 'script_cf_plot_time2nextRZ.m'
figdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\';
filename = 'Time2NextRZ_PostDataOnly_MedianSpeedSplit_PerTrial.mat';
load([figdatadir, filename])

axes(ax{iPlot})
hold on
vdat = []; cmat = []; datforstats = [];
for g = 1:2
    isGroup = strcmp(Time2NextRZ.group, gnames{g});
    for s = 1:2
        isPlotTrial = isGroup & (Time2NextRZ.ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = Time2NextRZ.time_s(isPlotTrial);
        datforstats{s} = Time2NextRZ.time_s(isPlotTrial);
        cmat = [cmat; c.(gnames{g}).(snames{s})];
    end
     cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'time2NextRZ', 'Ctrl_speed', snames, tablefilename)
    iStat = iStat + 1;
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
xticklabels({'slow', 'fast', 'slow', 'fast'})
ylim([0 120])
ylabel('Time (s)')
title('Time to next reward zone')












end

