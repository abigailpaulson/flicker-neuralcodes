function cf_plot_360decoding_unoccpiedzone(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_360decoding_unoccpiedzone
%
%ALP 4/16/24

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
filename = 'ThetaSeq_360decoding_trialdata_PCR_otherzonedecoding.mat';
load([seqdir, filename])

PlotData = OtherZonesTrial;
PlotDataDay = OtherZonesDay;

%%% change below when I add A
iPlot = 2;
iStat = 3;

%% get trial speed from unengaged vs. correct trials
figure(fh)

%%% lick latency
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.licklatency_s(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.licklatency_s(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'lick latency', 's', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([0 2])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick latency (s)')
title('Lick latency')
iPlot = iPlot+1;

%%% lick DI
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.LickDI(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.licklatency_s(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'lick discrimination', 'lick DI', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick discrimination index')
title('lick DI')
iPlot = iPlot+1;

%%% trial speed
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.Ctrl_speed(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.licklatency_s(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'trial speed', 'deg/s', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Speed (deg/s)')
title('trial speed')
iPlot = iPlot+1;


%%% per day information
% animal = DayData.animal;
% groups = DayData.group;
% avgDayDiff = DayData.avg_PCR_RRZ_difference;
% 
% axes(ax{iPlot})
% xplotvals = [1;2];
% for g = [1,2]
%     isGroup = strcmp(groups, gnames{g});
%     scatterdat{g} = avgDayDiff(isGroup);
%     mn = mean(avgDayDiff(isGroup), 'omitnan');
%     sem = std(avgDayDiff(isGroup), 'omitnan')./sqrt(sum(~isnan(avgDayDiff(isGroup))));
%     
%     plotprettypoints(fh, xplotvals(g,:), scatterdat(g), 'color', params.colors.(gnames{g}).(dnames{2}))
%     b = bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{2}));
%     b.FaceAlpha = 0.6;
%     er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
%     er.LineStyle = 'none';
%     
%     mn = []; sem = [];
% end
% cf_stats2txt2(scatterdat, statfid, panelL{iStat}, 'days', 'PCR difference engaged', 'PCR difference', gnames, tablefilename)
% xticks([1 2])
% xticklabels({'random', '40Hz'})
% ylim([-0.15 0.15])
% yticks([0-0.15 0 0.15])
% ylabel('prospective coding ratio difference')
% title('correct-unengaged within day')
% 





end

