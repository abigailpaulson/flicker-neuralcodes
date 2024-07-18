function cf_plot_PCR_engaged_unengaged(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_PCR_engaged_unengaged
%
%ALP 4/16/24

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
datadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ThetaSeq_Engaged_Unengaged_180Decoding_SupplementaryFig.mat';
load([datadir, filename])

PlotData = EngagedUnengagedData;

%% print helpful info about the analysis
greycolors = cbrewer('seq', 'Greys', 3);
greycolors = greycolors([2,3],:);

%% get trial speed from unengaged vs. correct trials
figure(fh)
iPlot = 1;
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for z = 1:2
    isPlotTrial = PlotData.isEngaged == (z-1); %set this up so that engaged and correct = 1 and unegaged = 0

    vdat.(['rewarded', num2str(z)]) = PlotData.Ctrl_speed(isPlotTrial);
    datforstats{z} = PlotData.Ctrl_speed(isPlotTrial);
end
cmat = greycolors; 
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([0 14])
xlim([0.5 1.5])
xticklabels({'Unengaged', 'Engaged'})
yticks([0 7 14])
ylabel('Speed (deg/s)')
title({'trial speed'})
cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'trials', 'speed', 'deg/s', {'unengaged', 'engaged'}, tablefilename)
iPlot = iPlot+1;

%%% PCR, both groups combined
axes(ax{iPlot})
%plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = []; datforstats = [];
for z = 1:2
    isPlotTrial = PlotData.isEngaged == (z-1);
    vdat.(['rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
    datforstats{z} = PlotData.PCR_loc_trial(isPlotTrial);
end
disp(['Both Groups - There are ', num2str(length(unique(PlotData.animal))), ' animals contributing to engaged/unengaged analysis'])
cmat = greycolors; 
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 1.5])
xticklabels({'Unengaged', 'Engaged'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'RRZ prospective coding'})
cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'trials', 'PCR engaged/un', 'prospective coding ratio', {'unengaged', 'engaged'}, tablefilename)
iPlot = iPlot + 1;

%%% both groups, engaged and unengaged

axes(ax{iPlot})
vdat = []; cmat = []; iStat = iPlot; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & PlotData.isEngaged == 0;
        else
            isPlotTrial = isGroup & PlotData.isEngaged == 1;
        end
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        datforstats{z} = PlotData.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'PCR engaged/un group', ...
        'prospective coding ratio', {[gnames{g}, ' unengaged'], [gnames{g}, ' engaged']}, tablefilename)
    
    disp([gnames{g} ' - There are ', num2str(length(unique(PlotData.animal(isGroup)))), ' animals contributing to engaged/unengaged analysis'])
    datforstats = [];
    iStat = iStat +1;
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'Unengaged', 'Engaged', 'Unengaged', 'Engaged'})
yticks([-0.75 0 0.75])
ylabel('Prospective coding ratio')
title({'RRZ prospective coding trials'})
iPlot = iPlot+1;

%%% per day information
animal = DayData.animal;
groups = DayData.group;
avgDayDiff = DayData.avg_PCR_RRZ_difference;

axes(ax{iPlot})
xplotvals = [1;2];
for g = [1,2]
    isGroup = strcmp(groups, gnames{g});
    scatterdat{g} = avgDayDiff(isGroup);
    mn = mean(avgDayDiff(isGroup), 'omitnan');
    sem = std(avgDayDiff(isGroup), 'omitnan')./sqrt(sum(~isnan(avgDayDiff(isGroup))));
    
    plotprettypoints(fh, xplotvals(g,:), scatterdat(g), 'color', params.colors.(gnames{g}).(dnames{2}))
    b = bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{2}));
    b.FaceAlpha = 0.6;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    disp([gnames{g} ' - There are ', num2str(length(unique(animal(isGroup)))), ' animals contributing to engaged/unengaged per day analysis'])

    
    mn = []; sem = [];
end
cf_stats2txt2(scatterdat, statfid, panelL{iStat}, 'days', 'PCR difference engaged', 'PCR difference', gnames, tablefilename)
xticks([1 2])
xticklabels({'random', '40Hz'})
ylim([-0.15 0.15])
yticks([0-0.15 0 0.15])
ylabel('prospective coding ratio difference')
title('correct-unengaged within day')






end

