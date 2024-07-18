function cf_plot_behaviorperformance_overall_perday(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_behaviorperformance_overall_perday
%   plot behavior performance (proportion of correct trials) per day
%ALP 5/16/24

%% load datastructures
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\flicker-neuralcodes\results\behavior\']; 
savefilename = 'Behavior_PropCorrect_PerDay_SupplementaryFigure.mat';
load([savedatadir, savefilename])


%%% get group information from the info structure
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};
groupcolors = [params.colors];


%% plot
%%% plots for supplement figure 1
xplotvals = [1 2];

figure(fh)
axes(ax{1})
scatterdat_propcorrect = []; mn_propcorrect = []; sem_propcorrect = []; colororder = [];
hold on
for g = 1:2
    isGroup = strcmp(PerfData.group, gnames{g});
    GroupData = PerfData(isGroup,:);
    scatterdat_propcorrect{g} = GroupData.nCorrect./GroupData.nTotCorr;
    mn_propcorrect(g) = nanmean(scatterdat_propcorrect{g});
    sem_propcorrect(g) = nanstd(scatterdat_propcorrect{g})./sqrt(sum(~isnan(scatterdat_propcorrect{g})));
    colororder = [colororder; params.colors.(gnames{g}).(dnames{2})];
end
plotprettypoints(fh, xplotvals, scatterdat_propcorrect, 'color', colororder)
b = bar(xplotvals, mn_propcorrect, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_propcorrect, sem_propcorrect, 'Color','k');
er.LineStyle = 'none';
title('proportion correct')
ylabel('proportion of correct trials')
yticks([0 0.5 1])
xticks([])
cf_stats2txt2(scatterdat_propcorrect, statfid, panelL{1}, 'days', 'proportion', 'proportion correct trials', gnames, tablefilename)




end

