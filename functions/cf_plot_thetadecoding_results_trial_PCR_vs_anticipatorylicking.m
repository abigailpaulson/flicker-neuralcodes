function cf_plot_thetadecoding_results_trial_PCR_vs_anticipatorylicking(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_thetadecoding_results_trial_PCR_vs_anticipatorylicking
%  PCR RRZ split by median anticipatory licking
%ALP 4/29/24

dir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ThetaSeq_TrialData_180Decoding__RRZ_by_AnticipatoryLickMedian_SupplementaryFigure.mat';

load([dir, filename])
PlotData = PCRAnticipatoryLicking;

%%% helpful plotting stuff 
iPlot = 1;
iStat = 1;
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};
subsetnames = {'low licks', 'high licks'};

%%% figures
figure(fh)
axes(ax{iPlot})

%plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for s = [1:2]
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & PlotData.Lick_number_subset == s;
        datforstats{s} = PlotData.PCR_loc_trial(isPlotTrial);
        vdat.([gnames{g}, num2str(s)]) = PlotData.PCR_loc_trial(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'prospective coding ratio', ...
        [gnames{g} ' low licks vs high licks'], subsetnames, tablefilename)
    cmat = [cmat; params.colors.(gnames{g}).(dnames{1}); params.colors.(gnames{g}).(dnames{2})];
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'ViolinAlpha', 0.4,'MedianSize', 25);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'Low licks', 'High licks', 'Low licks', 'High licks'})
yticks([-0.75 0 0.75])
title('Control vs. RRZ PCR')


end
