function cf_plot_thetadecoding_results_trial_vs_control(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_thetadecoding_results_trial_vs_control
%  PCR RRZ vs control zone
%ALP 4/16/2024

dir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ThetaSeq_TrialData_180Decoding__RRZ_vs_Control_SupplementaryFigure.mat';

load([dir, filename])
PlotData = ControlPCR;

%%% helpful plotting stuff 
iPlot = 1;
iStat = 1;
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};
PCRnames = {'RRZ', 'Control'};

%%% figures
figure(fh)
axes(ax{iPlot})

%plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for s = [2,1]
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & PlotData.PCRtype == s;
        datforstats{s} = PlotData.PCR_trial(isPlotTrial);
        vdat.([gnames{g}, num2str(s)]) = PlotData.PCR_trial(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'prospective coding ratio', [gnames{g} ' PCR Ctrl vs RRZ'], PCRnames, tablefilename)
    cmat = [cmat; params.colors.(gnames{g}).(dnames{1}); params.colors.(gnames{g}).(dnames{2})];
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'Control', 'Rewarded', 'Control', 'Rewarded'})
yticks([-0.75 0 0.75])
title('Control vs. RRZ PCR')


end
