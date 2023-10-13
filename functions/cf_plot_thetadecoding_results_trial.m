function cf_plot_thetadecoding_results_trial(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_thetadecoding_results
%   plot theta sequence results here, prospective coding ratio over
%   position, the average value for the area of interest
%ALP 1/18/2023

dir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'GROUPEDDATA_decoding_thetaSeq_180_PC_alltrials.mat';

load([dir, filename])

%%% helpful plotting stuff 
xvect_bar = [1 2]; 
iPanel = 1;
RRZ = params.PCR_positions{1}; 

%%% figures

figure(fh)
axes(ax{1})
hold on
f = 8; 
% add a patch for the reward zone
fill([RRZ(1) RRZ(1) RRZ(2) RRZ(2)], [-0.15 0.15 0.15 -0.15], 'm', 'EdgeColor', 'none')
for d = 2
    for g = 1:length(dnames)
        tmpdat = GD.(gnames{g}).(dnames{d}).(fnames{f});
        plotdat_mn(g,:) = mean(tmpdat,1,'omitnan'); 
        plotdat_sem(g,:) = std(tmpdat, 0,1, 'omitnan')./sqrt(sum(~isnan(tmpdat(:,1)))); 
        shadedErrorBar(params.dec_edges(1:end-1)+3, plotdat_mn(g,:), plotdat_sem(g,:), {'Color', groupcolors.(gnames{g}).(dnames{d}), 'LineWidth', 2}, 1)
    end
end
xlabel('distance to reward zone (deg)')
xlim([-81 99])
xticks([-81 0 99])
ylabel('prospective coding ratio')
ylim([-0.15 0.15])
yticks([-0.15 -0.1 -0.05 0 0.05 0.1 0.15])
title({'prospective coding ratio', 'over position'})
iPanel = iPanel+1;

scatterdat = []; plotdat_mn = []; plotdat_sem = []; colororder = []; 
figure(fh)
axes(ax{2})
hold on
f = 10; 
for d = 2
    for g = 1:length(gnames)
        nm = [gnames{g}, dnames{d}];
        scatterdat{g} = GD.(gnames{g}).(dnames{d}).(fnames{f});
        vdat.(nm) = scatterdat{g};

        colororder(g,:) = groupcolors.(gnames{g}).(dnames{d});
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', colororder, 'ShowData', false', 'BoxWidth', 0.018, 'MedianSize', 25,'CenterSpace', 0.04);
ylabel('prospective coding ratio')
xticks([])
yticks([-0.75 0 0.75])
ylim([-0.75 0.75])
title({'Prospective coding ratio', 'in the reward-related zone'})
cf_stats2txt2(scatterdat, statfid, panelL{iPanel}, 'trials', 'prospective coding ratio', 'prospective coding ratio', gnames, tablefilename)
iPanel = iPanel+1;


end
