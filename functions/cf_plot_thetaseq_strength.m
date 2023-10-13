function cf_plot_thetaseq_strength(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_thetaseq_strength
%
%ALP 4/13/23

%%% data dir and load
dir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\'];
load([dir, 'thetaseqdecoding_includeddays_group.mat'])

inclDay = strcmp(dayData.timepoint, 'post');
postData = dayData(inclDay,:); 

%%% helpful plotting stuff 
xvect_bar = [1 2]; 
iPanel = 1;
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 
groupcolors = params.colors;


%%% plot the data
figure(fh)
axes(ax{iPanel})
hold on
for d = 2
    for g = 1:length(gnames)
        isGroup = strcmp(postData.group, gnames{g});
        tmpdat = postData.mean_dayQRatio(isGroup);
        scatterdat{g} = tmpdat;
        plotdat_mn(g) = mean(scatterdat{g}, 'omitnan');
        plotdat_sem(g) = std(scatterdat{g}, 'omitnan')./sqrt(sum(~isnan(scatterdat{g})));
        colororder(g,:) = groupcolors.(gnames{g}).(dnames{d});
        
         %%% how many animals contribute to theta sequence strength?
            nAn = length(unique(postData.animal(isGroup)));
        disp(['There are ', num2str(nAn), ' animals contributing to theta sequence strength day for ', gnames{g}])
    
    end
end
plotprettypoints(fh, xvect_bar, scatterdat, 'color', colororder)
b = bar(xvect_bar, plotdat_mn, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
errorbar2(xvect_bar, plotdat_mn, plotdat_sem, 0.2, 'k-', 'LineWidth', 0.75);
ylabel('quadrant ratio')
xticks([])
yticks([0 0.05 0.1])
ylim([0 0.1])
title({'Overall theta sequence strength'})
cf_stats2txt2(scatterdat, statfid, panelL{iPanel}, 'days', 'quadrant ratio', 'quadrant ratio', gnames, tablefilename)
iPanel = iPanel+1; 



end

