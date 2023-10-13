function cf_plot_prospectivecoding_byspeed_supplement(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
% plot supplement info for supp figure 5
%ALP 06/30/23

%% set up colors, parameters, etc
datadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'figure5_data.mat';
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

%%
iPlot = 1; 

% violin plot of control speeds across groups
axes(ax{iPlot})
hold on
vdat = []; cmat = [];
for s = 1:2
    for g = 1:2
        isGroup = strcmp(postRPCR.group, gnames{g});
        isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
        vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.Ctrl_speed(isPlotTrial);
        cmat = [cmat; c.(gnames{g}).(snames{s})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
xticklabels({'low', 'high', 'low', 'high'})
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
ylabel('velocity (deg/s)')
title('control speed')


%I made this into its own function...
%cf_plot_RRZprospectivecoding_byspeed.m ALP 7/13/2023

% %   prospective coding in the reward zone comparison between groups split
% %   by speed
% 
% iPlot = iPlot+1;
% axes(ax{iPlot})
% hold on
% vdat = []; cmat = [];
% for s = 1:2
%     for g = 1:2
%         isGroup = strcmp(postRPCR.group, gnames{g});
%         isPlotTrial = isGroup & (postRPCR.Ctrl_speed_subset == s);
%         vdat.([gnames{g}, '_', num2str(s)]) = postRPCR.PCR_loc_trial(isPlotTrial);
%         cmat = [cmat; c.(gnames{g}).(snames{s})];
%     end
% end
% violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% % violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
% xticklabels({'low', 'high', 'low', 'high'})
% % ylim([-0.5 0.75])
% % yticks([-0.5 0 0.5])
% ylabel('prospective coding ratio')
% title('reward related zone prospective coding')


end

