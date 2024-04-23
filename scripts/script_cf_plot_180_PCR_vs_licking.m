% script_cf_plot_180_PCR_vs_licking
%
%ALP 3/14/24

clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table_240416.mat'];
load([filedir, filename])

%%% load theta seq decoding data
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
seqfilename = ['thetaseqdecoding_alldays_alltrials_230718.mat'];
load([seqdir, seqfilename])

%%% set up things for plotting
groupcolors = cbrewer('qual', 'Paired', 6);
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
scolors = groupcolors;
gnames = {'gamma', 'random'};
dnames = {'pre', 'post'}; %here post is used for the color for fast trials
params.dec_edges = 0:6:360;
params.RZ = [54 72; 234 252];
params.RRZ = [36 64; 216 244];

%% plot
%%% get rid of weird nans
isBadTrial = isnan(AllData.trialNum);
AllData = AllData(~isBadTrial,:); 

%%% get trials to include
inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
PlotData = AllData(inclData,:);

BehavData = allTMetrics(inclData,:);

%%% get median lick count
splitname = 'nLicksAZ';
qs = median(BehavData.(splitname), 'omitnan');
qs = [min(BehavData.(splitname)) qs max(BehavData.(splitname))];
iQ = discretize(BehavData.(splitname), qs);

%%% plot scatterplot, lick number and PCR
figure
hold on
plot(BehavData.nLicksAZ, PlotData.PCR_loc_trial,'.')

% figure
% hold on
% iPlot = 1;
% for g = [2,1]
% subplot(1,2,iPlot)
% hold on
% isGroup = strcmp(PlotData.group, gnames{g});
% 
% plot(BehavData.nLicksAZ(isGroup), PlotData.PCR_loc_trial(isGroup),'.')
% 
% 
% iPlot = iPlot+1;
% end

%%% plot PCR by median licking

d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & iQ == z;
        vdat.([gnames{g}, 'lick', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'low licks', 'high licks', 'low licks', 'high licks'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'180 decoding', 'RRZ PCR by Lick #'})
makefigurepretty(gcf)
% filename = 'RRZ_PCR_pertrial_180decoding_nLickMedian';
% figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
% savefigALP(figdir, filename, 'filetype', 'pdf')



%% save the data to run stats

StatsData.group = [PlotData.group];
StatsData.animal = [PlotData.animal];
StatsData.Lick_number_subset = iQ;
StatsData.PCR_loc_trial = [PlotData.PCR_loc_trial];
StatsData = struct2table(StatsData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_RRZPCR_180_perTrial_postDataOnly_LickNumber.txt';
writetable(StatsData, fullfile(statsdir, filename))

yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'180 decoding', 'RRZ PCR by Lick #'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_nLickMedian';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'pdf')





