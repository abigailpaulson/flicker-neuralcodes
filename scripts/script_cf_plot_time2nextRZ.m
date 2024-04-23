%%% Plot time to next reward zone - top and bottom speed thresholds
%
%ALP 3/6/2024

clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table.mat'];
load([filedir, filename])

%%% load theta seq data to get days to include
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
seqfilename = ['thetaseqdecoding_alldays_alltrials_230718.mat'];
load([seqdir, seqfilename])

%%% get rid of weird nans
isBadTrial = isnan(AllData.trialNum);
AllData = AllData(~isBadTrial,:); 

%%% get trials to include
inclData = AllData.significantSeq == 1 & allTMetrics.fullTrial==1 & allTMetrics.engaged ==1 & allTMetrics.rewarded == 1 & strcmp(allTMetrics.timepoint, 'post');
PlotData = allTMetrics(inclData,:);

%%% split data by median speed
ctrlSpeed = PlotData.Ctrl_speed;
qs = median(ctrlSpeed, 'omitnan');
qs = [min(ctrlSpeed) qs max(ctrlSpeed)];
iQ = discretize(ctrlSpeed, qs);

%%% plot
groupcolors = cbrewer('qual', 'Paired', 6); 
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
scolors = groupcolors;
gnames = {'gamma', 'random'};
dnames = {'pre', 'post'}; %here post is used for the color for fast trials

figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = 1:2
    isGroup = strcmp(PlotData.group, gnames{g});
    for s = 1:2
        isPlotTrial = isGroup & (iQ == s);
        vdat.([gnames{g}, '_', num2str(s)]) = PlotData.time_2_nextRZ(isPlotTrial);
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', scolors, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40);
xticklabels({'slow', 'fast', 'slow', 'fast'})
ylim([0 120])
ylabel('time to next reward zone (s)')
title('time to next reward zone')
makefigurepretty(gcf,1)
savefigALP(filedir, 'time2nextRZ_split_mediancontrolspeed')

%% save the data to run stats

StatsData.group = [PlotData.group];
StatsData.animal = [PlotData.animal];
StatsData.Ctrl_speed_subset = iQ;
StatsData.time_2_next_RZ = [PlotData.time_2_nextRZ];
StatsData = struct2table(StatsData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_Time2NextRZ_perTrial_postDataOnly_CtrlSpeed.txt';
writetable(StatsData, fullfile(statsdir, filename))

%% save the data for the figure
Time2NextRZ.group = [StatsData.group];
Time2NextRZ.animal = [StatsData.animal];
Time2NextRZ.ctrl_speed_subset = iQ;
Time2NextRZ.time_s = [StatsData.time_2_next_RZ]; 
Time2NextRZ = struct2table(Time2NextRZ); 

figdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\';
filename = 'Time2NextRZ_PostDataOnly_MedianSpeedSplit_PerTrial.mat';
save([figdatadir, filename], 'Time2NextRZ');
