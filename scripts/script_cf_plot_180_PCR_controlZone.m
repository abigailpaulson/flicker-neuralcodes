% script_cf_plot_180_PCR_vs_licking
%
%ALP 3/14/24

clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table_v2.mat'];
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
params.posEdges = -81:3:99;
params.posBins = 3;
dType = '180';
params.nDeg = 180;
posName = 'theta_d2r';
params.RZ = [0 18];
params.AZ = [-36 0];
params.dec_edges = -81:9:99; %could also try 9

%set up window for calculating the quadrant value
% theta sequence parameters
params.decodingBinsize = 0.02; %in s
params.timeAroundTrough = 0.2; %in s
params.speedThreshold = 1; %deg/s
params.minTrials = 5;
params.adjustedEdges = -params.nDeg/2:params.posBins:params.nDeg/2;
params.timeEdges = -params.timeAroundTrough:params.decodingBinsize:params.timeAroundTrough;
params.quadrantWindow = 0.16; %in s

midI = round(length(params.timeEdges)/2);
binRange = params.quadrantWindow/(2*params.decodingBinsize);
LBinRange = [midI-binRange midI-1];
RBinRange = [midI midI+binRange-1]; %this one should include 0
posMidI = round(length(params.adjustedEdges)/2);
binRange = (length(params.adjustedEdges)-1)/4;
FBinRange = [posMidI posMidI+binRange-1];
BBinRange = [posMidI-binRange posMidI-1];

% quadrant indices, following the convention of Farooq and Dragoi 2019
params.qX = [RBinRange; LBinRange; LBinRange; RBinRange]; %time
params.qY = [FBinRange; FBinRange; BBinRange; BBinRange]; %position


%% get PCR in the control zone

%%% get rid of weird nans
isBadTrial = isnan(AllData.trialNum);
AllData = AllData(~isBadTrial,:); 

%%% get trials to calculate on
inclData = strcmp(AllData.timepoint, 'post');
PostData = AllData(inclData,:);

days = unique(PostData.day);

CtrlData = table;
for d = 1:length(days)
    disp(['getting ctrl data for ', num2str(days(d))])
    isDay = PostData.day == days(d);
    DayData = PostData(isDay,:);
    
    TmpCtrlData = cf_getcodingratio_thetaseq_table(DayData, [], params);
    
    CtrlData = [CtrlData; TmpCtrlData];
end

%% get trials to include
inclData = PostData.fullTrial==1 & PostData.engaged ==1 & PostData.rewarded == 1 & strcmp(PostData.timepoint, 'post') & PostData.significantSeq == 1;

PlotData = [PostData(inclData,2:5) PostData(inclData,35) CtrlData(inclData,4)];

%% plot control zone PCR vs. RRZ PCR
figdir = seqdir;

figure('Position', [440 575 404 223]);
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    isPlotTrial = isGroup;
    vdat.([gnames{g}, 'Ctrl']) = PlotData.PCR_ctrl_trial(isPlotTrial);
    vdat.([gnames{g}, 'RRZ']) = PlotData.PCR_loc_trial(isPlotTrial);
    cmat = [cmat; params.colors.(gnames{g}).(dnames{1}); params.colors.(gnames{g}).(dnames{2})];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'control zone', 'RRZ', 'control zone', 'RRZ'})
yticks([-0.75 0 0.75])
title('Control vs. RRZ PCR')
makefigurepretty(gcf)
filename = 'RRZ_vs_Ctrl_PCR_pertrial_180decoding';
savefigALP(figdir, filename, 'filetype', 'png')

%% get within trial difference to plot
PlotData.PCR_zonediff = PlotData.PCR_loc_trial - PlotData.PCR_ctrl_trial; 

figure('Position', [440 573 242 225]);
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    isPlotTrial = isGroup;
    vdat.([gnames{g}, 'diff']) = PlotData.PCR_zonediff(isPlotTrial);
    cmat = [cmat; params.colors.(gnames{g}).(dnames{2})];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 1.5])
xticklabels({'Random', '40 Hz'})
yticks([-0.75 0 0.75])
ylabel('difference in prospective coding ratio')
title('Control vs. RR Zone PCR')
makefigurepretty(gcf)
filename = 'RRZ_vs_Ctrl_PCR_pertrial_180decoding_withintrialdifference';
savefigALP(figdir, filename, 'filetype', 'png')

%% save data for running stats
clear StatsData
%%% set up per trial data for stats
StatsData.group = PlotData.group;
StatsData.animal = PlotData.animal; 
StatsData = struct2table(StatsData);
StatsData = [StatsData; StatsData];
tmpPCR = [PlotData.PCR_loc_trial];
tmpPCR2 = [PlotData.PCR_ctrl_trial];
tmpPCR3 = [tmpPCR; tmpPCR2];
StatsData.PCR_trial = tmpPCR3;
StatsData.PCRtype = [ones(length(tmpPCR),1); 2*ones(length(tmpPCR2),1)]; % 1 for RRZ, 2 for control zone

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_TrialData_180Decoding_ThetaSeq_RRZ_vs_Control.txt';
writetable(StatsData, fullfile(statsdir, filename))

%% save data for adding to the supplementary figure
ControlPCR = StatsData;

seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
filename = 'ThetaSeq_TrialData_180Decoding__RRZ_vs_Control_SupplementaryFigure.mat';
save([seqdir, filename], 'ControlPCR')



