% script_cf_plot_180_PCR_unengaged
%
% plot unengaged trials, don't worry about incorrect
%
%   all groups together, and then split by group
%
% also prep data structure for running stats
%
%ALP 4/2/24

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
params.dec_edges = 0:6:360;
params.RZ = [54 72; 234 252];
params.RRZ = [36 64; 216 244];

%% prep the data for plotting - only days that have unengaged trials
%%% remove weird nan values in All Data
isnananimal = isnan(AllData.animal);
AllData = AllData(~isnananimal,:); 

%%% get trials to include
inclData = AllData.fullTrial==1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
PlotData = AllData(inclData,:);
FullTrialMetrics = allTMetrics(inclData,:);

%%% identify trials with licks > 0 in the reward zone
RZbins = 28:33; 
for t = 1:height(FullTrialMetrics)
     licksinRZ(t) = sum(FullTrialMetrics(t,:).lick_nh(RZbins)) > 0;
end

FullTrialMetrics.licksinRZ = licksinRZ'; 
PlotData.licksinRZ = licksinRZ';

%%% get days with both unengaged and correct trials - only want to include
%%% those days for plotting
days = unique(PlotData.day);
for d = 1:length(days)
    isDay = PlotData.day == days(d);
    nCorrect = sum(PlotData(isDay,:).engaged & PlotData(isDay,:).rewarded);
    nUnengaged = sum(~PlotData(isDay,:).engaged);
    groups{d} = unique(PlotData(isDay,:).group);
    hasBoth(d) = nCorrect > 0 & nUnengaged > 0;
end
hasBothDays = days(hasBoth); 
hasBothGroups = groups(hasBoth);
inclDays = ismember(PlotData.day, hasBothDays);

PlotData = PlotData(inclDays,:);


%% first plot unengaged trials vs correct trials - both groups combined
greycolors = cbrewer('seq', 'Greys', 3);
greycolors = greycolors([2,3],:);

figure('Position', [409 393 207 226])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for z = 1:2
    isPlotTrial = PlotData.engaged == (z-1) & PlotData.rewarded == (z-1);

    vdat.(['rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
    disp(['unengaged vs. correct trials - both groups ', 'z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
end
cmat = greycolors; 
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 1.5])
xticklabels({'unengaged', 'correct'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'Reward related prospective coding', 'unengaged trials'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_correct_unengaged_trials_bothdays_combinedgroups';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')


%% get trial speed from unengaged vs. correct trials
figure('Position', [409 393 207 226])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for z = 1:2
    isPlotTrial = PlotData.engaged == (z-1) & PlotData.rewarded == (z-1);

    vdat.(['rewarded', num2str(z)]) = PlotData.Ctrl_speed(isPlotTrial);
    disp(['unengaged vs. correct trials - both groups ', 'z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
end
cmat = greycolors; 
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
%ylim([-0.75 0.75])
xlim([0.5 1.5])
xticklabels({'unengaged', 'correct'})
%yticks([-0.75 0 0.75])
ylabel('trial speed (deg/s)')
title({'trial speed', 'unengaged trials'})
makefigurepretty(gcf)
filename = 'CtrlSpeed_180decoding_correct_unengaged_trials_bothdays_combinedgroups';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')

%% plot unengaged trials - both groups
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & PlotData.engaged == 0 & PlotData.rewarded == (z-1);
        else
            isPlotTrial = isGroup & PlotData.engaged == 1 & PlotData.rewarded == (z-1);
        end
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        disp(['unengaged trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'unengaged', 'correct', 'unengaged', 'correct'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'Reward related prospective coding', 'unengaged trials'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_correct_unengaged_trials_bothdays';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')

%% get within day difference
BothDayData = PlotData; 

days = unique(BothDayData.day); groups = [];
for d = 1:length(days)
    isDay = BothDayData.day == days(d);
    isCorrect = isDay & BothDayData.engaged == 1 & BothDayData.rewarded == 1;
    isUnengaged = isDay & BothDayData.engaged == 0;
    animal(d) = BothDayData(isDay,:).animal(1); 
    nCorrect = sum(isCorrect);
    avgCorrect(d) = mean(BothDayData(isCorrect,:).PCR_loc_trial, 'omitnan');
    nUnengaged = sum(isUnengaged);
    avgUnengaged(d) = mean(BothDayData(isUnengaged,:).PCR_loc_trial, 'omitnan');
    avgDayDiff(d) = avgCorrect(d) - avgUnengaged(d);
    groups{d} = BothDayData(isDay,:).group(1); 
end
groups = cellstr(groups);

fh = figure('Position', [353 329 316 255]);
hold on
xplotvals = [2;1];
for g = [2,1]
    isGroup = strcmp(groups, gnames{g});
    scatterdat{1} = avgDayDiff(isGroup);
    mn = mean(avgDayDiff(isGroup), 'omitnan');
    sem = std(avgDayDiff(isGroup), 'omitnan')./sqrt(sum(~isnan(avgDayDiff(isGroup))));
    
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', params.colors.(gnames{g}).(dnames{2}))
    b = bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{2}));
    b.FaceAlpha = 0.6;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    mn = []; sem = []; scatterdat = [];
end
xticks([1 2])
xticklabels({'random', '40Hz'})
ylim([-0.15 0.15])
yticks([0-0.15 0 0.15])
ylabel('prospective coding ratio difference')
title('correct-unengaged within day')
makefigurepretty(gcf)
filename = 'perday_difference_correct_unengaged_RRZ_PCR';
savefigALP(figdir, filename, 'filetype', 'png')

fh = figure('Position', [353 328 474 256]);
hold on
xplotvals = [4 5; 1 2];
for g = [2,1]
    isGroup = strcmp(groups, gnames{g});
    scatterdat{1} = avgUnengaged(isGroup);
    scatterdat{2} = avgCorrect(isGroup);
    
    mn = [mean(avgUnengaged(isGroup), 'omitnan') mean(avgCorrect(isGroup), 'omitnan')];
    sem = [std(avgUnengaged(isGroup), 'omitnan')./sqrt(sum(~isnan(avgUnengaged(isGroup)))) std(avgCorrect(isGroup), 'omitnan')./sqrt(sum(~isnan(avgCorrect(isGroup))))];
    cmat = [params.colors.(gnames{g}).(dnames{1}); params.colors.(gnames{g}).(dnames{2})];
    
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', cmat)
    b = bar(xplotvals(g,:), mn);
    b.FaceAlpha = 0.6;
    b.FaceColor = 'flat';
    b.CData = cmat;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    mn = []; sem = []; scatterdat = [];
end
xticks([1 2 4 5])
xticklabels({'unengaged', 'correct', 'unengaged', 'correct'})
ylim([-0.15 0.15])
yticks([-0.15 0 0.15])
ylabel('prospective coding ratio')
title('day avg reward related zone PCR')
makefigurepretty(gcf)
filename = 'perday_correct_unengaged_RRZ_PCR';
savefigALP(figdir, filename, 'filetype', 'png')

%% save data for stats

%%% set up per trial data for stats
isUnengaged = PlotData.engaged == 0 & PlotData.rewarded == 0;
isCorrect = PlotData.engaged == 1 & PlotData.rewarded == 1; 
inclTrial = isUnengaged | isCorrect; 
PlotData = PlotData(inclTrial,:);

StatsData.group = PlotData.group;
StatsData.animal = PlotData.animal;
StatsData.isUnengaged = PlotData.engaged == 0 & PlotData.rewarded == 0;
StatsData.isEngaged = PlotData.engaged == 1 & PlotData.rewarded == 1;
StatsData.PCR_loc_trial = PlotData.PCR_loc_trial;
StatsData.Ctrl_speed = PlotData.Ctrl_speed;
StatsData = struct2table(StatsData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_TrialData_180Decoding_ThetaSeq_Unengaged_Correct_Trials.txt';
writetable(StatsData, fullfile(statsdir, filename))

%%% set up per day data for stats
DayData.animal = animal';
DayData.group = groups';
DayData.avg_PCR_RRZ_Correct = avgCorrect';
DayData.avg_PCR_RRZ_Unengaged = avgUnengaged';
DayData.avg_PCR_RRZ_difference = avgDayDiff'; 

DayData = struct2table(DayData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_DayData_180Decoding_ThetaSeq_Unengaged_Correct.txt';
writetable(DayData, fullfile(statsdir, filename))

%% save data for figures

EngagedUnengagedData = StatsData;
filename = 'ThetaSeq_Engaged_Unengaged_180Decoding_SupplementaryFig.mat';
save([seqdir, filename], 'EngagedUnengagedData', 'DayData')











