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
params.dec_edges = 0:6:360;
params.RZ = [54 72; 234 252];
params.RRZ = [36 64; 216 244];

%% plot
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

%% plots
%%% plot PCR by incorrect trials

d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & PlotData.engaged == 1 & PlotData.rewarded == (z-1) & PlotData.licksinRZ == (z-1);
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        disp(['incorrect trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'incorrect', 'correct', 'incorrect', 'correct'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'Reward related prospective coding', 'incorrect trials'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_correct_incorrect_trials';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'pdf')

%%% plot PCR by unengaged trials

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
filename = 'RRZ_PCR_pertrial_180decoding_correct_enengaged_trials';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')

%%
%%% plot PCR by unengaged trials but only include days that have both
%%% engaged and unengaged trials

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

d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = inclDays & strcmp(PlotData.group, gnames{g});
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
filename = 'RRZ_PCR_pertrial_180decoding_correct_enengaged_trials_bothdays';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')

%% get within day difference
BothDayData = PlotData(inclDays,:); 

days = unique(BothDayData.day); groups = [];
for d = 1:length(days)
    isDay = BothDayData.day == days(d);
    isCorrect = isDay & BothDayData.engaged == 1 & BothDayData.rewarded == 1;
    isUnengaged = isDay & BothDayData.engaged == 0;
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
ylabel('difference')
title('correct-unengaged within day PCR')
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




%%
%%% plot PCR by incorrect trials, only days that have correct and incorrect
nCorrect = []; nIncorrect = []; 
days = unique(PlotData.day);
for d = 1:length(days)
    isDay = PlotData.day == days(d);
    nCorrect = sum(PlotData(isDay,:).engaged & PlotData(isDay,:).rewarded);
    nIncorrect = sum(PlotData(isDay,:).engaged & ~PlotData(isDay,:).rewarded & PlotData(isDay,:).licksinRZ == 0);
    groups{d} = unique(PlotData(isDay,:).group);
    hasBoth(d) = nCorrect > 0 & nIncorrect > 0;
end
hasBothDays = days(hasBoth); 
hasBothGroups = groups(hasBoth);
inclDays = ismember(PlotData.day, hasBothDays);

d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = inclDays & isGroup & PlotData.engaged == 1 & PlotData.rewarded == (z-1) & PlotData.licksinRZ == (z-1);
        vdat.([gnames{g}, 'rewarded', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        disp(['incorrect trials ', gnames{g}, ' z = ', num2str(z), ' # trials = ', num2str(sum(isPlotTrial))])
        
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'incorrect', 'correct', 'incorrect', 'correct'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'Reward related prospective coding', 'incorrect trials'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_correct_incorrect_trials_bothdays';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'png')

%% get within day for incorrect trials
BothDayData = PlotData(inclDays,:); 

days = unique(BothDayData.day); groups = [];
for d = 1:length(days)
    isDay = BothDayData.day == days(d);
    isCorrect = isDay & BothDayData.engaged == 1 & BothDayData.rewarded == 1;
    isIncorrect = isDay & BothDayData.engaged == 1 & BothDayData.rewarded == 0 & BothDayData.licksinRZ == 0;
    nCorrect = sum(isCorrect);
    avgCorrect(d) = mean(BothDayData(isCorrect,:).PCR_loc_trial, 'omitnan');
    nIncorrect = sum(isIncorrect);
    avgIncorrect(d) = mean(BothDayData(isIncorrect,:).PCR_loc_trial, 'omitnan');
    avgDayDiff(d) = avgCorrect(d) - avgIncorrect(d);
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
ylabel('PCR difference')
title('correct-incorrect within day PCR')
makefigurepretty(gcf)
filename = 'perday_difference_correct_incorrect_RRZ_PCR';
savefigALP(figdir, filename, 'filetype', 'png')

fh = figure('Position', [353 328 474 256]);
hold on
xplotvals = [4 5; 1 2];
for g = [2,1]
    isGroup = strcmp(groups, gnames{g});
    scatterdat{1} = avgIncorrect(isGroup);
    scatterdat{2} = avgCorrect(isGroup);
    
    mn = [mean(avgIncorrect(isGroup), 'omitnan') mean(avgCorrect(isGroup), 'omitnan')];
    sem = [std(avgIncorrect(isGroup), 'omitnan')./sqrt(sum(~isnan(avgIncorrect(isGroup)))) std(avgCorrect(isGroup), 'omitnan')./sqrt(sum(~isnan(avgCorrect(isGroup))))];
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
xticklabels({'incorrect', 'correct', 'incorrect', 'correct'})
ylim([-0.15 0.15])
yticks([-0.15 0 0.15])
ylabel('prospective coding ratio')
title('day avg reward related zone PCR')
makefigurepretty(gcf)
filename = 'perday_correct_incorrect_RRZ_PCR';
savefigALP(figdir, filename, 'filetype', 'png')

%% plot incorrect trials out of curiosity
BothTrialMetrics = FullTrialMetrics(inclDays,:);
isIncorrectTrial = BothDayData.engaged == 1 & BothDayData.rewarded == 0 & BothDayData.licksinRZ == 0;

IncorrectTrials = BothTrialMetrics(isIncorrectTrial,:); 
plotedges = -81:3:99;
plotedges = plotedges+1.5;

for t = 1:height(IncorrectTrials)
    figure
    hold on
    subplot(2,1,1)
    hold on
    plot(plotedges(1:end-1), IncorrectTrials(t,:).vel_h)
    subplot(2,1,2)
    hold on
    plot(plotedges(1:end-1), IncorrectTrials(t,:).lick_nh)
end










