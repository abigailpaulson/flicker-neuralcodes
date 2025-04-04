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
figure('Position', [443 471 557 327])
hold on
plot(BehavData.nLicksAZ, PlotData.PCR_loc_trial,'k.')
%%% fit a line to the data
incldat = ~isnan(PlotData.PCR_loc_trial) & ~isnan(BehavData.nLicksAZ);
Ydat = PlotData.PCR_loc_trial(incldat);
Xdat = BehavData.nLicksAZ(incldat);
c = polyfit(Xdat, Ydat, 1);
xvals_line = min(BehavData.nLicksAZ):1:max(BehavData.nLicksAZ);
yvals_line = polyval(c,Xdat);
% plot(Xdat, yvals_line, 'r-');
SStot = sum((Ydat-mean(Ydat)).^2);                    % Total Sum-Of-Squares
SSres = sum((Ydat-yvals_line).^2);                       % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                            % R^2
[~,~,~,~, stats] = regress(Ydat,[ones(length(Xdat),1), Xdat]);
title(['R2 = ', num2str(Rsq), ' N trials = ', num2str(sum(incldat))])
xlabel('Number of licks')
ylabel('Prospective coding ratio')
savefigname = 'revieweronly_AZlicknumber_PCR_trial_bothgroups_postonly';
savefigdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\flicker-neuralcodes\manuscriptfigures\';
makefigurepretty(gcf,1);
%savefigALP(savefigdir, savefigname, 'filetype', 'pdf', 'date', 1);

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
    
nAn = unique(PlotData.animal(isGroup));
nDays  = unique(PlotData.day(isGroup));
disp([gnames{g}, ' there are ', num2str(length(nAn)) ' animals and ', num2str(length(nDays)), ' days in the PCR vs. licking analysis'])
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'low licks', 'high licks', 'low licks', 'high licks'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'180 decoding', 'RRZ PCR by Lick #'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_nLickMedian';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
%savefigALP(figdir, filename, 'filetype', 'pdf')

%%% plot PCR by median lick count, compare groups
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
gnames = {'gamma', 'random'};
for z = 1:2
for g = [2,1]

        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & iQ == z;
        vdat.([gnames{g}, 'lick', num2str(z)]) = PlotData.PCR_loc_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
end
    
nAn = unique(PlotData.animal(isGroup));
nDays  = unique(PlotData.day(isGroup));
disp([gnames{g}, ' there are ', num2str(length(nAn)) ' animals and ', num2str(length(nDays)), ' days in the PCR vs. licking analysis'])
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'low licks', 'low licks', 'high licks', 'high licks'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'180 decoding', 'RRZ PCR by Lick #'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_pertrial_180decoding_nLickMedian_comparegroups';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'pdf')

%%% plot PCR by median licking mean and sem across sessions, ALP 7/18/24
d = 2; %cuz post only
xvect_bar = [1 2; 4 5];
gnames  = {'random', 'gamma'};
fh = figure('Unit', 'inches', 'Position', [7.8021 4.8438 3.1458 2.6667]);
daystats.group = []; daystats.lickthresh = []; daystats.animal = []; daystats.PCR = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & iQ == z;
        
        tmpDays = unique(PlotData.day(isPlotTrial));
        avgPCR_sess{g}{z} = []; tmpAnimals = []; scatterdat = [];
        for d = 1:length(tmpDays)
            isDay = PlotData.day == tmpDays(d);
            tmpPCRs = PlotData.PCR_loc_trial(isPlotTrial & isDay);
            avgPCR_sess{g}{z}(d) = mean(tmpPCRs, 'omitnan');
            tmpAnimals(d) = unique(PlotData.animal(isDay));
        end
        scatterdat{z} = avgPCR_sess{g}{z};
        plotdat_mn(z) = mean(avgPCR_sess{g}{z}, 'omitnan');
        plotdat_sem(z) = mean(avgPCR_sess{g}{z}, 'omitnan')./sqrt(length(avgPCR_sess{g}{z}));
        colororder = params.colors.(gnames{g}).(dnames{z});
        
        plotprettypoints(fh, xvect_bar(g,z), scatterdat(z), 'color', colororder)
        b = bar(xvect_bar(g,z), plotdat_mn(z), 'FaceColor', 'flat');
        b.CData = colororder;
        b.FaceAlpha = 0.6;
        errorbar2(xvect_bar(g,z), plotdat_mn(z), plotdat_sem(z), 0.2, 'k-', 'LineWidth', 0.75);
        
        daystats.group = [daystats.group; repmat(gnames(g), length(tmpAnimals),1)];
        daystats.lickthresh = [daystats.lickthresh; repmat(z,length(tmpAnimals),1)];
        daystats.animal = [daystats.animal; tmpAnimals'];
        daystats.PCR = [daystats.PCR; scatterdat{z}'];
        
        clear b
    end   
end
xticks([1,2,4,5])
xticklabels({'low licks', 'high licks', 'low licks', 'high licks'})
% ylim([0 0.4])
% yticks([0 0.1 0.2 0.3 0.4])
ylabel('prospective coding ratio')
title({'avg pcr', 'split by median AZ lick count'})
makefigurepretty(gcf)
filename = 'RRZ_PCR_avg_perday_180decoding_nLickMedian';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
savefigALP(figdir, filename, 'filetype', 'pdf')

%%%% assemble per data data for running stats
DayStats = struct2table(daystats);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_RRZPCR_180_perDay_postDataOnly_LickNumber.txt';
writetable(DayStats, fullfile(statsdir, filename))

%% get PCR by median lick count

%%% get median prospective coding
splitname = 'PCR_loc_trial';
qs = median(PlotData.(splitname), 'omitnan');
qs = [min(PlotData.(splitname)) qs max(PlotData.(splitname))];
iQ = discretize(PlotData.(splitname), qs);

%%% plot lick count by median PCR
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for g = [2,1]
for z = 1:2
        isGroup = strcmp(BehavData.group, gnames{g});
        isPlotTrial = isGroup & iQ == z;
        vdat.([gnames{g}, 'lick', num2str(z)]) = BehavData.nLicksAZ(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
end
    
% nAn = unique(PlotData.animal(isGroup));
% nDays  = unique(PlotData.day(isGroup));
% disp([gnames{g}, ' there are ', num2str(length(nAn)) ' animals and ', num2str(length(nDays)), ' days in the PCR vs. licking analysis'])
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
%ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'low PCR', 'high PCR', 'low PCR', 'high PCR'})
%yticks([-0.75 0 0.75])
ylabel('AZ lick count')
title({'180 decoding', 'Lick # by RRZ PCR'})
%makefigurepretty(gcf)
%filename = 'RRZ_PCR_pertrial_180decoding_nLickMedian';
%figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
%savefigALP(figdir, filename, 'filetype', 'pdf')


%% save the data to run stats

StatsData.group = [PlotData.group];
StatsData.animal = [PlotData.animal];
StatsData.Lick_number_subset = iQ;
StatsData.PCR_loc_trial = [PlotData.PCR_loc_trial];
StatsData = struct2table(StatsData);

statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_RRZPCR_180_perTrial_postDataOnly_LickNumber.txt';
writetable(StatsData, fullfile(statsdir, filename))

%%% save data for figure
PCRAnticipatoryLicking = StatsData; 
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ThetaSeq_TrialData_180Decoding__RRZ_by_AnticipatoryLickMedian_SupplementaryFigure.mat';
save([statsdir, filename], 'PCRAnticipatoryLicking')





