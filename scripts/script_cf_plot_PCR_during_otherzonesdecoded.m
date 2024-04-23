clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table_v2.mat'];
load([filedir, filename])

%%% load theta seq decoding data
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
seqfilename = ['thetaseqdecoding_alldays_alltrials_360dec_180trialsplit'];
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

%%% set up data to only include good trials
%%% get trials to include
isNaNTrial = isnan(AllData.trialNum);
AllData = AllData(~isNaNTrial,:);
inclData = AllData.fullTrial==1 & AllData.engaged ==1 & AllData.rewarded == 1 & strcmp(AllData.timepoint, 'post') & AllData.significantSeq == 1;
PlotData = AllData(inclData,:);

%%% do the same set up for behavior metrics
BehavData = allTMetrics(~isNaNTrial,:);
BehavData = BehavData(inclData,:); 

%% plot section
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\decoding_360_withsigotherzones\';


%%% plot --- first check if PCR in the RRZ is different between 40 Hz and Random
%%% with 360 decoding, both reward zones combined
d = 2; %cuz post only
figure('Position', [409 341 278 278])
hold on
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
iPlot = 1;
vdat = []; cmat = [];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    isPlotTrial = isGroup;
    vdat.([gnames{g}]) = PlotData.PCR_RRZ_trial(isPlotTrial);
    cmat = [cmat; params.colors.(gnames{g}).(dnames{d})];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 1.5])
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'360 decoding', 'reward related zone prospective coding'})
makefigurepretty(gcf,1)
filename = 'RRZ_PCR_pertrial_360decoding';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% plot --- first check if PCR in the RRZ is different between 40 Hz and Random
%%% with 360 decoding, split by reward zones
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
plot([0 10], [0 0], 'k-', 'LineWidth', 0.5)
vdat = []; cmat = [];
for z = 1:2
    for g = [2,1]
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [BehavData.zone] == z;
        vdat.([gnames{g}, '_zone', num2str(z)]) = PlotData.PCR_RRZ_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([-0.75 0.75])
xlim([0.5 2.5])
xticklabels({'zone1', 'zone1', 'zone2', 'zone2'})
yticks([-0.75 0 0.75])
ylabel('prospective coding ratio')
title({'360 decoding', 'reward related zone prospective coding'})
makefigurepretty(gcf,1)
filename = 'RRZ_PCR_pertrial_360decoding_perRZ';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)


%%% plot --- violin plot of PCR for trials with significant other zone
%%% decoding
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.PCR_RRZ_trial(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
ylabel('prospective coding ratio')
title('PCR for trials without and with sig other zone decoding')
makefigurepretty(gcf,1)
filename = 'RRZ_PCR_pertrial_360decoding_bySigOtherZone';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% violin plot
% per group, split by sig or non-sig other RZ decoding
% plot speed, plot lick DI
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.lickDI(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('lick discrimination index')
title('lick DI')
makefigurepretty(gcf)
filename = 'lickDI_pertrial_bySigOtherZone';
savefigALP(figdir, filename, 'filetype', 'png')

%%% violin plot
% per group, split by sig or non-sig other RZ decoding
% plot speed, plot lick DI
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.Ctrl_speed(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('speed (deg/s)')
title('trial speed')
makefigurepretty(gcf)
filename = 'speed_controlzone_pertrial_bySigOtherZone';
savefigALP(figdir, filename, 'filetype', 'png')

%%% violin plot
% per group, split by sig or non-sig other RZ decoding
% plot speed, plot lick DI
d = 2; %cuz post only
figure('Position', [409 402 303 217])
hold on
iPlot = 1;
vdat = []; cmat = [];
for g = [2,1]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.licklatency_s(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([0 2])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('lick latency (s)')
title('lick latency')
makefigurepretty(gcf)
filename = 'licklatency_pertrial_bySigOtherZone';
savefigALP(figdir, filename, 'filetype', 'png')

% %%% violin plot
% % per group, split by sig or non-sig other RZ decoding
% % plot speed, plot lick DI
% d = 2; %cuz post only
% figure('Position', [409 402 303 217])
% hold on
% iPlot = 1;
% vdat = []; cmat = [];
% for g = [2,1]
%     for z = 1:2
%         isGroup = strcmp(PlotData.group, gnames{g});
%         isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
%         vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.postRZ_speed(isPlotTrial);
%         cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
%     end
% end
% violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% % ylim([-0.5 0.75])
% % yticks([-0.5 0 0.5])
% ylabel('speed (deg/s)')
% title('post RZ speed - separated by sig non-local RRZ decoding')
% makefigurepretty(gcf)
% filename = 'speed_postRZ_pertrial_bySigOtherZone';
% savefigALP(figdir, filename, 'filetype', 'png')


%%% violin plot
% per group, split by correct or incorrect trials
% plot PCR in RRZ


%%% plot PCR over position for 360 deg trials
figure('Position', [242 487 362 185])
hold on
patch([params.RRZ(1,1) params.RRZ(1,2) params.RRZ(1,2) params.RRZ(1,1)], [-1 -1 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
patch([params.RRZ(2,1) params.RRZ(2,2) params.RRZ(2,2) params.RRZ(2,1)], [-1 -1 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
iPlot = 1;
for g = 1:2
    isGroup = strcmp(PlotData.group, gnames{g});
    mn = []; sem = [];
    mn = mean([PlotData.PCR_pos_trial(isGroup,:)], 'omitnan');
    sem = std([PlotData.PCR_pos_trial(isGroup,:)], 0,1, 'omitnan') ./sqrt(sum(isGroup)/2); %div by two bc each "trial" here is 1 approach tot he reward zone. going to plot all 360 deg
    shadedErrorBar(params.dec_edges(1:end-1), mn, sem, {'Color', params.colors.(gnames{g}).(dnames{z})},1);
    iPlot = iPlot+1;
end
xlabel('track position (deg)')
ylabel('prospective coding ratio')
ylim([-0.10 0.25])
yticks([-0.10 0 0.1 0.20])
xticks([0 180 360])
xlim([0 360])
makefigurepretty(gcf,1)
filename = 'PCR_pertrial_360decoding_overposition';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% plot PCR over position for 360 deg trials - split by significant other
%%% zone decoding
% plottitles = { 'no other zone decoding','sig other zone decoding'};
% figure('Position', [196 409 990 224])
% hold on
% for z = 1:2
%     subplot(1,2,z)
%     hold on
%     patch([params.RRZ(1,1) params.RRZ(1,2) params.RRZ(1,2) params.RRZ(1,1)], [-1 -1 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
%     patch([params.RRZ(2,1) params.RRZ(2,2) params.RRZ(2,2) params.RRZ(2,1)], [-1 -1 1.4 1.4], 'm', 'FaceAlpha', 0.1, 'EdgeColor', 'none')
%     iPlot = 1;
%     for g = 1:2
%         isGroup = strcmp(PlotData.group, gnames{g});
%         isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
%         mn = []; sem = [];
%         mn = mean([PlotData.PCR_pos_trial(isPlotTrial,:)], 'omitnan');
%         sem = std([PlotData.PCR_pos_trial(isPlotTrial,:)], 0,1, 'omitnan') ./sqrt(sum(isPlotTrial)/2); %div by two bc each "trial" here is 1 approach tot he reward zone. going to plot all 360 deg
%         shadedErrorBar(params.dec_edges(1:end-1), mn, sem, {'Color', params.colors.(gnames{g}).(dnames{z})},1);
%         iPlot = iPlot+1;
%     end
%     xlabel('track position (deg)')
%     ylabel('prospective coding ratio')
%     ylim([-0.10 0.25])
%     yticks([-0.10 0 0.1 0.20])
%     xticks([0 180 360])
%     xlim([0 360])
%     title(plottitles{z})
% end
% makefigurepretty(gcf,1)
% filename = 'PCR_pertrial_360decoding_overposition_bySigOtherZone';
% savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% plot --- get prop of trials/day with significant other zone deocoding by low/med/high RRZ PCR 
pbuckets = [-1.0 -0.1 0 0.1 1];
xplotvals = [1 2 3; 5 6 7];
plotday = 2; %plot only
figure('Position', [327 294 720 302])
hold on
for g = 1:2
    isGroup = strcmp(PlotData.group, gnames{g});
    grpData = PlotData(isGroup,:);
    seqType = discretize(grpData.PCR_RRZ_trial, pbuckets);
    
    days = unique(grpData.day);
    
    for d = 1:length(days)
        isDay = grpData.day == days(d);
        nTrials(d) = sum(isDay);
        for p = 1:3
            isP = isDay & seqType == p;
            nP(d,p) = sum(grpData.sig_otherpos_RRZ(isP));
            propP(d,p) = nP(d,p)./nTrials(d);
        end
    end
   
    mn = mean(propP,1, 'omitnan');
    sem = std(propP,0,1,'omitnan')./sqrt(length(days));
    
    bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{plotday}))
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    days = []; nTrials = []; nP = []; propP = [];
end
xticks([1 2 3 5 6 7])
xticklabels({'low PCR', 'near 0 PCR', 'high PCR', 'low PCR', 'near 0 PCR', 'high PCR'})
ylabel('proportion of trials')
title('prop trials with other zone decoding in RRZ by RRZ PCR subset')
makefigurepretty(gcf,1)
filename = 'prop_sigOtherZone_perday_360decoding_byRRZPCR';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% plot --- get prop of trials/day with significant other zone decoding by
%%% fast/slow trials
fh = figure('Position', [327 294 720 302]);
hold on
xplotvals = [1 2; 4 5];
split = median(BehavData.Ctrl_speed); 
buckets = [min(BehavData.Ctrl_speed) split max(BehavData.Ctrl_speed)];
allColors  = []; forStatsDay = []; forStatsGroup = []; forStatsPropRew = []; forStatsS = []; forStatsAn = [];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    grpData = PlotData(isGroup,:);
    speedType = discretize(grpData.Ctrl_speed, buckets);
    grpAn = [];
    
    days = unique(grpData.day);
    for d = 1:length(days)
        isDay = grpData.day == days(d);
        nTrials(d) = sum(isDay);
        tmpAn = grpData.animal(isDay);
        grpAn(d,1) = tmpAn(1);
        for s = 1:2
            isS = isDay & speedType == s;
            nS(d,s) = sum(grpData.sig_otherpos_RRZ(isS));
            propS(d,s) = nS(d,s)./nTrials(d);
        end
    end
    scatterdat{1} = propS(:,1);
    scatterdat{2} = propS(:,2);
    mn = mean(propS,1, 'omitnan');
    sem = std(propS,0,1,'omitnan')./sqrt(length(days));
    
    CData = [params.colors.(gnames{g}).pre; params.colors.(gnames{g}).post];
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', CData)
    b = bar(xplotvals(g,:), mn, 'FaceColor', 'flat');
    b.FaceAlpha = 0.6;
    b.CData = CData;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    allColors = [allColors; CData];
    
    forStatsDay = [forStatsDay; days; days];
    forStatsGroup = [forStatsGroup; repmat(gnames(g),length(days)*2,1)];
    forStatsPropRew = [forStatsPropRew; propS(:,1); propS(:,2)];
    forStatsS = [forStatsS; ones(length(days),1); 2*ones(length(days),1)];
    forStatsAn = [forStatsAn; grpAn; grpAn];
    
    days = []; nTrials = []; nS = []; propS = [];
end
xticks([1 2 4 5])
xticklabels({'slow', 'fast', 'slow', 'fast'})
ylabel('proportion of trials')
title('prop trials with significant other zone decoding in RRZ by fast/slow trials')
makefigurepretty(gcf,1)
filename = 'prop_sigOtherZone_perday_360decoding_byTrialSpeed';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

% figure('Position', [327 294 720 302])
% hold on
% xplotvals = [1 2; 4 5];
% split = median(BehavData.lickDI); 
% buckets = [min(BehavData.lickDI) split max(BehavData.lickDI)];
% for g = [2,1]
%     isGroup = strcmp(PlotData.group, gnames{g});
%     grpData = PlotData(isGroup,:);
%     speedType = discretize(grpData.Ctrl_speed, buckets);
%     
%     days = unique(grpData.day);
%     for d = 1:length(days)
%         isDay = grpData.day == days(d);
%         nTrials(d) = sum(isDay);
%         for s = 1:2
%             isS = isDay & speedType == s;
%             nS(d,s) = sum(grpData.sig_otherpos_RRZ(isS));
%             propS(d,s) = nS(d,s)./nTrials(d);
%         end
%     end
%     
%     mn = mean(propS,1, 'omitnan');
%     sem = std(propS,0,1,'omitnan')./sqrt(length(days));
%     
%     bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{plotday}))
%     er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
%     er.LineStyle = 'none';
%     
%     days = []; nTrials = []; nS = []; propS = [];
% end
% xticks([1 2 4 5])
% xticklabels({'slow', 'fast', 'slow', 'fast'})
% ylabel('proportion of trials')
% title('prop trials with significant other zone decoding in RRZ by fast/slow trials')
% makefigurepretty(gcf,1)
% filename = 'prop_sigOtherZone_perday_360decoding_byTrialSpeed';
%savefigALP(figdir, filename, 'filetype', 'pdf', 'date', 1)

%%% plot --- get proportion of trials/day with significant reward zone
%%% decoding
nS = []; propS = [];
fh = figure('Position', [353 329 316 255]);
hold on
xplotvals = [2;1];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    grpData = PlotData(isGroup,:);
    
    days = unique(grpData.day);
    for d = 1:length(days)
        isDay = grpData.day == days(d);
        nTrials(d) = sum(isDay);
        nS(d) = sum(grpData.sig_otherpos_RRZ(isDay));
        propS(d) = nS(d)./nTrials(d);
    end
        scatterdat{1} = propS;

    mn = mean(propS, 'omitnan');
    sem = std(propS,'omitnan')./sqrt(length(days));
    
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', params.colors.(gnames{g}).(dnames{plotday}))
    b = bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{plotday}));
    b.FaceAlpha = 0.6;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    days = []; nTrials = []; nS = []; propS = []; scatterdat = [];
end
xticks([1 2])
xticklabels({'random', '40Hz'})
ylim([0 1])
yticks([0 0.2 0.4 0.6 0.8 1])
ylabel('proportion of trials')
title('sig non-local RRZ decoding')
makefigurepretty(gcf)
filename = 'prop_sigOtherZone_perday_RRZ';
savefigALP(figdir, filename, 'filetype', 'png')

%%% plot --- get proportion of trials/day with significant reward zone
%%% decoding
nS = []; propS = [];
fh = figure('Position', [353 329 316 255]);
hold on
xplotvals = [2;1];
for g = [2,1]
    isGroup = strcmp(PlotData.group, gnames{g});
    grpData = PlotData(isGroup,:);
    
    days = unique(grpData.day);
    for d = 1:length(days)
        isDay = grpData.day == days(d);
        nTrials(d) = sum(isDay);
        nS(d) = sum(grpData.sig_otherpos_ctrl(isDay));
        propS(d) = nS(d)./nTrials(d);
    end
    scatterdat{1} = propS;
    mn = mean(propS, 'omitnan');
    sem = std(propS,'omitnan')./sqrt(length(days));
    
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', params.colors.(gnames{g}).(dnames{plotday}))
    b = bar(xplotvals(g,:), mn, 'FaceColor', params.colors.(gnames{g}).(dnames{plotday}));
    b.FaceAlpha = 0.6;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    
    days = []; nTrials = []; nS = []; propS = []; scatterdat = [];
end
xticks([1 2])
ylim([0 1])
xticklabels({'random', '40Hz'})
ylabel('proportion of trials')
yticks([0 0.2 0.4 0.6 0.8])
title('sig non-local control zone decoding')
makefigurepretty(gcf)
filename = 'prop_sigOtherZone_perday_Control';
savefigALP(figdir, filename, 'filetype', 'png')

%% save data for stats

StatsData.animal = PlotData.animal;
StatsData.group = PlotData.group;
StatsData.licklatency_s = PlotData.licklatency_s;
StatsData.sig_otherpos_RRZ = PlotData.sig_otherpos_RRZ;
StatsData.Ctrl_speed = PlotData.Ctrl_speed;
StatsData.LickDI = PlotData.lickDI; 

StatsData = struct2table(StatsData);
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_TrialData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt';
writetable(StatsData, fullfile(statsdir, filename))

StatsDataDay.animal = forStatsAn;
StatsDataDay.day = forStatsDay;
StatsDataDay.group = forStatsGroup;
StatsDataDay.propNonLocal = forStatsPropRew;
StatsDataDay.speedType = forStatsS;

StatsDataDay = struct2table(StatsDataDay);
statsdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\LMM_R\';
filename = 'TableData_DayData_360Decoding_ThetaSeq_PCR_OtherZoneDecoding.txt';
writetable(StatsDataDay, fullfile(statsdir, filename))

%% save data for figures

OtherZonesTrial = StatsData;
OtherZonesDay = StatsDataDay;

seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
filename = 'ThetaSeq_360decoding_trialdata_PCR_otherzonedecoding.mat';
save([seqdir, filename], 'OtherZonesTrial', 'OtherZonesDay')
























