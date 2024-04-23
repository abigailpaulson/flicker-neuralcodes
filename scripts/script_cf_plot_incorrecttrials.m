%%% Plot # or prop of incorrect trials per day
%
%ALP 3/12/2024

clear; close all
%%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%%% load behavior data
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\'];
filename = ['behavioranalysis_prepost_recordings_table_v2.mat'];
load([filedir, filename])

%%% get trials to include
inclData = allTMetrics.fullTrial==1 & strcmp(allTMetrics.timepoint, 'post');
PlotData = allTMetrics(inclData,:);

%%% helpful stuff for plots
groupcolors = cbrewer('qual', 'Paired', 6); 
gcolors = groupcolors([2,1,4,3],:); %rearranging it cuz of the way that the gramm does the plotting
scolors = groupcolors;
gnames = {'gamma', 'random'};
dnames = {'pre', 'post'}; %here post is used for the color for fast trials

%%% get proportion / number of incorrect trials per day
xplotvals = [1 2];
plotday = 2;
colororder = [];
iBar = 1;
for g = 1:2 
    isGroup = strcmp(PlotData.group, gnames{g});
    GroupData = PlotData(isGroup,:);
    tmpDays = unique(GroupData.day);
    for d = 1:length(tmpDays)
        isDay = GroupData.day == tmpDays(d);
        isEng = GroupData.engaged == 1;
        isRew = GroupData.rewarded == 1;
        isFull = GroupData.fullTrial == 1;
        
        nWrong(d) = sum(isDay&isEng&~isRew&isFull);
        nUnengaged(d) = sum(isDay&~isEng&~isRew&isFull);
        nCorrect(d) = sum(isDay&isEng&isRew&isFull);
        nTot(d) = sum(isDay&isFull);
        nTotCorr(d) = sum(isDay&isFull&isEng);
    end
    
    scatterdat_wrong{iBar} = nWrong;
    scatterdat_correct{iBar} = nCorrect;
    scatterdat_unengaged{iBar} = nUnengaged;
    scatterdat_propwrong{iBar} = nWrong./nTot;
    scatterdat_propcorrect{iBar} = nCorrect./nTotCorr;
    scatterdat_propuneng{iBar} = nUnengaged./nTot;
    
    mn_wrong(g) = mean(nWrong, 'omitnan');
    mn_unengaged(g) = mean(nUnengaged, 'omitnan');
    mn_correct(g) = mean(nCorrect, 'omitnan');
    mn_propwrong(g) = mean(nWrong./nTot, 'omitnan');
    mn_propuneng(g) = mean(nUnengaged./nTot, 'omitnan');
    mn_propcorrect(g) = mean(nCorrect./nTotCorr, 'omitnan');
    
    stde_wrong(g) = std(nWrong, 'omitnan')./sqrt(sum(~isnan(nWrong)));
    stde_unengaged(g) = std(nUnengaged, 'omitnan')./sqrt(sum(~isnan(nUnengaged)));
    stde_correct(g) = std(nCorrect, 'omitnan')./sqrt(sum(~isnan(nCorrect)));
    stde_propwrong(g) = std(nWrong./nTot, 'omitnan')./sqrt(sum(~isnan(nWrong)));
    stde_propuneng(g) = std(nUnengaged./nTot, 'omitnan')./sqrt(sum(~isnan(nUnengaged)));
    stde_propcorrect(g) = std(nCorrect./nTotCorr, 'omitnan')./sqrt(sum(~isnan(nCorrect)));
    
    
    nWrong = [];
    nUnengaged = [];
    nTot = [];
    nCorrect = [];
    iBar = iBar+1;
    colororder = [colororder; params.colors.(gnames{g}).(dnames{plotday})];
end

figdir = filedir; 
fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xplotvals, scatterdat_wrong, 'color', colororder)
b = bar(xplotvals, mn_wrong, 'FaceColor', 'flat');
er = errorbar(xplotvals, mn_wrong, stde_wrong, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
title('# wrong trials per day')
ylabel('number of trials')
title('wrong trials')
xticks([1 2])
xticklabels({'40 Hz', 'Random'})
makefigurepretty(gcf)
figname = 'count_wrongtrials';
savefigALP(figdir, figname, 'filetype', 'png')

fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xplotvals, scatterdat_unengaged, 'color', colororder)
b = bar(xplotvals, mn_unengaged, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_unengaged, stde_unengaged, 'Color','k');
er.LineStyle = 'none';
title('unengaged trials')
ylabel('number of trials')
xticks([1 2])
xticklabels({'40 Hz', 'Random'})
makefigurepretty(gcf)
figname = 'count_unengagedtrials';
savefigALP(figdir, figname, 'filetype', 'png')

fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xplotvals, scatterdat_propcorrect, 'color', colororder)
b = bar(xplotvals, mn_propcorrect, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_propcorrect, stde_propcorrect, 'Color','k');
er.LineStyle = 'none';
title('prop correct trials')
ylabel('proportion of trials')
xticks([1 2])
xticklabels({'40 Hz', 'Random'})
makefigurepretty(gcf)
figname = 'prop_correcttrials';
savefigALP(figdir, figname, 'filetype', 'png')

fh = figure('Position', [165 480 257 219]);
hold on
plotprettypoints(fh, xplotvals, scatterdat_correct, 'color', colororder)
b = bar(xplotvals, mn_correct, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_correct, stde_correct, 'Color','k');
er.LineStyle = 'none';
title('# correct trials')
ylabel('number of trials')
xticks([1 2])
xticklabels({'40 Hz', 'Random'})
makefigurepretty(gcf)
figname = 'count_correcttrials';
savefigALP(figdir, figname, 'filetype', 'png')


fh = figure;
hold on
plotprettypoints(fh, xplotvals, scatterdat_propwrong, 'color', colororder)
b = bar(xplotvals, mn_propwrong, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_propwrong, stde_propwrong, 'Color','k');
er.LineStyle = 'none';
title('proportion wrong trials per day')

fh = figure;
hold on
plotprettypoints(fh, xplotvals, scatterdat_propuneng, 'color', colororder)
b = bar(xplotvals, mn_propuneng, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
er = errorbar(xplotvals, mn_propuneng, stde_propuneng, 'Color','k');
er.LineStyle = 'none';
title('proportion unengaged trials per day')



