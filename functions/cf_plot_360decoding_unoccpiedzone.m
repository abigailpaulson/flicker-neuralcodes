function cf_plot_360decoding_unoccpiedzone(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_360decoding_unoccpiedzone
%
%ALP 4/16/24

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
filename = 'ThetaSeq_360decoding_trialdata_PCR_otherzonedecoding.mat';
load([seqdir, filename])

PlotData = OtherZonesTrial;
PlotDataDay = struct2table(OtherZonesDay);
PlotDataDaySpeed = OtherZonesDaySpeed;

%%% change below when I add A
iPlot = 1;
iStat = 1;

%% get trial speed from unengaged vs. correct trials
figure(fh)

axes(ax{iPlot})
xplotvals = [1 2; 4 5];
for g = 1:2
    isGroup = strcmp(PlotDataDay.group, gnames{g});
  
    tmpDat = PlotDataDay.propNonLocal_RZ(isGroup);
    scatterdat{2} = tmpDat;
    mn(2) = mean(tmpDat, 'omitnan');
    sem(2) = std(tmpDat, 'omitnan')./sqrt(sum(~isnan(tmpDat)));
    
    tmpDat = PlotDataDay.propNonLocal_Ctrl(isGroup);
    scatterdat{1} = tmpDat;
    mn(1) = mean(tmpDat, 'omitnan');
    sem(1) = std(tmpDat, 'omitnan')./sqrt(sum(~isnan(tmpDat)));
    
    CData = [params.colors.(gnames{g}).(dnames{1}); params.colors.(gnames{g}).(dnames{2})];
    plotprettypoints(fh, xplotvals(g,:), scatterdat, 'color', CData)
    b = bar(xplotvals(g,:), mn, 'FaceColor', 'flat');
    b.FaceAlpha = 0.6;
    er = errorbar(xplotvals(g,:), mn, sem, 'Color','k');
    er.LineStyle = 'none';
    b.CData = CData;
    
    cf_stats2txt2(scatterdat, statfid, panelL{iStat}, 'days', 'proportion', ...
        'proportion of trials with decoding of zone', {'RZ', 'Ctrl'}, tablefilename)
    iStat = iStat+1;
    
    mn = []; sem = []; tmpDat = []; scatterdat = [];
end
xticks([1 2 4 5])
ylim([0 1])
xticklabels({'Ctrl', 'RZ', 'Ctrl', 'RZ'})
ylabel('Proportion of trials')
iPlot = iPlot + 1;

%%% lick latency
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.licklatency_s(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.licklatency_s(isPlotTrial);
        disp([gnames{g} ' - There are ', num2str(length(unique(PlotData.animal(isPlotTrial)))), ' animals contributing to unoccupied zone per day analysis; z = ', num2str(z)])
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'lick latency', 's', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
ylim([0 2])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick latency (s)')
title('Lick latency')
iPlot = iPlot+1;

%%% lick DI
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.LickDI(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.LickDI(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'lick discrimination', 'lick DI', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick discrimination index')
title('lick DI')
iPlot = iPlot+1;

%%% trial speed
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        isPlotTrial = isGroup & [PlotData.sig_otherpos_RRZ] == (z-1);
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.Ctrl_speed(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.Ctrl_speed(isPlotTrial);
    end
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'trial speed', 'deg/s', {'local', 'nonlocal'}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
% ylim([-0.5 0.75])
% yticks([-0.5 0 0.5])
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Speed (deg/s)')
title('trial speed')
iPlot = iPlot+1;


%% per day information
axes(ax{iPlot})
xplotvals = [1 2; 4 5]; CData = []; tmpDat = [];
for g = [1,2]
    for s = 1:2
        isGroup = strcmp(PlotDataDaySpeed.group, gnames{g});
        isSpeed = PlotDataDaySpeed.speedType == s;
        
        tmpDat = PlotDataDaySpeed.propNonLocal(isGroup&isSpeed);
        scatterdat{s} = tmpDat;
        mn = mean(tmpDat, 'omitnan');
        sem = std(tmpDat, 'omitnan')./sqrt(sum(~isnan(tmpDat)));
        
        plotprettypoints(fh, xplotvals(g,s), scatterdat(s), 'color', params.colors.(gnames{g}).(dnames{s}))
        b = bar(xplotvals(g,s), mn, 'FaceColor', 'flat');
        b.FaceAlpha = 0.6;
        er = errorbar(xplotvals(g,s), mn, sem, 'Color','k');
        er.LineStyle = 'none';
        b.CData = params.colors.(gnames{g}).(dnames{s});
        %CData = [CData; params.colors.(gnames{g}).(dnames{s})];
        
        mn = []; sem = []; tmpDat = [];
    end
    
    cf_stats2txt2(scatterdat, statfid, panelL{iStat}, 'days', 'proportion of trials', ...
        'proportion of trials with unoccupied zone decoding', {[gnames{g}, ' slow'], [gnames{g}, ' fast']}, tablefilename)
    disp([gnames{g} ' - There are ', num2str(length(unique(PlotDataDaySpeed.animal(isGroup)))), ' animals contributing to unoccupied zone per day analysis'])

    scatterdat = [];
    iStat = iStat + 1;
end
% b.CData = CData;
xticks([1 2 4 5])
xticklabels({'slow', 'fast', 'slow', 'fast'})
ylabel('Proportion of trials')
title('significant other zone decoding by trial speed')







end

