function cf_plot_360decoding_unoccpiedzone_ripples(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_360decoding_unoccpiedzone_ripples
%
%ALP 4/16/24

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
%from script_cf_plot_SWRdecoding_360_otherzone
savedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\';
filename = 'RippleDecoding_360Decoding_otherzone_pertrial_reviewerfiguredata.mat';
load([savedir, filename])

PlotData = RippleDecodingUnoccupied;

%%% change below when I add A
iPlot = 1;
iStat = 1;

figure(fh)

%%% proportion of trials that have unoccupied zone decoding
axes(ax{iPlot})
xvals = [1 2]; colororder = [];
for g = 1:2
    isGroup = strcmp({RippleData.group}, gnames{g})';
    scatterdat{g} = [RippleData(isGroup).propSigRZ];
    mn(g) = mean([RippleData(isGroup).propSigRZ], 'omitnan');
    sem(g) = std([RippleData(isGroup).propSigRZ], 'omitnan')./sqrt(sum(~isnan([RippleData(isGroup).propSigRZ])));
    colororder = [colororder; params.colors.(gnames{g}).post];
end
plotprettypoints(fh, xvals, scatterdat, 'color', colororder)
b = bar(xvals, mn, 'FaceColor', 'flat');
er = errorbar(xvals, mn, sem, 'Color','k');
b.CData = colororder;
b.FaceAlpha = 0.6;
er.LineStyle = 'none';
cf_stats2txt2(scatterdat, statfid, panelL{iStat}, 'days', 'proportion of ripples', ...
    'prop ripples with unoccupied zone decoding', gnames, tablefilename)
title({'RZ ripples', 'non-local zone decoding'})
ylabel('proportion of ripples')
xticks([1 2])
xticklabels(gnames)

iStat = iStat+1;
iPlot = iPlot+1;

%%% speed on trials where >=1 ripple has unoccupied zone decoding. only
%%% trials with > 1 ripple at all are included.
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & [PlotData.isNonLocal] == (z-1);
        else
            isPlotTrial = isGroup & [PlotData.isNonLocal] >= (z-1);
        end
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.Ctrl_speed(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.Ctrl_speed(isPlotTrial);
    end
    
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'deg/s', ...
        'ctrl speed', {[gnames{g} ' local'], [gnames{g}, ' nonlocal']}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Speed (deg/s)')
title('ctrl speed')
iPlot = iPlot+1;

%%% LICK LATENCY --- on trials where >=1 ripple has unoccupied zone decoding. only
%%% trials with > 1 ripple at all are included.
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & [PlotData.isNonLocal] == (z-1);
        else
            isPlotTrial = isGroup & [PlotData.isNonLocal] >= (z-1);
        end
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.licklatency_s(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.licklatency_s(isPlotTrial);
        
        disp([gnames{g}, 'z = ', num2str(z), ' - lick latency there are ', num2str(sum(datforstats{z} > 2)), ' pts above 1.5'])
    end
    
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 's', ...
        'lick latency', {[gnames{g} ' local'], [gnames{g}, ' nonlocal']}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick latency (s)')
ylim([0 1.5])
iPlot = iPlot+1;

%%% LICK DISCRIMINATION --- on trials where >=1 ripple has unoccupied zone decoding. only
%%% trials with > 1 ripple at all are included.
axes(ax{iPlot})
vdat = []; cmat = []; datforstats = [];
for g = [1,2]
    for z = 1:2
        isGroup = strcmp(PlotData.group, gnames{g});
        if z == 1
            isPlotTrial = isGroup & [PlotData.isNonLocal] == (z-1);
        else
            isPlotTrial = isGroup & [PlotData.isNonLocal] >= (z-1);
        end
        vdat.([gnames{g}, '_', num2str(z-1)]) = PlotData.lickDI(isPlotTrial);
        cmat = [cmat; params.colors.(gnames{g}).(dnames{z})];
        datforstats{z} = PlotData.lickDI(isPlotTrial);
    end
    
    cf_stats2txt2(datforstats, statfid, panelL{iStat}, 'trials', 'lick DI', ...
        'lick DI', {[gnames{g} ' local'], [gnames{g}, ' nonlocal']}, tablefilename)
    iStat = iStat+1;
    datforstats = [];
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'ShowData', false, 'BoxWidth', 0.018, 'MedianSize', 40, 'ViolinAlpha', 0.4);
xticklabels({'local', 'non-local', 'local', 'non-local'})
ylabel('Lick discrimination index')
iPlot = iPlot+1;



end

