function cf_plot_rippledecoding_allripples(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_rippledecoding_allripples
%
%ALP 3/22/2023

iPlot = 1;

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
datadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\';
filename = 'rippledecoding_halfratiozones_180_PC_AllData.mat'; 
load([datadir, filename])
filename = 'rippledecoding_halfratiozones_180_PC_DayData.mat'; 
load([datadir, filename])

%%% get only the data were interested in
isPost = strcmp(AllData.timepoint, 'post');
inclRipples = AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos);
PlotData = AllData(isPost&inclRipples,:); 

%%% get contribution of diff animals
nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

isPostDay = strcmp(DayData.timepoint, 'post');
DayData = DayData(isPostDay,:);

%%% plot ripple properties 
xvect = [1:2; 4:5];

%% plot per day number of ripples
fnames = {'nRipples', 'nRipples_RZ'};
analysisname = {'ripple count', 'reward zone ripple count'};

axes(ax{iPlot})
hold on
for f = 1:2
    for d = 2
        scatterdat = []; bdat = []; bstde = []; cmat = [];
        for g = 1:length(gnames)
            tmpdat = [];
            isDay = strcmp(DayData.group, gnames{g});
            tmpdat = DayData.(fnames{f})(isDay);
            scatterdat{g} = tmpdat;
            bdat(g) = nanmean(tmpdat);
            bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
        plotprettypoints(fh, xvect(f,:), scatterdat, 'color', cmat)
        b = bar(xvect(f,:), bdat, 'FaceColor', 'flat');
        b.CData = cmat;
        b.FaceAlpha = 0.6;
        errorbar2(xvect(f,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
    end
    cf_stats2txt2(scatterdat, statfid, [panelL{iPlot} num2str(f)], 'days', 'ripples', analysisname{f}, gnames, tablefilename);
end
ylabel('ripple count')
xticks([])
    xticks([])
ylim([0 500])
title('number of Ripples - all/RZ')
iPlot = iPlot+1;

%% plot prospective coding ratio for all ripples
fnames = {'rippleRatio'}; 
yy = {'prospective coding ratio'};
analysisname = {'ripple prospective coding ratio'};

axes(ax{iPlot})
hold on
vdat = []; cmat = []; stat = [];
for d = 2
    datforstats = [];
    for g = 1:2
        nm = [gnames{g}, dnames{d}];
        isGroup = strcmp(PlotData.group, gnames{g});
        tmpdat = PlotData.(fnames{1})(isGroup);
        vdat.(nm) = tmpdat;
        datforstats{g} = tmpdat;
        cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        
        clear tmpdat isGroup
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', false);
ylabel(yy{1})
title([fnames{1}])
cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'ripples', 'prospective coding ratio', analysisname{1}, gnames, tablefilename)

clear datforstats
iPlot = iPlot+1;

%% plot per day prospective coding ratio for all ripples
fnames = {'prospRipples'};
analysisname = {'prospective coding ripples'};

axes(ax{iPlot})
hold on

for d = 2
    scatterdat = []; bdat = []; bstde = []; cmat = [];
    for g = 1:length(gnames)
        tmpdat = [];
        isDay = strcmp(DayData.group, gnames{g});
        tmpdat = DayData.(fnames{1})(isDay);
        scatterdat{g} = tmpdat;
        bdat(g) = nanmean(tmpdat);
        bstde(g) = nanstd(tmpdat)./sqrt(sum(~isnan(tmpdat)));
        cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
    end
    plotprettypoints(fh, xvect(d,:), scatterdat, 'color', cmat)
    b = bar(xvect(d,:), bdat, 'FaceColor', 'flat');
    b.CData = cmat;
    b.FaceAlpha = 0.6;
    errorbar2(xvect(d,:), bdat, bstde, 0.2, 'k-', 'LineWidth', 0.75);
end
cf_stats2txt2(scatterdat, statfid, panelL{iPlot}, 'days', 'fraction of ripples', 'prospective ripples', gnames, tablefilename);

ylabel('fraction of ripples')
xticks([])
title('prospective ripples')

end

