function cf_plot_rippledecoding_RZripples(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
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
inclRipples = AllData.RZRipple & AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos);
PlotData = AllData(isPost&inclRipples,:); 

isPostDay = strcmp(DayData.timepoint, 'post');
DayData = DayData(isPostDay,:);

nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

%%% plot ripple properties 
xvect = [1:2; 4:5];

%% plot prospective coding ratio for all ripples
fnames = {'rippleRatio'}; 
yy = {'prospective coding ratio'};
analysisname = {'ripple prospective coding ratio'};

figure(fh)
axes(ax{iPlot})
axis square
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

%% plot distribution of relative decoded positions for all ripples
err_edges = -90:18:90; 
plot_err_edges = err_edges(1:end-1)+9;
axes(ax{iPlot})
hold on
d = 2;
for g = 1:length(gnames)
    tmpdat = [];
    isGroup = strcmp(PlotData.group, gnames{g});
    tmpdat = PlotData.rippleRelPos(isGroup);
    %correct x axis as needed
    tmpdat(tmpdat > 90) = tmpdat(tmpdat > 90) - 180;
    tmpdat(tmpdat < -90) = tmpdat(tmpdat < -90) + 180;
    
    %get normalized histogram
    err_h = histcounts(tmpdat, err_edges);
    err_nh = err_h./sum(err_h);
    
    plot(plot_err_edges, err_nh, 'Color', groupcolors.(gnames{g}).(dnames{d}), 'LineWidth', 2)
end
xlim([-90 90])
xticks([-90 0 90])
ylim([0 0.3])
yticks([0 0.1 0.2 0.3])
ylabel('Fraction of ripples')
xlabel('Relative decoded position (deg)')
iPlot = iPlot + 1;

%% plot per day prospective coding ratio for all ripples
fnames = {'prospRipples_RZ'};
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
% cf_stats2txt2(scatterdat, statfid, panelL{iPlot}, 'days', 'fraction of ripples', 'prospective ripples', gnames, tablefilename);

ylabel('fraction of ripples')
xticks([])
title('prospective ripples')

end

