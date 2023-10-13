function cf_plot_rippleproperties_RZ(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_rippleproperties_RZ
% based on cf_plot_rippleproperties
%ALP 9/28/23

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
datadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\'];
filename = 'rippledecoding_halfratiozones_180_PC_AllData.mat';
load([datadir, filename])

%%% get only the data were interested in
isPost = strcmp(AllData.timepoint, 'post');
isExcludeLong = AllData.rippleDurationS < 0.5; %less than 500ms
isRZripple = ~isnan(AllData.rippleposition) & AllData.rippleposition >=0 & AllData.rippleposition < 18; 
PlotData = AllData(isPost&isExcludeLong&isRZripple,:); 

%%% get contribution of animals to ripple props
nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

disp(['There are ', num2str(nAnGamma) ' animals in 40 Hz ripple props'])
disp(['There are ', num2str(nAnRandom) ' animals in Random ripple props'])

%%% plot ripple properties 
fnames = {'rippleDurationS', 'rippleSize'}; 
yy = {'Ripple duration (s)', 'Ripple power (sd above mean)'};
analysisname = {'ripple duration', 'ripple power'};
yUnit = {'seconds', 's.d.'};
xxl = {[0 0.3], [2 25]};

figure(fh)
iPlot = 1;

for f = 1:length(fnames)
    axes(ax{iPlot})
    hold on
 
    vdat = []; cmat = []; stat = [];
    for d = 2
        datforstats = [];
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            isGroup = strcmp(PlotData.group, gnames{g}); 
            tmpdat = PlotData(isGroup,:).(fnames{f}); 
            vdat.(nm) = tmpdat;
            datforstats{g} = tmpdat;
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
            
            clear tmpdat isGroup
        end
    end
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', false, 'DataPointSize', 2);
    ylabel(yy{f})
    ylim(xxl{f})
        xticks([])
    title([fnames{f}])
    cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'ripples', yUnit{f}, analysisname{f}, gnames, tablefilename)

    clear datforstats
    iPlot = iPlot+1;
end

%% prep the ripple rate data
%%% load the data structure
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\ripplerate_RZ.mat'])

xvect = [1 2];
groupcolors = params.colors; 
isPost = strcmp(RateDataRZ.timepoint, 'post');
isEnough = RateDataRZ.nRipples > 5; 
RateDataRZ = RateDataRZ(isPost&isEnough,:);
PlotData = RateDataRZ; 

fnames = {'rippleRateS'};
analysisname = {'ripple rate'};

nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

disp(['There are ', num2str(nAnGamma) ' animals in 40 Hz RZ ripple rate'])
disp(['There are ', num2str(nAnRandom) ' animals in Random RZ ripple rate'])


figure(fh)
axes(ax{iPlot})
    hold on
hold on
for f = 1
    for d = 2
        scatterdat = []; bdat = []; bstde = []; cmat = [];
        for g = 1:length(gnames)
            tmpdat = [];
            isGroup = strcmp(RateDataRZ.group, gnames{g});
            tmpdat = RateDataRZ.(fnames{f})(isGroup);
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
    cf_stats2txt2(scatterdat, statfid, [panelL{iPlot} num2str(f)], 'days', 'ripple abundance', analysisname{f}, gnames, tablefilename);
end
xticks([])
yticks([0 0.04 0.08])
ylabel('Ripple abundance (Hz)')
title('ripple rate - RZ')








end

