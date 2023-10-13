function cf_plot_rippleproperties(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%cf_plot_rippleproperties
%
%ALP 3/21/23

%% stuff for plotting
groupcolors = params.colors; 
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};

%% load the data structure
datadir = '\\ad.gatech.edu\bme\lbs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\';
filename = 'rippleproperties.mat'; 
load([datadir, filename])

%%% get only the data were interested in
isPost = strcmp(AllSWRData.timepoint, 'post');
isExcludeLong = AllSWRData.rippleDurS < 0.5; 
PlotData = AllSWRData(isPost&isExcludeLong,:); 

%%% get contribution of animals to ripple props
nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

disp(['There are ', num2str(nAnGamma) ' animals in 40 Hz ripple props'])
disp(['There are ', num2str(nAnRandom) ' animals in Random ripple props'])

%%% plot ripple properties 
fnames = {'rippleDurS', 'rippleSize'}; 
yy = {'Ripple duration (s)', 'Ripple power (sd above mean)'};
analysisname = {'ripple duration', 'ripple power'};
yUnit = {'seconds', 's.d.'};
xxl = {[0 0.5], [2 25]};

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
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', true, 'DataPointSize', 2);
    ylabel(yy{f})
    ylim(xxl{f})
        xticks([])
    title([fnames{f}])
    cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'ripples', yUnit{f}, analysisname{f}, gnames, tablefilename)

    clear datforstats
    iPlot = iPlot+1;
end

%% prep the ripple rate data
inclNTData = strcmp(AllNTData.timepoint, 'post') & AllNTData.nonThetaDurS > 5; 
PlotData = AllNTData(inclNTData,:);

%%% ripple abundance
axes(ax{iPlot})
hold on

nAnGamma = length(unique(PlotData.animal(strcmp(PlotData.group, 'gamma'))));
nAnRandom = length(unique(PlotData.animal(strcmp(PlotData.group, 'random'))));

disp(['There are ', num2str(nAnGamma) ' animals in 40 Hz nontheta ripple props'])
disp(['There are ', num2str(nAnRandom) ' animals in Random nontheta ripple props'])

vdat = []; cmat = []; stat = [];
for d = 2
    datforstats = [];
    for g = 1:2
        nm = [gnames{g}, dnames{d}];
        isGroup = strcmp(PlotData.group, gnames{g});
        tmpdat = PlotData(isGroup,:).rippleRatePeriod;
        vdat.(nm) = tmpdat;
        datforstats{g} = tmpdat;
        cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        
        clear tmpdat isGroup
    end
end
violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'BoxWidth', 0.018, 'MedianSize', 40, 'ShowData', false);
ylabel('Ripple abundance (Hz)')
ylim([0 1.5])
    xticks([])
title(['Ripple abundance'])
cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'non-theta periods', 'SWR abundance (Hz)', 'ripple abundance', gnames, tablefilename)

iPlot = iPlot+1;







end

