function cf_plot_thetaspikingexample(dirs, params, allindex, metadata, fh, ax)
%cf_plot_thetaspikingexample
%
%ALP 4/18/2023

%%% example day
% exday = 200822;
% exfile = 2;
% extrials = 4:6;
% 
% exday = 200831;
% exfile = 1;
% extrials = 8:10;

exday = 210211;
exfile = 1;
extrials = 14:16;

tempindex = allindex(allindex(:,2) == exday,:);
index = tempindex(tempindex(:,3) == exfile,:);

%%% directories
anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(index(1,1)), '_', ...
    num2str(index(1,2)), '\'];

%%% load spiking data. use old PC data for now but I should update...
[singleunits, clusterMetrics, spikeVect] = cf_getspikes(anprocesseddatadir, dirs, params, ...
    index, params.brainReg);
brainReg = {clusterMetrics.brainReg};

%%% load ephys data
bestchdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\';
load([bestchdir, 'CA1\bestChannel_', num2str(index(1,2)), '.mat'])
ca1ch = bestRippleChan.all; 
ca3ch = bestRippleChan.all; 

ca1 = load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\eeg', num2str(index(3)), '.mat']);
filteredTheta = load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\EEG\theta', num2str(index(3)), '.mat']);
load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\thetas', num2str(index(3)), '.mat'])
load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\ripples', num2str(index(3)), '.mat'])
LFPData = ca1.eeg{index(1)}{index(2)}{index(3)}.data;
LFPFilt = filteredTheta.theta{index(1)}{index(2)}{index(3)}.data(:,1);
LFPtime = 0:1/2000:length(LFPData(1,:))/2000;

%%% get theta periods (should change this to be running times later ALP
%%% 4/18/23
thetaper = [thetas{index(1)}{index(2)}{index(3)}.startind thetas{index(1)}{index(2)}{index(3)}.endind];
thetamid = thetas{index(1)}{index(2)}{index(3)}.midind;
thetaTs = thetas{index(1)}{index(2)}{index(3)}.samprate; 

%%% load rate maps for sorting the spiking
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\spatialmaps_all_500_full_', num2str(index(1,2)), '.mat'])
peakpos = ratemaps.peakpos; 
isPC = ratemaps.includeCell;
PC_peakpos = peakpos(isPC);
PCID = ratemaps.cellIDs(isPC);
[~, iSort] = sort(PC_peakpos);
PCID_sorted = PCID(iSort);

%%% pick some theta times to to plot


figure(fh)
axes(ax{1})
cla
hold on
tTime = 0.5; %inS
timeAroundMid = [thetamid thetamid+2*tTime*thetaTs];
for t = size(thetamid,1) 
    plotTime = timeAroundMid(t,1):timeAroundMid(t,2);
    plotTime = plotTime./thetaTs;
    %%% get spikes in this period
    isPerSpike = logical(isExcluded(spikeVect.spikeTimes, timeAroundMid(t,:)./thetaTs));
    perSpikeTimes = spikeVect.spikeTimes(isPerSpike);
    perSpikeIDs = spikeVect.spikeIDs(isPerSpike);
    isPerPC = ismember(perSpikeIDs, PCID);
    PCperSpikeIDs = perSpikeIDs(isPerPC);
    PCperSpikeTimes = perSpikeTimes(isPerPC);
    
    iPlot = 1;
    for c = 1:length(PCID_sorted)
        isCell = PCperSpikeIDs == PCID_sorted(c);
        nSpikes = sum(isCell);
        if nSpikes == 0
            continue
        end
        yVect = iPlot*ones(1,nSpikes);
        plot_raster(PCperSpikeTimes(isCell), yVect)
        iPlot = iPlot+1;
    end
    LFPscale = 20;
    iPlot = 30;
    plot(plotTime, iPlot+LFPFilt(timeAroundMid(t,1):timeAroundMid(t,2))./LFPscale)
    %scale bar
    plot([plotTime(50) plotTime(50)], [iPlot iPlot+100/LFPscale], 'r') 
    iPlot = 50;
    plot(plotTime, iPlot+LFPData(timeAroundMid(t,1):timeAroundMid(t,2))./LFPscale, 'k')
    %scale bar
    plot([plotTime(50) plotTime(50)], [iPlot iPlot+100/LFPscale], 'r') 
    xlim([plotTime(1) plotTime(end)])
    plot([plotTime(100) plotTime(500)], [50 50], 'r')
    xticks([])
    yticks([])
    title('x scale 200ms y scale 100uV')
end


end

