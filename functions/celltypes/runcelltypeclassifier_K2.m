function runcelltypeclassifier_K2(dirs, indices, params, makenewfiles)
%SP 3.14.19
%Updated ALP 1/15/2020 for K2 datastructures

%% initialize variables
savedatadir = [dirs.savedatadir 'celltype\' dirs.processedappend];
dirs.savedfiguresdir = savedatadir; 
spikedatadir = [savedatadir 'cellTypeProps\'];
dirs.spikedatadir = spikedatadir;
if ~exist(savedatadir); mkdir(savedatadir); end
if ~exist(spikedatadir); mkdir(spikedatadir); end
processeddatadir = dirs.processeddatadir;
iden = params.iden;
allindex = indices.allindex;
dayindex = indices.dayindex;

if strcmp(params.brainReg, 'CA1')
    dayindex = dayindex(dayindex(:,2) ~= 210514,:); %missing kilosort files, no spikes
elseif strcmp(params.brainReg, 'CA3')
    dayindex = dayindex(dayindex(:,2) ~= 200928,:); %missing kilosort files, no spikes
    %dayindex = dayindex(dayindex(:,1) < 38,:); %for now ALP 11/18/21
end

makenewfiles = [];
makenewfiles.spikewidth = 0; 
makenewfiles.acorr = 0; 
makenewfiles.all = 1;


%% make filter for eeg waveforms
samprate = 20000;
hpFilt = designfilt('highpassiir', 'StopbandFrequency', 100, ...
     'PassbandFrequency', 500, 'StopbandAttenuation', 60, ...
     'PassbandRipple', 1, 'SampleRate', samprate, 'DesignMethod', 'butter');
 
%% load up clusters
for i = 1:size(dayindex,1)    
    clusterdatadir = fullfile(processeddatadir, [iden num2str(dayindex(i,1)) '_' num2str(dayindex(i,2))], dirs.processedappend, dirs.clusterappend);
    recs = allindex(allindex(:,2) == dayindex(i,2),3);
    [clusterInfo(i)] = getClusterInfo_K2(clusterdatadir, dayindex(i,:), recs);
end

%% get FR data
ALL_meanFR = [];
for i = 1:size(dayindex,1)
    allrecProps(i).dayindex = dayindex(i,:);
    allrecProps(i).meanFR = clustermetrics2mat(clusterInfo(i).metrics, 'stable', 'meanFR');
    allrecProps(i).peakFR = clustermetrics2mat(clusterInfo(i).metrics, 'stable', 'peakFR');
    allrecProps(i).stabletimes = clustermetrics2mat(clusterInfo(i).metrics, 'stable', 'times');
    allrecProps(i).waveforms = clustermetrics2mat(clusterInfo(i).metrics, 'WF', 'mn');
    
    %plot FR data for each recording day
    figure; clf; hold on;
    edges = 0:0.5:30;
    histogram(allrecProps(i).meanFR, edges)
    ylabel('Number of Clusters'); xlabel('Mean FR');
    title(['FR distribution - ', iden,  num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
    ALL_meanFR = [ALL_meanFR allrecProps(i).meanFR'];
    filename = 'distFR_mean_';
    savefigSP(dayindex(i,:), [savedatadir, 'FRDistributions\'], filename, iden)
end

%plot the figure of all recs
figure; clf; hold on;
edges = 0:0.5:30;
histogram(ALL_meanFR, edges); hold on;
ylabel('Number of Clusters'); xlabel('Mean Firing Rate')
title(['Mean Firing Rate Dist - ALL RECS']);
filename = [savedatadir 'FRDistributions\ALLRECS_meanFR_dist'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%% get the spikewidth data
ALL_spikewidth = []; ALL_spikewidth2 = []; 
for i = 1:size(dayindex,1)
    [allrecProps(i).spikewidth, badWFs{i}] = getspikewidthforcelltypeclassification_K2(dayindex(i,:),...
        clusterInfo(i).metrics, spikedatadir, dirs.savedfiguresdir, iden, makenewfiles.spikewidth);
    
    %plot peak2trough for each recording day
    figure; clf; hold on;
    edges = 0:0.01:1.5;
    histogram(allrecProps(i).spikewidth.peak2troughDiff, edges)
    %histogram(allrecProps(i).spikewidth.SW2, edges)
    ylabel('Number of Clusters'); xlabel('Spikewidth')
    title(['Spikewidth Dist - ', iden,  num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
    ALL_spikewidth = [ALL_spikewidth allrecProps(i).spikewidth.peak2troughDiff];
    ALL_spikewidth2 = [ALL_spikewidth2 allrecProps(i).spikewidth.SW2];
    filename = 'peak2troughDiff_dist_';
    savefigSP(dayindex(i,:), [savedatadir, 'SpikewidthDistributions\'], filename, iden)
end

figure; clf; hold on;
edges = 0:0.025:1.5;
histogram(ALL_spikewidth, edges); hold on;
plot([0.5-0.0125 0.5-0.0125],[0 40], 'k--')
ylabel('Number of Clusters'); xlabel('Spikewidth')
title(['Spikewidth Dist - ALL RECS']);
filename = [savedatadir 'SpikewidthDistributions\ALLRECS_peak2troughDiff_dist'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%% get the autocorrelogram data
ALL_meanAC = [];
for i = 1:size(dayindex,1)
    recs = allindex(allindex(:,2) == dayindex(i,2),3);
    [allrecProps(i).autocorr_mean, allrecProps(i).autocorr]...
        = getautocorrelogramforcelltypeclassification_K2(dayindex(i,:),...
        dirs, allrecProps(i).spikewidth.fullWFDiff, clusterInfo(i).metrics, iden, makenewfiles.acorr);
    
    %get the AC distribution
    figure; clf; hold on;
    edges = 0:0.2:10;
    histogram(cell2mat(allrecProps(i).autocorr_mean), edges);
    ylabel('Number of Clusters'); xlabel('Autocorrelogram Mean')
    title(['Autocorrelogram distribution - ', iden,  num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
    ALL_meanAC = [ALL_meanAC cell2mat(allrecProps(i).autocorr_mean)];
    filename = 'distAC_';
    savefigSP(dayindex(i,:), [savedatadir, 'AutocorrelogramDistributions\'], filename, iden)
end

figure; clf; hold on;
edges = 0:0.2:10;
histogram(ALL_meanAC, edges); hold on;
ylabel('Number of Clusters'); xlabel('Mean of the AutoCorrelogram')
title(['Mean AC Dist - ALL RECS']);
filename = [savedatadir 'AutocorrelogramDistributions\ALLRECS_meanAC_dist'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%% plot the 3 parameters
for i = 1:size(dayindex,1)
    plotparametersforcelltypeclassification_K2(dayindex(i,:), savedatadir, iden, allrecProps(i));
end
close all

%% plot strange cells for examination
% ALP 1/22/2021
% plot cells in the lower left quadrant for examination
%
binsize = 5;
for i = 1:size(dayindex,1)
    isShortAC = cell2mat(allrecProps(i).autocorr_mean) < params.thresholds.acthresh_IN;
    isSmallWF = allrecProps(i).spikewidth.peak2troughDiff < params.thresholds.widththresh_IN;
    isCell2Plot = isShortAC & isSmallWF;
    
    for unit = 1:length(isCell2Plot)
        if isCell2Plot(unit)
            wf = allrecProps(i).spikewidth.fullWFDiff(unit,:); 
            peakIdx = allrecProps(i).spikewidth.peakIdxDiff(unit); 
            troughIdx = allrecProps(i).spikewidth.troughIdxDiff(unit);

            acorr = cell2mat(allrecProps(i).autocorr(unit));
            acorr_mean = cell2mat(allrecProps(i).autocorr_mean(unit));
            time = (binsize*(length(acorr)-1))/2;

            figure; hold on;
            subplot(1,2,1); hold on;
            plot(wf)
            plot(peakIdx, wf(peakIdx), 'rs')
            plot(troughIdx, wf(troughIdx), 'ks')
            subplot(1,2,2); hold on;
            bar(-time:binsize:time, acorr)
            plot(acorr_mean*binsize, 2, '*g')
            xlim([-time time])
            sgtitle(['lower L quad cell - ', num2str(dayindex(i,2)), ' unitIdx ', num2str(unit)])
        end
    end
    
end
% pause

%% separate the groups and then plot the traces
%individual waveforms
for i = 1:size(dayindex,1)
    allrecProps(i).cellTypeInfo = plotclassifiedwaveforms_K2(dayindex(i,:),...
        spikedatadir,savedatadir,allrecProps(i),params.thresholds,iden);  
end

%plot average waveforms for all recordings
ALL_PYRwfs = [];
ALL_INwfs = [];
for i = 1:size(dayindex,1)
    %get aligned WFs
    PYRwf = allrecProps(i).cellTypeInfo.alignedWF_PYR;
    INwf = allrecProps(i).cellTypeInfo.alignedWF_IN;    
    ALL_PYRwfs = [ALL_PYRwfs; PYRwf];
    ALL_INwfs = [ALL_INwfs; INwf];  
end      

%get mean and SEM
meanINwf = nanmean(ALL_INwfs,1); % mean across rows
meanPYRwf = nanmean(ALL_PYRwfs,1); %mean across rows
semINwf = nanstd(ALL_INwfs,0,1);%/sqrt(size(ALL_INwfs,1)); %mean across rows
semPYRwf = nanstd(ALL_PYRwfs,0,1);%/sqrt(size(ALL_PYRwfs,1)); %mean across rows

figure; hold on;
time = -1.2375:(1/160):2.5125-(1/160); %need to update this this is not going to work for multiple sampling rates 
%time = 0:0.125:75-0.125; %160 samples = 1ms, 199 sample = 0 point so sample1 is -1.2375ms from start, end sample is 2.5125s, step size is 1/160 or 0.0063
plot(time,meanINwf, 'b-','LineWidth', 3)
plot(time,meanPYRwf, 'r-','LineWidth', 3)
plot(time, meanINwf-semINwf, 'b--', 'LineWidth', 1.5)
plot(time, meanINwf+semINwf, 'b--', 'LineWidth', 1.5)


plot(time, meanPYRwf-semPYRwf, 'r--', 'LineWidth', 1.5)
plot(time, meanPYRwf+semPYRwf, 'r--', 'LineWidth', 1.5)
% ciplot(meanINwf-semINwf, meanINwf+semINwf, time, 'r')
% % alpha(0.3)
% ciplot(meanPYRwf-semPYRwf, meanPYRwf+semPYRwf, time, 'b')
% % alpha(0.3)
title(['Average waveforms']); legend('Interneurons', 'Pyramidal')
ylim([-1 1]);
filename = [savedatadir 'Classification\ALLRECS_avgwaveforms'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%% disp info about # cells
disp([num2str(size(ALL_INwfs,1)), ' interneurons classified'])
disp([num2str(size(ALL_PYRwfs,1)), ' pyramidal cells classified'])

%% plot all the recordings at once
ALL_ac_mean = []; ALL_spikewidth = []; ALL_FR_mean = [];
for i = 1:size(dayindex,1)
    ac_mean = cell2mat(allrecProps(i).autocorr_mean);
    spikewidth = allrecProps(i).spikewidth.peak2troughDiff;
    FR_mean = allrecProps(i).meanFR;
    if size(FR_mean,1) > 1
        FR_mean = FR_mean';
    end
    
    ALL_ac_mean = [ALL_ac_mean ac_mean];
    ALL_spikewidth = [ALL_spikewidth spikewidth];
    ALL_FR_mean = [ALL_FR_mean FR_mean];
end
figure; clf; hold on;
plot3(ALL_ac_mean, ALL_spikewidth, ALL_FR_mean, 'o')
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title(['ALL RECS - Mean AC vs. Peak2Trough Diff  vs. Mean FR'])

%save figures
filename = [savedatadir 'Spikewidth_MeanFR_MeanAC\ALLRECS_peak2troughDiff'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');

%separate the plots by color for different cell types
ALL_ac_mean_PYR = []; ALL_ac_mean_IN = [];
ALL_spikewidth_PYR = []; ALL_spikewidth_IN = [];
ALL_FR_mean_PYR = []; ALL_FR_mean_IN = [];
for i = 1:size(dayindex,1)
    INIdx = allrecProps(i).cellTypeInfo.INidx;
    PYRIdx = allrecProps(i).cellTypeInfo.PYRidx;

    ac_mean = cell2mat(allrecProps(i).autocorr_mean);
    spikewidth = allrecProps(i).spikewidth.peak2troughDiff;
    FR_mean = allrecProps(i).meanFR;
    if size(FR_mean,1) > 1
        FR_mean = FR_mean';
    end
    
    ALL_ac_mean_IN = [ALL_ac_mean_IN ac_mean(INIdx)];
    ALL_spikewidth_IN = [ALL_spikewidth_IN spikewidth(INIdx)];
    ALL_FR_mean_IN = [ALL_FR_mean_IN FR_mean(INIdx)];
    
    ALL_ac_mean_PYR = [ALL_ac_mean_PYR ac_mean(PYRIdx)];
    ALL_spikewidth_PYR = [ALL_spikewidth_PYR spikewidth(PYRIdx)];
    ALL_FR_mean_PYR = [ALL_FR_mean_PYR FR_mean(PYRIdx)];
end
figure; clf; hold on;
plot3(ALL_ac_mean, ALL_spikewidth, ALL_FR_mean, 'ko')
plot3(ALL_ac_mean_PYR, ALL_spikewidth_PYR, ALL_FR_mean_PYR, 'ro'); hold on;
plot3(ALL_ac_mean_IN, ALL_spikewidth_IN, ALL_FR_mean_IN, 'bo'); hold on;
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title(['ALL RECS CLASSIFIED - Mean AC vs. Peak2Trough Diff  vs. Mean FR'])

%save figures
filename = [savedatadir 'Spikewidth_MeanFR_MeanAC\ALLRECS_classifiedtypes'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');

%plot all recs at once by mean of AC and spikewidth
figure; clf; hold on;
plot(ALL_ac_mean, ALL_spikewidth, 'ko')
plot(ALL_ac_mean_PYR, ALL_spikewidth_PYR, 'bo'); hold on;
plot(ALL_ac_mean_IN, ALL_spikewidth_IN, 'ro'); hold on;
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
title(['All recordings - mean AC vs spikewidth'])
filename = [savedatadir 'Spikewidth_MeanFR_MeanAC\allrecs_spikewidthvsmeanac'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%plot all recs at once by spikewidth and mean FR
figure; clf; hold on;
plot(ALL_spikewidth, ALL_FR_mean,'ko')
plot(ALL_spikewidth_PYR, ALL_FR_mean_PYR, 'bo'); hold on;
plot(ALL_spikewidth_IN, ALL_FR_mean_IN, 'ro'); hold on;
xlabel('Spikewidth - peak2trough diff')
ylabel('Mean FR (Hz)')
title(['All recordings - spikewidth vs meanFR'])
filename = [savedatadir 'Spikewidth_MeanFR_MeanAC\allrecs_spikewidthvsmeanFR'];
saveas(gcf,filename,'fig'); saveas(gcf,filename,'png');
close all

%% save the structure
%celltype data struct
for i = 1:size(dayindex,1)
    if makenewfiles.all || ~exist(filename, 'file')
        cellTypeProps = [];
        cellTypeProps{dayindex(i,1)}{dayindex(i,2)} = allrecProps(i);
        filename = [spikedatadir 'cellTypeProps_' iden num2str(dayindex(i,1)) '_' num2str(dayindex(i,2))];
        save(filename,'cellTypeProps');
    end
end
end