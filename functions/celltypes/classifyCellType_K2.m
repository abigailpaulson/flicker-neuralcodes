function cellTypeInfo = classifyCellType_K2(allindex, processeddatadir, savedfiguresdir, ...
    clusterdirappend, params, iden)
%SP 4.27.18, classify as PYR/IN
%Edited ALP 6.28.18
%Updated ALP 1/14/2020 for Kilosort2 pipeline datastructures
% INPUTS:
%   params.raweegsamprate - in Hz
%   params.tBeforeSpike - in ms, usually 1ms
%   params.tAfterSpike - in ms, usually 3ms
%
% OUTPUTS:

%% setup
allindex2 = allindex(:, 1:2);
dayindex = unique(allindex2,'rows');
spikedatadir = [savedfiguresdir 'celltypeProps\'];

if ~exist(spikedatadir, 'dir')
    mkdir(spikedatadir)
end

% make filter for eeg waveforms
samprate = params.raweegsamprate;
hpFilt = designfilt('highpassiir', 'StopbandFrequency', 100, ...
     'PassbandFrequency', 500, 'StopbandAttenuation', 60, ...
     'PassbandRipple', 1, 'SampleRate', samprate, 'DesignMethod', 'butter');
 
 %% load up clusters and cluster metrics
 for i = 1:size(dayindex,1)
     if ~exist([spikedatadir 'cellTypeProps_' iden(i) num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'], 'file')
         % will want to write in a "if ~exist(cellTypeProps)" here, eventually
         clusterdatadir = fullfile(processeddatadir, iden(i), num2str(dayindex(i,1)),...
             '_', num2str(dayindex(i,2)), clusterdirappend);
         recs = allindex(allindex(:,2) == dayindex(i,2),3);
         recinfo.index = dayindex(i,:);
         recinfo.iden = iden;
         clusterInfo(i).metrics = getClusterInfo_K2(clusterdatadir, dayindex, recs);
         
         allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR = clusterInfo(i).metrics.meanFR;
         allrecProps{dayindex(i,1)}{dayindex(i,2)}.peakFR = clusterInfo(i).metrics.peakFR;
         allrecProps{dayindex(i,1)}{dayindex(i,2)}.stabletimes = clusterInfo(i).metrics.stabletimes;
         allrecProps{dayindex(i,1)}{dayindex(i,2)}.waveforms = clustermetrics2mat(clusterInfo(i).metrics, 'WF', 'mn');
         
         for clusIdx = 1:length(clusterInfo(i).metrics)
             spikewidthInfo = calcSpikewidth_K2(allrecProps{dayindex(i,1)}{dayindex(i,2)}.waveforms(clusIdx,:),...
                 recinfo, clusterInfo(i).metrics(clusIdx).ID, clusterInfo(i).metrics(1).samprate, savedfiguresdir);
             allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff(clusIdx) = spikewidthInfo.peak2troughDiff;
             allrecProps{dayindex(i,1)}{dayindex(i,2)}.peakIdxDiff(clusIdx) = spikewidthInfo.peakIdxDiff;
             allrecProps{dayindex(i,1)}{dayindex(i,2)}.troughIdxDiff(clusIdx) = spikewidthInfo.troughIdxDiff;
             allrecProps{dayindex(i,1)}{dayindex(i,2)}.fullWFDiff(clusIdx,:) = spikewidthInfo.fullWFDiff;
         end
     else
         load([spikedatadir 'cellTypeProps_' iden(i) num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '.mat'], 'cellTypeProps')
         allrecProps{dayindex(i,1)}{dayindex(i,2)} = cellTypeProps{dayindex(i,1)}{dayindex(i,2)}; 
         clear cellTypeProps
     end
 end
 
 %% get FR distribution
 for i = 1:size(dayindex,1)
     figure; clf; hold on;
     edges = 0:0.5:30;
     histogram(allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR, edges)
     ylabel('Number of Clusters')
     xlabel('Mean FR')
     title(['FR distribution- ', iden(i), num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
     
     %save figures
     datadir = [savedfiguresdir 'FRDistributions/'];
     savefigSP(dayindex(i,:), datadir, 'distFR_mean', iden(i))
 end
 close all
 
%% remove bad waveforms/clusters

%remove bad waveforms/clusters
%bad WFs (found by manually going through all the waveforms)
 
for i = 1:size(dayindex,1)
     WF2exclude = badWFs{dayindex(i,1)}{dayindex(i,2)}{probe};
    allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff(WF2exclude) = nan;
end

%spikewidth distribution
for i = 1:size(dayindex,1)
    figure; clf; hold on;
    edges = 0:0.01:1.5;
    histogram(allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff, edges)
    ylabel('Number of Clusters')
    xlabel('Spikewidth')
    title(['Spikewidth Dist - ', iden(i), num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
    
    %save figures
    datadir = [savedfiguresdir 'SpikewidthDistributions/'];
    savefigSP(dayindex(i,:), datadir, 'peak2troughDiff_dist', iden(i));
end

%spikewidth distribution for all recordings
ALL_spikewidth = [];
for i = 1:size(dayindex,1)
     spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
    ALL_spikewidth = [ALL_spikewidth spikewidth];
end
figure; clf; hold on;
edges = 0:0.025:1.5;
histogram(ALL_spikewidth, edges); hold on;
% plot([0.5-0.0125 0.5-0.0125],[0 40], 'k--')
ylabel('Number of Clusters')
xlabel('Spikewidth')
title(['Spikewidth Dist - ALL RECS']);

%save figures
datadir = [savedfiguresdir 'SpikewidthDistributions/'];
filename = [datadir 'ALLRECS_peak2troughDiff_dist'];
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

%% get the autocorrelogram data
% can use is field to check and see if I need to calculate for each one 
filename = [spikedatadir 'WFstructure_', iden(i), num2str(dayindex(end,1)) '_' num2str(dayindex(end,2)) '.mat']; 
if ~exist(filename)
for i = 1:size(dayindex,1) 
    clusterdatadir = [processeddatadir iden(i) num2str(dayindex(i,1)) '_' num2str(dayindex(i,2)) '/Port ', probeletter,' Sorted/'];
    tempindex = allindex(allindex(:,2) == dayindex(i,2),:);
    [autocorr{dayindex(i,1)}{dayindex(i,2)}, centerofmass] = calcAutocorrelogram(tempindex, clusterdatadir);
    autocorr_mean{dayindex(i,1)}{dayindex(i,2)} = centerofmass;
    
    %store in allrecProps data structure
    allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean = centerofmass;
    allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr = autocorr{dayindex(i,1)}{dayindex(i,2)};
end
end

% plot the individual autocorrelograms and their respective WFs
for i = 1:size(dayindex,1)
    for unit = 1:length(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr)
        acorr = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr(unit));
        acorr_mean = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean(unit));
        wf = allrecProps{dayindex(i,1)}{dayindex(i,2)}.fullWFDiff(unit,:);
        
        if ~exist(filename)
            plotAutocorrWF(dayindex(i,:), acorr, acorr_mean, wf, unit, savedfiguresdir);
        end
    end
    close all
end

%get the AC distribution
ALL_autocorr = [];
for i = 1:size(dayindex,1)
    figure; clf; hold on;
    edges = 0:0.2:10;
    histogram(cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean), edges)
    ylabel('Number of Clusters')
    xlabel('Autocorrelogram Mean')
    title(['Autocorrelogram distribution - ', iden(i), num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);
    
    %save figures
    datadir = [savedfiguresdir 'AutocorrelogramDistributions/'];
    savefigSP(dayindex(i,:), datadir, 'distAC', iden(i))
    
    ALL_autocorr = [ALL_autocorr cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean)];
end

figure; clf; hold on;
edges = 0:0.2:10;
histogram(ALL_autocorr, edges)
ylabel('Number of Clusters')
xlabel('Autocorrelogram Mean')
title(['Autocorrelogram distribution - ', iden(i), num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))]);

datadir = [savedfiguresdir 'SpikewidthDistributions/'];
filename = [datadir 'ALLRECS_meanautocorr_dist'];
saveas(gcf,filename,'fig');

%% plot the 3 parameters
for i = 1:size(dayindex,1)
    ac_mean = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean);
    spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
    FR_mean = allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR;
    
    figure; clf; hold on;
    plot3(ac_mean, spikewidth, FR_mean, 'o')
    xlabel('Mean ac (ms)')
    ylabel('Spikewidth - peak2trough diff')
    zlabel('Mean FR (Hz)')
    title(['Mean AC vs. Peak2Trough Diff  vs. Mean FR - ', iden(i), num2str(dayindex(i,1)) ' ' num2str(dayindex(i,2))])
    
    %save figures
    datadir = [savedfiguresdir 'Spikewidth_MeanFR_MeanAC/'];
    savefigSP(dayindex(i,:), datadir, 'MeanFR_MeanAC_peak2troughDiff', iden(i));
end

%% plot all the recordings at once
ALL_ac_mean = [];
ALL_spikewidth = [];
ALL_FR_mean = [];
for i = 1:size(dayindex,1)
    ac_mean = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean);
    spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
    FR_mean = allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR;
    
    ALL_ac_mean = [ALL_ac_mean ac_mean];
    ALL_spikewidth = [ALL_spikewidth spikewidth];
    ALL_FR_mean = [ALL_FR_mean FR_mean];
end
figure; clf; hold on;
plot3(ALL_ac_mean, ALL_spikewidth, ALL_FR_mean, 'o')
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title(['ALL RECS - Mean AC vs. Peak2Trough Diff vs. FR'])

%save figures
datadir = [savedfiguresdir 'Spikewidth_MeanFR_MeanAC/'];
if ~exist(datadir); mkdir(datadir); end
filename = [datadir 'ALLRECS_peak2troughDiff'];
saveas(gcf,filename,'fig');

%% separate the groups manually and then plot the traces
%set thresholds
if strcmp(probeletter, 'B')
    thresholds.widththresh_IN = 0.5; % long AC - 0.6
    thresholds.acthresh_IN = 4; %8;
    thresholds.widththres_PYR = 0.5; %short and long AC
    thresholds.acthresh_PYR = 4; %9;
else
    thresholds.widththresh_IN = 0.6; % long AC - 0.6
    thresholds.acthresh_IN = 4; %8;
    thresholds.widththres_PYR = 0.6; %short and long AC
    thresholds.acthresh_PYR = 5; %9;
end

%individual waveforms
for i = 1:size(dayindex,1)
    %load up parameter info
    WFprops.ac_mean = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean);
    WFprops.spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
    WFprops.FR_mean = allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR;
    WFprops.wf = allrecProps{dayindex(i,1)}{dayindex(i,2)}.fullWFDiff;
    
    %get indices and plot separation
    cellTypeInfo{dayindex(i,1)}{dayindex(i,2)} = getCellTypeIdx(dayindex(i,:),WFprops,thresholds, iden(i), savedfiguresdir);
    allrecProps{dayindex(i,1)}{dayindex(i,2)}.cellTypeInfo = cellTypeInfo{dayindex(i,1)}{dayindex(i,2)};
    allrecProps{dayindex(i,1)}{dayindex(i,2)}.thresholds = thresholds;
end


%plot average waveforms for all recordings
ALL_PYRwfs = [];
ALL_INwfs = [];
for i = 1:size(dayindex,1)
    %get indices
    WFprops.wf = allrecProps{dayindex(i,1)}{dayindex(i,2)}.fullWFDiff;
    INIdx = cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.INidx;
    PYRIdx = cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.PYRidx;
    %get mean wfs
    INwf = WFprops.wf(INIdx,:);
    PYRwf = WFprops.wf(PYRIdx,:);
    ALL_PYRwfs = [ALL_PYRwfs; PYRwf];
    ALL_INwfs = [ALL_INwfs; INwf];  
end      
%get mean and SEM
meanINwf = nanmean(INwf);
meanPYRwf = nanmean(PYRwf);
stdINwf = nanstd(INwf);
stdPYRwf = nanstd(PYRwf);
stdINwf_line1 = meanINwf+stdINwf;
stdINwf_line2 = meanINwf-stdINwf;
stdPYRwf_line1 = meanPYRwf+stdPYRwf;
stdPYRwf_line2 = meanPYRwf-stdPYRwf;

figure; hold on;
time = 0:0.125:45;
plot(time,meanINwf, 'r-','LineWidth', 3)
plot(time,meanPYRwf, 'b-','LineWidth', 3)
plot(time,stdINwf_line1,'r--')
plot(time,stdINwf_line2,'r--')
plot(time,stdPYRwf_line1,'b--')
plot(time,stdPYRwf_line2,'b--')
title(['Average waveforms (width threshold - ' num2str(thresholds.widththres_PYR) ')'])
legend('Interneurons','Pyramidal')
ylim([-1 1])
%save the figure
datadir = [savedfiguresdir 'Classification/'];
filename = [datadir 'ALLRECS_avgwaveforms'];
saveas(gcf,filename,'fig');
saveas(gcf,filename,'png');

%% Plot all with cell classifications
ALL_ac_mean = [];
ALL_spikewidth = [];
ALL_FR_mean = [];

figure; clf; hold on;
for i = 1:size(dayindex,1)
    ac_mean = cell2mat(allrecProps{dayindex(i,1)}{dayindex(i,2)}.autocorr_mean);
    spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
    FR_mean = allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR;
    
    INIdx = cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.INidx;
    PYRIdx = cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.PYRidx;
    EXLIdx = find(isnan(cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.celltype));
    
    plot3(ac_mean(INIdx), spikewidth(INIdx), FR_mean(INIdx), 'bo')
    plot3(ac_mean(PYRIdx), spikewidth(PYRIdx), FR_mean(PYRIdx), 'ro')  
    plot3(ac_mean(EXLIdx), spikewidth(EXLIdx), FR_mean(EXLIdx), 'ko', 'MarkerEdgeColor', [0.85 0.85 0.85])  
end

% edges = 0:0.2:10;
% a = histc(ALL_autocorr, edges);
% bar(edges, (-1).*a)
% edges = 0:0.025:1.5;
% b = histc(ALL_spikewidth, edges);
% bar((-1).*edges, b)
xlabel('Mean ac (ms)')
ylabel('Spikewidth - peak2trough diff')
zlabel('Mean FR (Hz)')
title('ALL RECS - Mean AC vs. Peak2Trough Diff vs. FR')

%save figures
datadir = [savedfiguresdir 'Spikewidth_MeanFR_MeanAC/'];
if ~exist(datadir); mkdir(datadir); end
filename = [datadir 'ALLRECS_peak2troughDiff_2'];
saveas(gcf,filename,'fig');


%% Plot excluded waveforms 

% ALL_PYRwfs = [];
% ALL_INwfs = [];
% time = 0:0.125:45;
% datadir = [savedfiguresdir 'ExcludedCells/'];
% if ~exist(datadir); mkdir(datadir); end
% for i = 1:size(dayindex,1)
%     %get indices
%     WFprops.wf = allrecProps{dayindex(i,1)}{dayindex(i,2)}.fullWFDiff;
%     EXLIdx = find(isnan(cellTypeInfo{dayindex(i,1)}{dayindex(i,2)}.celltype));
%     
%     ac_mean = cell2mat(autocorr_mean{dayindex(i,1)}{dayindex(i,2)});
%     spikewidth = allrecProps{dayindex(i,1)}{dayindex(i,2)}.peak2troughDiff;
%     FR_mean = allrecProps{dayindex(i,1)}{dayindex(i,2)}.meanFR;
%     
%     for cell = 1:length(EXLIdx)
%         figure; hold on; 
%         subplot(1,2,1); hold on; 
%         plot3(ac_mean(EXLIdx(cell)), spikewidth(EXLIdx(cell)), FR_mean(EXLIdx(cell)), 'ko')
%         xlabel('Mean ac (ms)')
%         ylabel('Spikewidth - peak2trough diff')
%         zlabel('Mean FR (Hz)')
%         
%         subplot(1,2,2); hold on; 
%         plot(time, WFprops.wf(EXLIdx(cell),:), 'k-', 'LineWidth', 3)
%         title(['Excluded Cell ', num2str(EXLIdx(cell))])
%         
%         filename = [datadir 'EXLCELLS_', num2str(dayindex(i,1)),'_', num2str(dayindex(i,2)),'_Unit', num2str(EXLIdx(cell))];
%         saveas(gcf,filename,'png');
%         saveas(gcf,filename,'fig');
%     end
% 
% end      

%% save the structure
spikedatadir = [savedfiguresdir 'CellTypeProps/'];
if ~exist(spikedatadir)
    mkdir(spikedatadir); 
end

for i = 1:size(dayindex,1)
    cellTypeProps = {}; 
    filename = [spikedatadir 'cellTypeProps_' iden(i) num2str(dayindex(i,1)) '_' num2str(dayindex(i,2))];
    cellTypeProps{dayindex(i,1)}{dayindex(i,2)} = allrecProps{dayindex(i,1)}{dayindex(i,2)};
    save(filename,'cellTypeProps');
end

