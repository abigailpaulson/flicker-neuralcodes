function [allrecProps, WF2exclude] = getspikewidthforcelltypeclassification_K2(dayindex, clusterinfo,...
    spikedatadir, savedfiguresdir, iden, makenewfiles)
%SP 3.14.19
%updated ALP for Kilosort2 structures

disp(['Calculating spikewidths for sess ' num2str(dayindex(1)) num2str(dayindex(2))])
filename = [spikedatadir 'cellTypeDataSpikewidth_' iden num2str(dayindex(1)), ...
    '_' num2str(dayindex(2)) '.mat'];
recinfo.iden = iden; 
recinfo.index = dayindex;
samprate = clusterinfo(1).samprate;

if ~exist(filename) || makenewfiles
    %% calculate peak to trough values using differential to find trough
    for clusIdx = 1:length(clusterinfo)
        clusterID = clusterinfo(clusIdx).ID; 
        
        %get waveform data
        WF = clusterinfo(clusIdx).WF.mn;
        if ~isnan(WF(1))
            spikewidthInfo = calcSpikewidth_K2(WF, recinfo, clusterID, clusIdx, samprate, [savedfiguresdir, 'ExampleWaveforms\', num2str(dayindex(2)), '\']);
        else
            spikewidthInfo.peak2troughDiff = nan;
            spikewidthInfo.peakIdxDiff = nan;
            spikewidthInfo.troughIdxDiff = nan;
            spikewidthInfo.fullWFDiff = nan;
            spikewidthInfo.sw2 = nan; 
        end
        allrecProps.peak2troughDiff(clusIdx) = spikewidthInfo.peak2troughDiff;
        allrecProps.peakIdxDiff(clusIdx) = spikewidthInfo.peakIdxDiff;
        allrecProps.troughIdxDiff(clusIdx) = spikewidthInfo.troughIdxDiff;
        allrecProps.fullWFDiff(clusIdx,:) = spikewidthInfo.fullWFDiff;
        allrecProps.SW2(clusIdx) = spikewidthInfo.sw2;
    end
    cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)} = allrecProps;
    save(filename, 'cellTypeDataSpikewidth');
else
    load(filename, 'cellTypeDataSpikewidth')
    allrecProps.peak2troughDiff = cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)}.peak2troughDiff;
    allrecProps.peakIdxDiff = cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)}.peakIdxDiff;
    allrecProps.troughIdxDiff = cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)}.troughIdxDiff;
    allrecProps.fullWFDiff = cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)}.fullWFDiff;
    allrecProps.SW2 = cellTypeDataSpikewidth{dayindex(1)}{dayindex(2)}.SW2;
end

%remove bad waveforms/clusters
filename = [spikedatadir 'badWFs_' iden num2str(dayindex(1)) '_' num2str(dayindex(2)) '.mat'];
if ~exist(filename, 'file')
    prompt = {['Enter MatIdx of bad units for:' num2str(dayindex(1)) ' ' num2str(dayindex(2))]};
    prompttitle = 'Exclude bad WFs for exclusion';
    definput = {''}; dims = [1 35];
    answer = inputdlg(prompt,prompttitle,dims,definput);
    badWFs{dayindex(1)}{dayindex(2)} = str2num(answer{1});
    save(filename, 'badWFs');
else
    load(filename)
end
WF2exclude = badWFs{dayindex(1)}{dayindex(2)};
allrecProps.peak2troughDiff(WF2exclude) = nan;
end