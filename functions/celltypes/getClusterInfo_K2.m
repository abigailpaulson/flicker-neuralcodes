function [clusterInfo] = getClusterInfo_K2(datadir, dayindex, recs)
% this function returns gets cluster info from Mountain Lab files
% inputs: 
%        datadir of sorted unit info
%        dayindex - [animal day]
%        recs - list of recording file numbers
% outputs:
%        clusterInfo: structure of cluster metrics, waveforms, main chans,
%        recInfo with spikes, vector with included file numbers, duration
%        vector with recording durations of the specified files in seconds.
%
% SP 12.19.17
% Updated ALP 1/14/2020 for Kilosort 2 datastructures

%% load sorting props for recording duration
if exist([datadir, 'kilosort\sortingprops.mat'], 'file')
    load([datadir, 'kilosort\sortingprops.mat'])
else
    %I know at least one day is missing sorting props bc Kilosort failed
    disp(['cluster files missing for ', num2str(dayindex)])
    
    for i = 1:length(recs)
        clusterInfo.recInfo(i).singleunits = [];
        clusterInfo.recInfo(i).file = recs(i); 
        clusterInfo.recInfo(i).duration = [];
        clusterInfo.metrics = [];
    end
    return
end

%% get spikes
for i = 1:length(recs)
    load([datadir, 'clusters', num2str(recs(i)), '.mat'], 'clusters')    
    clusterInfo.recInfo(i).singleunits = clusters{dayindex(1,1)}{dayindex(1,2)}{recs(i)}.data;
    clusterInfo.recInfo(i).file = recs(i); 
    clusterInfo.recInfo(i).duration = props.recLength(props.fileNums == recs(i))/props.sampRate; %duration in seconds
end

%% get cluster metrics
load([datadir 'clustermetrics.mat'], 'clustermetrics');

%only metrics for included clusters
rawclustermetricsind = [clustermetrics.ID];
finalclustersind = [clusterInfo.recInfo(end).singleunits.ID]; 
bothclus = ismember(rawclustermetricsind, finalclustersind); 
newclustermetrics = clustermetrics(bothclus); 
clusterInfo.metrics = newclustermetrics;


end