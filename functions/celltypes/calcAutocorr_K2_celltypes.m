function [autocorr, centerofmass] = calcAutocorr_K2_celltypes(clusterdir, clusterInfo)
%this function calculates the autocorrelogram of clusters from a recording
%SP 4.16.18
%modified to use spike indicies to calculate the autocorrelograms and not
%spiketimes NJ 08.14.19

%INPUT
%clusterdir: directory to rawclusters.mat files

%OUTPUT (for each cluster)
%spiketrain
%autocorrelogram
%center of mass
%cell ID

%fixed a bug related to indices, now uses spike indicies of all recs per
%unit 09.06.19
% updated for use in cell type classifier ALP 1/16/2020

%load sorting props
load(fullfile(clusterdir, 'kilosort', 'sortingprops.mat'))
load(fullfile(clusterdir, 'rawclusters_allrec.mat')) %edit from rawclusters(recNum) bc indices all restart from 1: NJ 09.06.19

goodunits = ismember([rawclusters_allrec.ID], [clusterInfo.ID]);
rawclusters_allrec = rawclusters_allrec(goodunits); 

for unit = 1:size(rawclusters_allrec,2)
    totalspikes{unit} = rawclusters_allrec(unit).spikeInds'; %change vertical to horizontal structure
end

%total length of recording day, same for all units
totalsamples = sum(props.recLength);

for unit = 1:size(totalspikes,2)
    %make the spike train - using indices instead of times
    stepsize = 5 * props.sampRate / 1000; %number of samples for 5ms binsize
    spiketrainedges = 0:stepsize:totalsamples; %5ms bins
    spiketrain{unit} = histc(totalspikes{unit}, spiketrainedges);
    
    %get the autocorr
    lag_num = 50 * props.sampRate / 1000; %number of samples for 50ms
    lag = lag_num/stepsize; %in bins
    autocorr{unit} = xcorr(spiketrain{unit},lag);
    autocorr{unit}(lag+1) = 0;
    
    %eliminate ones with not enough spikes
    if max(autocorr{unit}) < 10
        autocorr{unit} = nan(1,length(autocorr{unit}));
    end
end

%get the first moment of the autocorr
sampN = stepsize:stepsize:lag_num;
for unit = 1:length(autocorr)
    centerofmass{unit} = (sum(autocorr{unit}(lag+2:end).*sampN)/sum(autocorr{unit}(lag+2:end)))/stepsize;
end

%save structure
% save(fullfile(clusterdir, 'autocorr.mat'), 'fr') %choosing to save
% outside of this function ALP 1/16/2020
end

