function [singleunits, clustermetricsout] = cf_loadclusterfiles(clusterdir, index, brainReg)
%cf_loadclusterfiles
%   modified from loadClusterFiles_2reg Load Clusters#.mat for all files specified in index. 
% inputs:
%       clusterdir - location of cluster files
%       index - nx3
% outputs:
%       singleunits - structure
%       newclustermetrics - clustermetrics that is the same size as the new
%           clusters .mat file
%       brainReg - list of brainRegions in a cell array
% ALP 5/7/2020
% ALP 12/14/21 edited to skip loading spikes from one brain region if the
% cluster files don't exist
%ALP 12/19/22 edit to index into the structure with the file number instead
%   of the file index

clustermetricsout = [];

for br = 1:length(brainReg)
    if exist([clusterdir, brainReg{br}, '\sorted\clustermetrics.mat'], 'file')
        load([clusterdir, brainReg{br}, '\sorted\clustermetrics.mat'], 'clustermetrics')
    else
        if length(brainReg) == 1
            singleunits = []; 
            clustermetricsout = []; 
        end
        continue
    end
    
    for i = 1:size(index,1)
        if br == 1 || ~exist('singleunits', 'var')
            singleunits(index(i,3)).data = [];
        end
        load([clusterdir, brainReg{br}, '\sorted\clusters', num2str(index(i,3)), '.mat'], 'clusters')
       
        if ~isfield(clusters{index(i,1)}{index(i,2)}{index(i,3)}.data, 'spikeTimes')
            for c = 1:length(clusters{index(i,1)}{index(i,2)}{index(i,3)}.data)
                clusters{index(i,1)}{index(i,2)}{index(i,3)}.data(c).spikeTimes =  clusters{index(i,1)}{index(i,2)}{index(i,3)}.data(c).spikeInds./clusters{index(i,1)}{index(i,2)}{index(i,3)}.samprate;
            end
        end
        
        if ~isfield(clusters{index(i,1)}{index(i,2)}{index(i,3)}.data, 'brainReg')
            [clusters{index(i,1)}{index(i,2)}{index(i,3)}.data.brainReg] = deal(brainReg{br});
        end
        
        singleunits(index(i,3)).data = [singleunits(index(i,3)).data clusters{index(i,1)}{index(i,2)}{index(i,3)}.data];
        singleunits(index(i,3)).index = index(i,:); 
    end
    
    %only metrics for included clusters
    rawclustermetricsind = [clustermetrics.ID];
    finalclustersind = [clusters{index(i,1)}{index(i,2)}{index(i,3)}.data.ID];
    bothclus = ismember(rawclustermetricsind, finalclustersind);
    newclustermetrics = clustermetrics(bothclus);
    
    if ~isfield(newclustermetrics, 'brainReg')
        [newclustermetrics.brainReg] = deal(brainReg{br}); 
    end
    if ~isfield(newclustermetrics, 'index')
        [newclustermetrics.index] = deal(index(1,1:2)); 
    end
        
    clustermetricsout = [clustermetricsout newclustermetrics];   
end


end

