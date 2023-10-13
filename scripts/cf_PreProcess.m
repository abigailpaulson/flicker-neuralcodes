%cf_PreProcess
%   Run this first! 
%ALP 12/1/2022
close all; clear 

%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');
dirs.clusterappend = 'sorted\';
%% set run flags
run.posInfo = 0;
run.trialStructure = 0;
run.cellTypeClass = 1;
run.placecellclass = 0; 
run.trialStructure_ephys = 0; 
run.zscoregamma = 0; 

%% add other params
params.speedThreshold = 1; 
params.place_edges = 0:2:360;

%% make positionInfo structure
if run.posInfo
    [~,~,~, allindexfull] = projectInfo2('chronicflicker_annulartrack');
    cf_getpositioninfo(allindexfull, dirs.virmendir, dirs.processeddatadir, params)
end

%% make per trial behavior structure
if run.trialStructure
    cf_makebehaviortrialstructure(dirs, params, allindex, metadata, '180')
end

%% cell type classification
if run.cellTypeClass
    %[~,~,~, allindexfull] = projectInfo2('chronicflicker_annulartrack'); %all files, including flicker, for classifying stuff
    dayindex = unique(allindex(:,1:2), 'rows');
    isPost = metadata.FlickerDay > 5;
    inclDay = dayindex(isPost,:);
    isPost = ismember(allindex(:,2), inclDay);
    allindexpost = allindex(isPost,:);
    cf_getcelltypes(dirs, params, allindexpost, metadata)
end

%% place cell classification
[~,~,~, allindexfull] = projectInfo2('chronicflicker_annulartrack'); %all files, including flicker, for classifying stuff
if run.placecellclass
    params.placecells_postype = 'full';
    params.condition = 'all';
    params.nShuff = 500;
    cf_getplacecells(dirs, params, allindexfull, metadata)
    
    params.placecells_postype = 'd2r';
    params.condition = 'all';
    params.nShuff = 0;
    cf_getplacecells(dirs, params, allindexfull, metadata)
    
    params.placecells_postype = 'full';
    params.condition = 'preflicker';
    params.nShuff = 0;
    cf_getplacecells(dirs, params, allindex, metadata)
    
    params.placecells_postype = 'd2r';
    params.condition = 'preflicker';
    params.nShuff = 0;
    cf_getplacecells(dirs, params, allindex, metadata)
end

%% make ephys trial structure
%add spiking information to the behavior trial structure. also include
%helpful info like cell types and place cells
if run.trialStructure_ephys
    cf_maketrialstructure_ephys(dirs, params, allindex, metadata, '360')
end

%% get z-score gamma power for each recording day
if run.zscoregamma
    cf_zscorelfp(dirs, 'CA3', allindex, 'lowgamma')
end