%chronic flicker behavior
%
%
%ALP 12/5/2022

%% load metadata
clear; clc; close all;
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%% run flags
run.behaviorExample = 0; 
run.fullTrackBehavior = 1;
run.behaviorAnalysis = 0;

%% params

%% per trial behavior example with ephys
if run.behaviorExample
    cf_getbehaviorexample(dirs, params, allindex, metadata)
end

%% full track average licking and velocity
if run.fullTrackBehavior
    cf_getbehavior_fulltrack(dirs, params, allindex, metadata)
end

%% calculate behavior differences between the groups
if run.behaviorAnalysis
    cf_getbehavior_recordings(dirs, params, allindex, metadata)
end

