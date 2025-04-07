%cf_RunAnalyses
%
%ALP 12/13/2022

clear; close all
%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%% set run flags
run.behaviorAnalysis = 1;
run.spatialMapAnalysis = 0; 
run.currPositionDecoding = 0;
run.thetaSeqDecoding = 0; 
run.rippleDecoding = 0;
run.flickerSequences = 0; 
run.rippleProperties = 0;
run.rippleActivation = 0;
run.rippleRateRZ = 0;

%% run stuff

if run.behaviorAnalysis 
    cf_getbehavior_recordings_table(dirs, params, allindex, metadata)
end

if run.spatialMapAnalysis
    cf_spatialmaps_fieldproperties_table(dirs, params, allindex, metadata, 'd2r')
end

if run.currPositionDecoding
    cf_decoding_currentposition(dirs, params, allindex, metadata, 'full', 'PC')
    cf_decoding_currentposition(dirs, params, allindex, metadata, 'd2r', 'PC')
end

if run.thetaSeqDecoding
    cf_decoding_thetaseq_360dec_180trialsplit_seqlength(dirs, params, allindex, metadata, 'full', 'PC') % correct control speed
end

if run.rippleDecoding
%     allindex = allindex(allindex(:,2) == 200831,:); %example day
    %cf_decoding_ripples(dirs, params, allindex, metadata, 'd2r', 'PC')
    %cf_decoding_ripples_table(dirs, params, allindex, metadata, 'd2r', 'PC')
    cf_decoding_ripples_table(dirs, params, allindex, metadata, 'full', 'PC')
end

if run.rippleProperties
    cf_ripple_properties(dirs, params, allindex, metadata)
end

%activation and coactivation of place cells during ripples 
if run.rippleActivation 
    cf_placecell_rippleactivation_RZripples(dirs, params, allindex, metadata)
end

if run.rippleRateRZ
    cf_rippleproperties_rate_RZ(dirs, params, allindex, metadata)
end