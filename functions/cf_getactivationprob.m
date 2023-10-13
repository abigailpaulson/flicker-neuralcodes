function out = cf_getactivationprob(unitIDs, spikes, ...
    rippleTimes, pairs, regions)
%getactivationprob
%   get activation and coactivation probability for spike times given
%   ripple times
%   inputs:
%       unitIDs - [nx1] or [1xn] vector of unitIDs
%       spikes - structure of length nUnits, where one field is spikeTimes
%           in seconds
%       rippleTimes - [nx2] where n=nRipples, column1 is startTime in
%           seconds and column2 is end time in seconds
%       pairs - [nx2] where n=nUnits, pair indices. input pair indices when
%           you wish to define the pairs used to calculate coactivation.
%           LEAVE EMPTY to calculate coactivation for all possible pairs
% 
%   outputs:
%       out.activity - [nx2] where n=nUnits, column1 is the number of
%           active and stable ripples, and column2 is the number of stable
%           ripples
%       out.coactivity - [nx2] where n=nPairs, column1 is the number of
%           ripples for which both cells in the pair were stable and active, 
%           and column2 is the number of ripples both cells in the pair were 
%           stable for
%       out.pairs - list of pairs unitID(pair(1,1)) and unitID(pair(1,2))
%           cell IDs of the first pair
%       out.unitIDs - list of unitIDs that were input and used to caluclate
%           pairs
%   To get probability of activation or coactivation, simply divide the
%   first column by the second. 
%      probActivity = out.activity(:,1)/out.activity(:,2) 
%
% ALP 9/27/21 chronicflicker-ephys
% ALP 3/19/23 based on getactivationprob, moved to chronicflicker-ephys-prospectivecoding
%   edited to remove stability critera

%% check inputs
if size(rippleTimes,2) > 2
    rippleTimes = rippleTimes';
    disp('check dimension of rippleTimes')
end

if length(spikes) ~= length(unitIDs)
    error('the size of the spike time input is not equivalent to the unitID input')
end

if ~isfield(spikes, 'spikeTimes')
    error('your spike structure does not contain a field with spike times in seconds')
end

%% create activity matrices 
activity = zeros(length(unitIDs), size(rippleTimes,1)); 

for iRip = 1:size(rippleTimes,1)
    for iCell = 1:length(unitIDs)
         if sum(isExcluded(spikes(iCell).spikeTimes, rippleTimes(iRip,:))) > 0
             activity(iCell, iRip) = 1; 
         end
    end
end

%% get activity prob
probActivity = zeros(length(unitIDs), 2);
nRipples = size(rippleTimes,1); 

for iCell = 1:length(unitIDs)
    isActive = logical(activity(iCell,:));
    nActive = sum(isActive); %cell is active during ripple and is stable during the ripple time
    
    probActivity(iCell,:) = [nActive nRipples]; 
    
    activeRegions{iCell} = regions{iCell}; 
    clear isActive nActive
end

%% get coactivity prob
if isempty(pairs)
    if length(unitIDs) > 1
        pairs = nchoosek(1:length(unitIDs),2); %work in cell indices
    else
        pairs = [];
    end
end

probCoactivity = zeros(size(pairs,1), 2); 
for iPair = 1:size(pairs,1)
    bothactivity = sum([activity(pairs(iPair,1),:); activity(pairs(iPair,2),:)],1);
    isActive = bothactivity == 2; 
    
    nActive = sum(isActive); 
    probCoactivity(iPair,:) = [nActive nRipples]; 
    
    coactiveRegions{iPair} = strcat(regions{pairs(iPair,1)}, regions{pairs(iPair,2)});
    
    clear bothactivity isActive nActive
end

%% organize outputs
out.activity  = probActivity; 
out.coactivity = probCoactivity; 
out.nActive = probActivity(:,1); 
out.nAcRipples = probActivity(:,2); 
out.nCoactive = probCoactivity(:,1); 
out.nCoRipples = probCoactivity(:,2); 
out.nAllRipples = size(rippleTimes,1); 
out.activeRegions = activeRegions;
out.coactiveRegions = coactiveRegions;
out.pairs = pairs;
out.unitIDs = unitIDs; 
out.distance = [];
end

