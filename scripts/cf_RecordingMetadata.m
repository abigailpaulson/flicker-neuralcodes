%RecordingMetadata
%   create a giant table of information about each recording day that
%   includes any outlier information like channel issues, bad channels,
%   etc. and also compiles information about the animal's experimental
%   group, flicker day, etc. 
% ALP 6/7/2022
%modified from script_createDataProperties
%ALP 12/01/2022

clear
%% get index information first
[params, dirs, indices, groups] = projectInfo('chronicflicker_annulartrack', ...
    'datethresh', 200800);

dayind = indices.dayindex;
allind = indices.allindex; 
fday = indices.flickerDay; 
recType = indices.recType; 

%% create variables for table
%%% index information
AnimalID = dayind(:,1); 
Date = dayind(:,2); 
Files = arrayfun(@(x) allind(allind(:,2) == dayind(x,2),3), 1:size(dayind,1), 'UniformOutput', false);
Files = Files';

%%% group information
isGroupGamma = ismember(AnimalID, groups.animals{1});
Groups = cell(length(AnimalID), 1); 
Groups(isGroupGamma) = groups.names(1); 
Groups(~isGroupGamma) = groups.names(2); 

%%% flicker day information
FlickerDay = indices.flickerDay; 

%%% recording type information
RecordingType = arrayfun(@(x) recType(allind(:,2) == dayind(x,2)), 1:size(dayind,1), 'UniformOutput', false);
RecordingType = RecordingType'; 

%%% hemisphere 
hemiNum = indices.hemisphere; 
Hemisphere = cell(length(AnimalID),1);
Hemisphere(hemiNum == 1) = {'R'};
Hemisphere(hemiNum == 2) = {'L'};

%%% brain regions
Regions = repmat({'CA1', 'CA3'}, [length(AnimalID),1]);

%%% channel information
Channels = repmat({0:63, 0:63}, [length(AnimalID), 1]);
Channels{Date == 201001,1} = 32:63;

%%% bad channel information
BadChannels = repmat({24, 7}, [length(AnimalID), 1]);

%%% any major notes? 
Notes = cell(length(AnimalID),1); 
Notes(Date == 201001) = {'No port C'};

metadata = table(AnimalID, Date, Files, RecordingType, Groups, FlickerDay, Hemisphere, Regions, Channels, BadChannels, Notes);
metadata.Properties.Description = 'Metadata for 2020-2021 chronic flicker experiment. Abigail Paulson, 2022'; 
metadata.Properties.VariableUnits = {'', 'YYMMDD', '', '0 VR only, 1 VR pre, 2 flicker, 3 VR post', '', 'days', '', '', '0 based', '0 based', ''}; 

%% save! 
save([dirs.savedatadir, 'experimentinfo\experimentmetadata.mat'], 'metadata');


