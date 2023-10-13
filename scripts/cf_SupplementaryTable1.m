%cf_SupplementaryTable1
%
%ALP 4/11/2023
clear; 

%% load general info
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');
dayindex = unique(allindex(:,1:2), 'rows');

%% initialize table
STable1 = table; 

%% directories
anprocessedatadir = '';
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
trialdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\'; 
ripplechandir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\CA1\';
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\celltypes\';


%% general animal and day information
recInfo = metadata(:,[1,2,5,6,7]);
recInfo.day = recInfo.Date;

%% load cell type information
for d = 1:size(dayindex,1)
    index = allindex(dayindex(d,2) == allindex(:,2),:);
    anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\'];
    [~, unitInfo, ~] = cf_getspikes(anprocesseddatadir, dirs, params, index, params.brainReg);
    
    nUnits_CA1(d,1) = sum(strcmp({unitInfo.brainReg}, 'CA1'));
    nUnits_CA3(d,1) = sum(strcmp({unitInfo.brainReg}, 'CA3'));
    nUnits(d,1) = length(unitInfo);
    nPYR(d,1) = sum(strcmp({unitInfo.cellType}, 'PYR'));
    nIN(d,1) = sum(strcmp({unitInfo.cellType}, 'IN'));
    day(d,1) = dayindex(d,2);
    
    unitInfo = [];
end

cellType = table(nUnits, nUnits_CA1, nUnits_CA3, nPYR, nIN, day);

%% load place cell information
for d = 1:size(dayindex,1)
    PCfilename = ['spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'];
    load([ratemapdir, PCfilename])
    
    nPC(d,1) = sum(ratemaps.includeCell);
    day(d,1) = dayindex(d,2);
    ratemaps = [];
end
placeCells = table(nPC, day);

%% load behavior information
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\behavior\behavioranalysis_prepost_recordings_table.mat']);
allTData = struct2table(allTData); 
inclTrials = allTData.engaged == 1 & allTData.fullTrial == 1 & allTData.rewarded == 1;
behaviorData = allTData(inclTrials,:);

for d = 1:size(dayindex,1)
    isDay = behaviorData.day == dayindex(d,2);
    nTrials(d,1) = sum(isDay);
    day(d,1) = dayindex(d,2);
end

Trials = table(nTrials,day);

%% load ripple information
savedatadir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\ripple_properties\'];
load([savedatadir, 'rippleproperties.mat'])

Ripples = AllDayData(:,1:2);

%% make supp table
STable1 = outerjoin(recInfo, cellType, 'MergeKeys', true);
STable1 = outerjoin(STable1, placeCells, 'MergeKeys', true);
STable1 = outerjoin(STable1, Trials, 'MergeKeys', true);
STable1 = outerjoin(STable1, Ripples, 'MergeKeys', true);
STable1 = STable1(STable1.FlickerDay > 5,:);

STableAnimal = groupsummary(STable1(:,[1,7:end]), "AnimalID", 'sum'); 
savedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\supplementTable1.csv'; 
writetable(STableAnimal, savedir)


