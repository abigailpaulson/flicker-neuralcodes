function [params, dirs, indices, groups, allindex, dayindex] = projectInfo(project, varargin)
%projectInfo Define project parameters and directories for the input
%project.
%   inputs:     project     string with project name ex 'chronicflicker'
%
%   outputs:    params      structure with 
%
%   example 
%   call:       [params, dirs] = projectInfo('chronicflicker'); 
%
% ALP 4/15/2020
%
% Change Log
%   Date   Who     What Changed 
%   12/1/22 Abby    copied to manuscript folder

for v = 1:2:length(varargin)-1
    switch varargin{v}
        case 'datethresh'
            datethresh = varargin{v+1};
        case 'index'
            indextoget = varargin{v+1};
        case 'day'
            daytoget = varargin{v+1};
        otherwise
            error('option input not understood')
    end
end


if strcmp(project, 'chronicflicker_annulartrack')
    %%% define directories
    dirs.rawdatadir = '\\ad.gatech.edu\bme\labs\singer\RawData\Flicker_Chronic_VR\';
    dirs.processeddatadir = '\\ad.gatech.edu\bme\labs\singer\ProcessedData\Flicker_7Day_VR\'; 
    dirs.virmendir = '\\ad.gatech.edu\bme\labs\singer\Virmen Logs\ChronicFlicker_AnnularTrack\';
    dirs.savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\'; 
    dirs.spreadsheetdir = '\\ad.gatech.edu\bme\labs\singer\Abby\experimentspreadsheets\chronicflicker_annulartrack_ephys.xls';
%     dirs.celltypedir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\celltype\';
%     dirs.placefielddir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\';
%     dirs.projectdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\'; 
%     dirs.xcorrdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\monosynconnex\'; 
%     
    %%% define some common parameters
    params.brainReg = {'CA1', 'CA3'};
    params.nShanks = 2;
    params.iden = 'A'; 
    params.nChannels = 64; 
    params.badChannels{1} = 24; 
    params.badChannels{2} = 7; 

    
    %%% get all index
    %!!! data threshold, include all dates after this one!!!%
    if ~exist('datethresh', 'var') %auto fill this if another isn't put in
        datethresh = 200800;
    end
    
    [allindex, data] = getallindexALP(dirs.processeddatadir, dirs.spreadsheetdir, 0);
    allindex = allindex(allindex(:,2) > datethresh,:);
    
    if exist('daytoget', 'var')
        allindex = allindex(allindex(:,2) == daytoget,:); 
    end
    %%% for old data
   % allindex = allindex(allindex(:,2) < 200800,:);
    %%% end for old data
    %allindex = allindex(allindex(:,2) ~= 200911,:); % gap
    %allindex = allindex(allindex(:,2) ~= 200917,:); % gap
    %allindex = allindex(allindex(:,2) ~= 201001,:); % no port C
    %allindex = allindex(allindex(:,2) ~= 210127,:); %port c and d switched
    %allindex = allindex(allindex(:,2) ~= 201103,:); % gap
    %allindex = allindex(allindex(:,2) < 210101,:); %don't include 2021 data yet
    allindex = allindex(allindex(:,2) ~= 210212,:); %exclude bc only one cell this day, random post
    allindex = allindex(allindex(:,2) ~= 210514,:); %exclude bc no CA1 cells this day, gamma pre
    allindex = allindex(allindex(:,2) ~= 201029,:); %trigger issue first rec, excluding for now, gamma post
%     allindex = allindex(allindex(:,2) ~= 210206,:); %exclude bc high amplitude noise 
    dayindex = unique(allindex(:,1:2), 'rows');
    [~, flickerDayI, ~] = unique(allindex(:,1:2), 'rows');
    
    %get hemisphere
    datIn = cell2mat(data(:,2:3)); 
    hemisphere = []; 
    count = 1; 
    [uniqueDays, iUniqueDay, ~] = unique(datIn, 'rows');
    for d = 1:size(uniqueDays,1)
        if ismember(uniqueDays(d,2), dayindex(:,2))
            if strcmp(data{iUniqueDay(d), 13}, 'R')
                hemisphere(count) = 1; 
            elseif strcmp(data{iUniqueDay(d),13} , 'L')
                hemisphere(count) = 2;
            end
            count = count+1;
        end
    end
    
    [behaviorindex] = getallindexALP(dirs.processeddatadir, dirs.spreadsheetdir,...
        'stimulation', 'v'); 
    behaviorindex = behaviorindex(behaviorindex(:,2) > datethresh,:); 
    
    if exist('indextoget', 'var')
        allindex = indextoget; 
        dayindex = indextoget(:,1:2); 
        behaviorindex = behaviorindex(behaviorindex(:,2) == allindex(1,2),1:3); 
        flickerDayI = 1; 
    end
    
    indices.allindex = allindex(:,1:3); 
    indices.dayindex = dayindex;
    indices.flickerDay = allindex(flickerDayI, 4); %same length as dayindex
    indices.behaviorindex = behaviorindex(:,1:3);
    indices.recType = allindex(:,5); %0 is behavior only day, 1 is pre-flicker files, 2 is flicker, 3 is post-flicker files
    indices.hemisphere = hemisphere;
   
    %%% define groups
    groups.names = {'gamma', 'random'}; %40, random
    %groups.animals = {[22, 25], [21, 24, 26]}; %old data. check these
    %   groups correspond to gamma and random if using these animal number
    groups.animals = {[28, 30, 32, 35, 39, 41, 45, 47], [29, 31, 36, 37, 40, 42, 46, 48]}; 
    groups.colors = {[166 206 227; 31 120 180] , ...
        [178 223 128; 51 160 44]}; %blue and green
    groups.colors{1} = groups.colors{1}./256; 
    groups.colors{2} = groups.colors{2}./256; 
    groups.timepoint = {'pre', 'postchronic'};
    
    
    %%% make session info table...trying this out ALP 7/16/2020
    indices.SessionInfo = table(); 
    indices.SessionInfo.animal = allindex(:,1); 
    indices.SessionInfo.date = allindex(:,2); 
    indices.SessionInfo.file = allindex(:,3);
    indices.SessionInfo.flickerDay = allindex(:,4); 
    indices.SessionInfo.recType(allindex(:,5) <=1 ) = {'b1'};
    indices.SessionInfo.recType(allindex(:,5) == 2) = {'f'};
    indices.SessionInfo.recType(allindex(:,5) == 3) = {'b2'};
    isGroup1 = ismember(allindex(:,1), groups.animals{1});
    indices.SessionInfo.group(isGroup1) = groups.names(1); 
    indices.SessionInfo.group(~isGroup1) = groups.names(2); %if its not group 1, must be group2
end
end

