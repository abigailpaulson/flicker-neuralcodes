function [params, dirs, metadata, allindex] = projectInfo2(project, varargin)
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
% 12/3/22   Abby    updates for manuscript folders

for v = 1:2:length(varargin)-1
    switch varargin{v}
        case 'datethresh'
            datethresh = varargin{v+1};
        case 'index'
            indextoget = varargin{v+1};
        case 'day'
            daytoget = varargin{v+1};
        case 'condition'
            condition = varargin{v+1};
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
    dirs.celltypedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\celltypes\';
    dirs.placefielddir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\';
    dirs.projectdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker-ephys-prospectivecoding\results\'; 
    dirs.xcorrdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\monosynconnex\'; 
    
    %%% define some common parameters
    params.brainReg = {'CA1', 'CA3'};
    params.nShanks = 2;
    params.iden = 'A';
    params.nChannels = 64;
    params.badChannels{1} = 24;
    params.badChannels{2} = 7;
    
    load('\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\experimentinfo\experimentmetadata.mat')
    
    if exist('condition', 'var')
        if strcmp(condition, 'prepostVR')
            allindex = []; 
            dayindex = [metadata.AnimalID metadata.Date];
            for d = 1:height(metadata)
                tmpdayindex = dayindex(d,:); 
                tmpfiles = cell2mat(metadata{d,'Files'});
                tmprectype = cell2mat(metadata{d,'RecordingType'}); 
                tmpday = metadata{d, 'FlickerDay'};
                
                if tmpday < 5
                    inclFiles = tmprectype < 2;
                elseif tmpday > 5
                    inclFiles = tmprectype ~= 2;
                end
                
                newfiles = tmpfiles(inclFiles);
                newtype = tmprectype(inclFiles);
                
                index = repmat(tmpdayindex, [length(newfiles), 1]);
                index = [index newfiles];
                allindex = [allindex; index];
                clear index tmp*
            end
        end
    else
        dayindex = [metadata.AnimalID metadata.Date];
        allindex = [];
        for d = 1:height(metadata)
            tmpdayindex = dayindex(d,:);
            tmpfiles = cell2mat(metadata{d,'Files'});
            index = repmat(tmpdayindex, [length(tmpfiles), 1]);
            index = [index tmpfiles];
            
            allindex = [allindex; index];
            clear tmp* index
        end
    end
    
    %%% define colors
    colors.gamma.pre = [166 206 227]./255;
    colors.gamma.post = [31 120 180]./255;
    colors.random.pre = [178 223 128]./255; 
    colors.random.post = [51 160 44]./255;
    colors.behavior.RZ = [];
    colors.behavior.AZ = [];
    colors.behavior.half(1,:) = hex2rgb('#a62d4c');
    colors.behavior.half(2,:) = hex2rgb('#903a81');
    colors.celltypes.PYR = [];
    colors.celltypes.IN = [];
    %colors.regions.CA1 = hex2rgb('#9D9BCC');
    colors.regions.CA1 = hex2rgb('#035243');
    colors.regions.CA3 = hex2rgb('#5C2D9A');
    
    params.colors = colors; 
    
    
end
end

