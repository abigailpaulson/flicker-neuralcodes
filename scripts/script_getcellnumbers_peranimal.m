%% get cell number information for different analyses

clear all; close all

%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

basedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\'; 
isPost = [metadata.FlickerDay] > 5;
PostDays = metadata.Date(isPost);
isIncluded = ismember(allindex(:,2), PostDays);

allindex = allindex(isIncluded,:);
brainReg = {'CA3', 'CA1'};
gnames = {'random', 'gamma'};

isGamma = strcmp(metadata.Groups, 'gamma');
groupanimals.gamma = unique(metadata.AnimalID(isGamma));
groupanimals.random = unique(metadata.AnimalID(~isGamma)); 

%% cell types

for br = 1:2
    for g = 1:2
        ganimals = groupanimals.(gnames{g});
        cellNum = zeros(length(ganimals),2);
        for a = 1:length(ganimals)
            anday = unique(allindex(allindex(:,1) == ganimals(a),2));
            for d = 1:length(anday)
                dir = [basedir, 'celltypes\', brainReg{br}, '\cellTypeProps\'];
                filename = ['cellTypeProps_A', num2str(ganimals(a)), '_', num2str(anday(d)), '.mat'];
                load([dir, filename])
                CellType = cellTypeProps{ganimals(a)}{anday(d)};
                cellNum(a,1) = cellNum(a,1)+length(CellType.cellTypeInfo.PYRidx);
                cellNum(a,2) = cellNum(a,2)+length(CellType.cellTypeInfo.INidx);
            end
        end
        disp(['------- ', gnames{g}, ' ---------'])
        disp(['There are ', num2str(min(cellNum(:,1))), ' - ', num2str(max(cellNum(:,1))), ' PYR cells in ', brainReg{br}])
        disp(['There are ', num2str(min(cellNum(:,2))), ' - ', num2str(max(cellNum(:,2))), ' IN cells in ', brainReg{br}])
    end
end

%% place cells

for g = 1:2
    ganimals = groupanimals.(gnames{g});
    cellNum = zeros(length(ganimals));
    for a = 1:length(ganimals)
        anday = unique(allindex(allindex(:,1) == ganimals(a),2));
        for d = 1:length(anday)
            dir = [basedir, 'placecells\'];
            filename = ['spatialmaps_all_500_full_', num2str(anday(d)), '.mat'];
            load([dir, filename])

            cellNum(a,1) = cellNum(a,1)+sum(ratemaps.includeCell);
        end
    end
    disp(['------- ', gnames{g}, ' ---------'])
    disp(['There are ', num2str(min(cellNum(:,1))), ' - ', num2str(max(cellNum(:,1))), ' place cells'])
end

%% number of cells in theta decoding analyses
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
seqfilename = ['thetaseqdecoding_alldays_alltrials_230718.mat'];
load([seqdir, seqfilename])
sigDays = unique(AllData.day(AllData.significantSeq == 1));

for g = 1:2
    ganimals = groupanimals.(gnames{g});
    cellNum = zeros(length(ganimals),1);
    for a = 1:length(ganimals)
        anday = unique(allindex(allindex(:,1) == ganimals(a),2));
        isSigDay = ismember(anday, sigDays);
        anday = anday(isSigDay);
        if isempty(anday)
            cellNum(a,1) = NaN;
        end
        for d = 1:length(anday)
            dir = [basedir, 'placecells\'];
            filename = ['spatialmaps_all_500_full_', num2str(anday(d)), '.mat'];
            load([dir, filename])

            cellNum(a,1) = cellNum(a,1)+sum(ratemaps.includeCell);
        end
    end
    disp(['------- ', gnames{g}, ' ---------'])
    disp(['There are ', num2str(min(cellNum(:,1))), ' - ', num2str(max(cellNum(:,1))), ' place cells included in theta seq decoding'])
end

%% number of cells in ripple decoding analysis
filedir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_ripples\'];
filename = ['rippledecoding_halfratiozones_360_PC_AllData.mat'];
load([filedir, filename])
inclRipples = AllData.includeDay == 1 & ~isnan(AllData.rippleposition) & ~isnan(AllData.rippleRatio) & ~isnan(AllData.rippleDecPos) & ~isnan(AllData.rippleRelPos) & strcmp(AllData.timepoint, 'post');
sigDays = unique(AllData.day(inclRipples));

for g = 1:2
    ganimals = groupanimals.(gnames{g});
    cellNum = zeros(length(ganimals),1);
    for a = 1:length(ganimals)
        anday = unique(allindex(allindex(:,1) == ganimals(a),2));
        isSigDay = ismember(anday, sigDays);
        anday = anday(isSigDay);
        if isempty(anday)
            cellNum(a,1) = NaN;
        end
        for d = 1:length(anday)
            dir = [basedir, 'placecells\'];
            filename = ['spatialmaps_all_500_full_', num2str(anday(d)), '.mat'];
            load([dir, filename])

            cellNum(a,1) = cellNum(a,1)+sum(ratemaps.includeCell);
        end
    end
    disp(['------- ', gnames{g}, ' ---------'])
    disp(['There are ', num2str(min(cellNum(:,1))), ' - ', num2str(max(cellNum(:,1))), ' place cells included in ripple decoding decoding'])
end

%% number of cells in unengaged/engaged analysis
seqdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\']; 
filename = 'ThetaSeq_Engaged_Unengaged_180Decoding_SupplementaryFig.mat';
load([seqdir, filename])

sigDays = DayData.days;
for g = 1:2
    ganimals = groupanimals.(gnames{g});
    cellNum = zeros(length(ganimals),1);
    for a = 1:length(ganimals)
        anday = unique(allindex(allindex(:,1) == ganimals(a),2));
        isSigDay = ismember(anday, sigDays);
        anday = anday(isSigDay);
        if isempty(anday)
            cellNum(a,1) = NaN;
        end
        for d = 1:length(anday)
            dir = [basedir, 'placecells\'];
            filename = ['spatialmaps_all_500_full_', num2str(anday(d)), '.mat'];
            load([dir, filename])

            cellNum(a,1) = cellNum(a,1)+sum(ratemaps.includeCell);
        end
    end
    disp(['------- ', gnames{g}, ' ---------'])
    disp(['There are ', num2str(min(cellNum(:,1))), ' - ', num2str(max(cellNum(:,1))), ' place cells included in unengaged/engaged theta seq'])
end

