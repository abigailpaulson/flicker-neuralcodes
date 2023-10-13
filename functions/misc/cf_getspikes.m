function [units, unitInfo, spikingvect] = cf_getspikes(anprocesseddatadir, ...
    dirs, params, index, brainReg)
%cf_getspikes
%   modified from loadsingleunitdata.m
%ALP 12/19/2022

%% check stuff
if ~iscell(brainReg)
    brainReg = {brainReg}; 
end

%% load cluster files
[units, unitInfo] = cf_loadclusterfiles(anprocesseddatadir, index, ...
    brainReg);

if isempty(units)
    return
end

%% check length unit vector
if length(units) ~= index(end,3)
    error('check units length, not equal to file #')
end

%% make spiking vector
%ALP 10/06/22 

for i = 1:length(units) %files
    tmpTimes = [];
    tmpIDs = [];
    iSort = []; 
    for j = 1:length(units(i).data) %# of units
        tmpTimes = [tmpTimes; units(i).data(j).spikeTimes];
        tmpIDs = [tmpIDs; j.*ones(length(units(i).data(j).spikeTimes),1)];
    end
    
    [~, iSort] = sort(tmpTimes);
    spikingvect(i).spikeTimes = tmpTimes(iSort);
    spikingvect(i).spikeIDs = tmpIDs(iSort);
end

%% get number of cells in each brain regions
for br = 1:length(brainReg)
    nCells(br) = sum(strcmp({unitInfo.brainReg}, brainReg{br}));
end

%% add cell type information 
celltype = [];
for br = 1:length(brainReg)
    %load cell type info
    if nCells(br) > 0 %don't try to load if there are no cells on this day
        cellTypeFileName = [dirs.celltypedir, brainReg{br}, '\cellTypeProps\cellTypeProps_', ...
            params.iden, num2str(index(1,1)), '_', num2str(index(1,2)), '.mat'];
        if exist(cellTypeFileName, 'file')
            load(cellTypeFileName, 'cellTypeProps')
            celltypeappend = cellTypeProps{index(1,1)}{index(1,2)}.cellTypeInfo.celltype;
            if sum(isnan(cellTypeProps{index(1,1)}{index(1,2)}.cellTypeInfo.celltype)) == length(cellTypeProps{index(1,1)}{index(1,2)}.cellTypeInfo.celltype)
                if length(cellTypeProps{index(1,1)}{index(1,2)}.cellTypeInfo.celltype) > length(cellTypeProps{index(1,1)}{index(1,2)}.meanFR)
                    celltypeappend = cellTypeProps{index(1,1)}{index(1,2)}.cellTypeInfo.celltype(1:end-1);
                end
            end
        else
            celltypeappend = NaN(nCells(br),1); %for cases where the spikes exist but cell types haven't been run yet
        end
    else
        celltypeappend = [];
    end
    celltype = [celltype; celltypeappend];
end

pyrinds = celltype == 1;
ininds = celltype == 2;
celltype = num2cell(celltype); 
celltype(pyrinds) = {'PYR'}; 
celltype(ininds) = {'IN'}; 

[unitInfo(1:length(celltype)).cellType] = celltype{:}; 


end

