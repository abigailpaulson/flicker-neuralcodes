function [units, unitInfo, spikingvect] = loadsingleunitdata(anprocesseddatadir, ...
    dirs, params, index, brainReg)
%loadsingleunitdata load single unit files
%inputs:
%   anproceseddatadir - base directory of the recording day
%   dirs - structure with directories pointing to celltypes and placecells
%   index - n x 3 matrix, an date file
%   brainReg - cell array of brain regions to load clusters files. 
%outputs:
%   units - structure with length(nRecs). units(i).data contains the
%   following subfields with data from each recording.
%       ID - vector of clusterID tags
%       maxChan - vector of channels of largest cluster amplitude
%       spikeInds - vector of spike indices
%       brainReg - vector of brainRegion
%   clusterMetrics - structure with length(nUnits) and fields as defined
%   in Singer Lab kilosort 2 pipeline     
%saved structures: none
%example call: 
%   [singleUnits, clusterMetrics] = loadsingleunitdata(datadir, ...
%       [22 200202], [1:3], {'CA1', 'CA3'}); 
%detailed description: n/a
%ALP 5/8/2020 chronicflicker-ephys
%change log (date, who, what, why)
%ALP 7/11/22 add distance to reward fields

%% check stuff
%check that brainReg is a cell array
if ~iscell(brainReg)
    brainReg = {brainReg}; 
end

%% load cluster files
[units, unitInfo] = loadClusterFiles_2reg(anprocesseddatadir, index, ...
    brainReg); 

if isempty(units)
    return
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

%% append cell type information
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

%% append place cell information
%this should work
PF = [];
MAP = [];
for br = 1:length(brainReg)
    if strcmp(brainReg{br}, 'CA1')
        placefieldfile = [dirs.placefielddir, ...
            'CA1\data\placefields_PYR_95threshold_randomshuffle_', ...
            params.iden, num2str(index(1,1)), '_', num2str(index(1,2)), '.mat'];
    else %otherwise should be CA3
        placefieldfile = [dirs.placefielddir, ...
            'CA3\data\placefields_PYR_95threshold_randomshuffle_', ...
            params.iden, num2str(index(1,1)), '_', num2str(index(1,2)), '.mat'];
    end
    
    if exist(placefieldfile, 'file')
        load(placefieldfile)
        
        PFappend = placefields{index(1,1)}{index(1,2)}.include;
        mapappend = placefields{index(1,1)}{index(1,2)}.ratenormocc;
        
        if size(PFappend,2) > 1
            PFappend = PFappend';
        end
        
    else
        PFappend = nan(nCells(br),1);
        mapappend = nan(nCells(br),180);
    end
    
    PF = [PF; PFappend];
    MAP = [MAP; mapappend];
end
MAP2 = num2cell(MAP, 2);
PF = num2cell(PF);
[unitInfo(1:length(PF)).isPlaceField] = PF{:};
[unitInfo(1:length(PF)).ratenormocc] = MAP2{:};

%% append distance to reward PC information

d2rdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\placefields\distance2reward\dat_d2r_ratenormocc\';
d2rfname = ['d2r_ratenormocc_', num2str(index(1,2)), '.mat'];

if exist([d2rdir, d2rfname], 'file')
    load([d2rdir, d2rfname]);
    
    if length(brainReg) < 2
        reg = brainReg{1};
        isReg = strcmp(spikerate_d2r.brainReg, reg);
        d2rMAP = spikerate_d2r.ratenormocc(isReg,:);
    else
        d2rMAP = spikerate_d2r.ratenormocc;
    end
    
    cdMAP = num2cell(d2rMAP,2);
    nc = size(d2rMAP,1);
    [unitInfo(1:nc).d2r_ratenormocc] = cdMAP{:};
else
    append = NaN(1, length(PF));
    append = num2cell(append);
    [unitInfo(1:length(PF)).d2r_ratenormocc] = append{:};
end



end

