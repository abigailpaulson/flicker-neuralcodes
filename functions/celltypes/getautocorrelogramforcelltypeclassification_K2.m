function [autocorr_mean, autocorr] = getautocorrelogramforcelltypeclassification_K2(dayindex, dirs,...
    fullWFDiff, clusterInfo, iden, makenewfiles)
%SP 3.14.19
% Updated ALP 1/16/2020 for Kilosort2 datastructures 
disp(['Calculating autocorrelograms for sess ' num2str(dayindex(1)) num2str(dayindex(2))])

clusterdatadir = fullfile(dirs.processeddatadir, [iden num2str(dayindex(1)) '_' num2str(dayindex(2))], dirs.processedappend, dirs.clusterappend);
filename = [dirs.spikedatadir 'cellTypeDataAC_' iden num2str(dayindex(1)) '_' num2str(dayindex(2)) '.mat'];
if ~exist(filename) || makenewfiles
    [allrecProps.autocorr allrecProps.autocorr_mean] = calcAutocorr_K2_celltypes(clusterdatadir, clusterInfo);
    
    % plot the individual autocorrelograms and their respective WFs
    for unit = 1:length(allrecProps.autocorr)
        acorr = cell2mat(allrecProps.autocorr(unit));
        acorr_mean = cell2mat(allrecProps.autocorr_mean(unit));
        wf = fullWFDiff(unit,:);
        
        plotAutocorrWF(dayindex, acorr, acorr_mean, wf, unit, dirs.savedfiguresdir, iden);
    end
    close all
    
    cellTypeDataAC{dayindex(1)}{dayindex(2)} = allrecProps;
    save(filename, 'cellTypeDataAC');
else
    load(filename)
    allrecProps.autocorr = cellTypeDataAC{dayindex(1)}{dayindex(2)}.autocorr;
    allrecProps.autocorr_mean = cellTypeDataAC{dayindex(1)}{dayindex(2)}.autocorr_mean;
end

%% get new struct for output
autocorr = allrecProps.autocorr;
autocorr_mean = allrecProps.autocorr_mean;

%% save
cellTypeDataAC{dayindex(1)}{dayindex(2)} = allrecProps;
save(filename, 'cellTypeDataAC')

end