function cf_plot_currPosDecoding(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
%cf_plot_currPosDecoding
%   plot current position decoding for manuscript figures
%ALP 3/23/3023

iPlot = 1; 
%% load
savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\';
savefilename = ['Group_decoding_currpos_180', '_PC_test20']; 
load([savedatadir, savefilename])

%% stuff for plotting
params.posEdges = -81:3:99;
params.posBins = 3;
dType = '180';
nDeg = 180;
posName = 'theta_d2r';
params.RZ = [0 18];

RZxplot = sort(reshape(params.RZ, [numel(params.RZ),1]));
RZxplot = repmat(RZxplot, [1, length(params.posEdges(1:end-1))]);
RZyplot = repmat(params.posEdges(1:end-1), [numel(params.RZ),1]);


isGamma = strcmp(metadata.Groups, 'gamma');
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5;

%% get data of interest
inclDay = [groupData.include]';
cellPosEst = {groupData.mnEstPos};
gammaPosEst = cat(3,cellPosEst{isGamma&isPost&inclDay});
randomPosEst = cat(3,cellPosEst{isRandom&isPost&inclDay});

 %%% how many animals contribute to decoding?
nAnGamma = length(unique(metadata.AnimalID(isGamma&isPost&inclDay)));
nAnRandom = length(unique(metadata.AnimalID(isRandom&isPost&inclDay)));

disp(['There are ', num2str(nAnGamma), ' animals contributing to decoding for gamma'])
disp(['There are ', num2str(nAnRandom), ' animals contributing to decoding for random'])


%% figure('Position', [440 613 644 185])
axes(ax{iPlot})
hold on
title('probability/chance - post flicker sess - all sess')
imagesc(params.posEdges, params.posEdges, mean(randomPosEst,3,'omitnan').*(nDeg/params.posBins), [0.25 3])
plot(RZxplot', RZyplot', 'w-')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([min(params.posEdges) max(params.posEdges)])
xlabel('actual position (deg)')
ylabel('estimated position (deg)')
colorbar
title('random')

iPlot = iPlot+1;

axes(ax{iPlot})
hold on
imagesc(params.posEdges, params.posEdges, mean(gammaPosEst,3,'omitnan').*(nDeg/params.posBins), [0.25 3])
plot(RZxplot', RZyplot', 'w-')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([min(params.posEdges) max(params.posEdges)])
xlabel('actual position (deg)')
ylabel('estimated position (deg)')
title('gamma')
colorbar
iPlot = iPlot+1; 

axes(ax{iPlot})
hold on
gammaBinError = {groupData(isGamma&isPost&inclDay).binError};
gammaBinError = cell2mat(gammaBinError');
randomBinError = {groupData(isRandom&isPost&inclDay).binError};
randomBinError = cell2mat(randomBinError');

gammaMnError = mean(gammaBinError,1,'omitnan');
randomMnError = mean(randomBinError,1,'omitnan');
stdeGammaError = std(gammaBinError, [], 1, 'omitnan')./sqrt(sum(~isnan(gammaBinError(:,1))));
stdeRandomError = std(randomBinError, [], 1, 'omitnan')./sqrt(sum(~isnan(randomBinError(:,1))));

shadedErrorBar(params.posEdges(1:end-1), gammaMnError, stdeGammaError, {'Color', params.colors.gamma.post},1)
shadedErrorBar(params.posEdges(1:end-1), randomMnError, stdeRandomError, {'Color', params.colors.random.post},1)
plot(RZxplot', RZyplot', 'r--')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([-nDeg/8 nDeg/8])
xlabel('position (deg)')
ylabel('decoding error (deg)')
title(['decoding error - post flicker - random ', num2str(size(gammaBinError,1)), ' days - gamma ', num2str(size(gammaBinError,1))])





end

