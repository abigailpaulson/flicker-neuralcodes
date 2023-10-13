function cf_plot_placecellmaps(dirs, params, allindex, metadata, fh, ax)
%cf_plot_placecellmaps
%
%ALP 1/22/2023

%%% helpful vectors for groups, etc
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 
dayID = 1:height(metadata);
gID.gamma.pre = dayID(isGamma & isPre);
gID.gamma.post = dayID(isGamma & isPost);
gID.random.pre = dayID(isRandom & isPre);
gID.random.post = dayID(isRandom & isPost);
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 
groupcolors = params.colors;

%% 360 degree rate maps
%load full position directory
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
mapfilename = 'ALLDAYS_spatialmaps_preflicker_0_full_';

load([ratemapdir, mapfilename])

%%% plot rate maps
figure(fh)
hold on
d = 2;
for g = 1:2
    axes(ax{g})
    hold on
    
    tmp = groupData.(gnames{g}).(dnames{d}).peakpos;
    fr = groupData.(gnames{g}).(dnames{d}).ratemaps;
    incl = groupData.(gnames{g}).(dnames{d}).isPC;
    dat = groupData.(gnames{g}).(dnames{d}).ratemaps./mean(groupData.(gnames{g}).(dnames{d}).ratemaps,2,'omitnan');
    
    %%% how many animals contribute to PC?
    nAn = length(unique(groupData.(gnames{g}).(dnames{d}).animal(incl)));
    disp(['There are ', num2str(nAn), ' animals contributing to place cells for ', gnames{g}])
    
    [~, iSort] = sort(tmp);
    iSort = iSort(incl);
    plotDat = dat(iSort,:);
    
    plotRZ = ones(4, size(plotDat,1));
    plotRZ = reshape(RZ,1,4).*plotRZ'; 
    
    imagesc([0:2:358]+1, 1:size(plotDat,1), plotDat, [0.5 1.5])
    xlim([0 360])
    xticks([0 360])
    yticks([0 round(size(plotDat,1))])
    ylim([0 size(plotDat,1)])
    plot(plotRZ, repmat(1:size(plotDat,1), 4,1)', 'm-')
    title([gnames{g} ' ' dnames{d}])
    c = colorbar;
    c.Label.String = 'normalized firing rate'; 
    c.Ticks = [0.5 1 1.5];  
    colormap(flipud(cmocean('gray')));
    xlabel('position (deg)')
    ylabel('place cell count')
    
    clear tmp fr incl dat iSort plotDat
end
clear groupData allData params

% the below figure was moved to the main figure 3
%
%   ALP 7/13/2022
% %% 180 degree rate maps
% %load full position directory
% ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
% mapfilename = 'ALLDAYS_spatialmaps_preflicker_0_d2r_';
% 
% load([ratemapdir, mapfilename])
% 
% %%% plot rate maps
% figure(fh)
% hold on
% d = 2;
% for g = 1:2
%     axes(ax{g+2})
%     hold on
%     
%     tmp = groupData.(gnames{g}).(dnames{d}).peakpos;
%     fr = groupData.(gnames{g}).(dnames{d}).ratemaps;
%     incl = groupData.(gnames{g}).(dnames{d}).isPC;
%     dat = groupData.(gnames{g}).(dnames{d}).ratemaps./mean(groupData.(gnames{g}).(dnames{d}).ratemaps,2,'omitnan');
%     
%     [~, iSort] = sort(tmp);
%     iSort = iSort(incl);
%     plotDat = dat(iSort,:);
%     
%     plotRZ = ones(2, size(plotDat,1));
%     plotRZ = reshape(RZ,1,2).*plotRZ'; 
%     
%     imagesc([-81:3:97]+1, 1:size(plotDat), plotDat, [0.5 1.5])
%     xlim([-81 99])
%     xticks([-81 0 99])
%     yticks([0 round(size(plotDat,1))])
%     ylim([0 size(plotDat,1)])
%     plot(plotRZ, repmat(1:size(plotDat,1), 2,1)', 'w-')
%     title([gnames{g} ' ' dnames{d}])
%     c = colorbar;
%     c.Label.String = 'normalized firing rate'; 
%     c.Ticks = [0.5 1 1.5];  
%     xlabel('position (deg)')
%     ylabel('place cell count')
% end




end

