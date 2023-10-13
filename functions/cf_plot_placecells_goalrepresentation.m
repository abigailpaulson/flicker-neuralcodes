function cf_plot_placecells_goalrepresentation(dirs, params, allindex, metadata, ...
    fh, ax, statfid, panelL, tablefilename)
% plot prop reward zone cells per day
%   only cells included in the theta cycle analysis
%ALP 1/21/23

%%% helpful vectors for plotting
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5;
% gID.gamma.pre = dayID(isGamma & isPre);
% gID.gamma.post = dayID(isGamma & isPost);
% gID.random.pre = dayID(isRandom & isPre);
% gID.random.post = dayID(isRandom & isPost);
gnames = {'random', 'gamma'};
dnames = {'pre', 'post'}; 

%%% load stuff
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
mapfilename = ['spatialmaps_preflicker_0_d2r_'];
PCfilename = ['spatialmaps_all_500_full_'];

dayindex = unique(allindex(:,1:2), 'rows');

for d = 1:size(dayindex,1)
    dayfilename = [mapfilename num2str(dayindex(d,2)), '.mat']; 
    dayPCs = [PCfilename, num2str(dayindex(d,2)), '.mat'];
    
    load([ratemapdir, dayfilename])
    tmpPC = load([ratemapdir, dayPCs]); 
    
    isPC = tmpPC.ratemaps.includeCell; 
    nPC = sum(isPC); 
    
    peakPos = ratemaps.peakpos(isPC); 
    isRRZPC = peakPos >= -18 & peakPos <9; 
    
    nRRZPC = sum(isRRZPC); 
    
    propRRZ_PC = nRRZPC/nPC; 
    
    plotdata(d).propRRZ_PC = propRRZ_PC; 
    
    clear tmpPC nPC ratemaps peakPos isRRZPc nRRZPC 
end

%%% load theta sequence data for inclusion variable
dir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ALL_decoding_thetaSeq_180_PC_alltrials.mat';
load([dir, filename], 'allData')
for d = 1:size(dayindex,1)
    if ~isempty(allData(d).sigDay)
        inclDay(d) = allData(d).sigDay;
    else
        inclDay(d) = 0;
    end
end

inclDay = logical(inclDay');

%%% get group data for plotting
groupdata.gamma.pre.propRRZ_PC = [plotdata(isGamma&isPre&inclDay).propRRZ_PC];
groupdata.gamma.post.propRRZ_PC = [plotdata(isGamma&isPost&inclDay).propRRZ_PC];
groupdata.random.pre.propRRZ_PC = [plotdata(isRandom&isPre&inclDay).propRRZ_PC];
groupdata.random.post.propRRZ_PC = [plotdata(isRandom&isPost&inclDay).propRRZ_PC];


%%% plot
iPanel = 1; 
figure(fh)
axes(ax{iPanel})
hold on;
xvect_bar = [1 2]; 
clear scatterdat plotdat_mn plotdat_sem
for d = 2
    for g = 1:2
        scatterdat{g} = groupdata.(gnames{g}).(dnames{d}).propRRZ_PC;
        plotdat_mn(g) = mean(scatterdat{g}, 'omitnan');
        plotdat_sem(g) = mean(scatterdat{g}, 'omitnan')./sqrt(length(scatterdat{g})); 
        colororder(g,:) = params.colors.(gnames{g}).(dnames{d}); 
    end
end
plotprettypoints(fh, xvect_bar, scatterdat, 'color', colororder)
b = bar(xvect_bar, plotdat_mn, 'FaceColor', 'flat');
b.CData = colororder;
b.FaceAlpha = 0.6;
errorbar2(xvect_bar, plotdat_mn, plotdat_sem, 0.2, 'k-', 'LineWidth', 0.75);
xticks([])
ylim([0 0.4])
yticks([0 0.1 0.2 0.3 0.4])
ylabel('proportion of place cells')
title({'reward related zone', 'place cells'})
cf_stats2txt2(scatterdat, statfid, panelL{iPanel}, 'days', 'proportion of place cells', 'RRZ_PC', gnames, tablefilename)






end

