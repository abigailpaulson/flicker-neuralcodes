function cf_plot_placecellproperties(dirs, params, allindex, metadata, fh, ax, statfid, panelL, tablefilename)
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

%% plot place cell properties for distance to reward

%load full position directory
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';
mapfilename = 'ALLDAYS_spatialmaps_preflicker_0_d2r_';

load([ratemapdir, mapfilename])

fnames = fieldnames(groupData.gamma.pre);

yy = {'spatial information (bits/spike)', 'sparsity', 'peak position (deg)', 'peak firing rate (hz)', 'mean firing rate (hz)'};
xxl = {[0 1], [0 1], [-81 99], [0 10], [0 10]};
analysisnames = {'spatial information', 'sparsity', 'peak position', 'peak firing rate', 'mean firing rate'};
yUnits = {'bits/spike', 'sparsity', 'deg', 'Hz', 'Hz'};

figure(fh)
iPlot = 1;

for f = 2:length(fnames)-4
    axes(ax{iPlot})
    hold on
    
    vdat = []; cmat = []; stat = [];
    for d = 2
        datforstats = [];
        for g = 1:2
            nm = [gnames{g}, dnames{d}];
            tmpdat = groupData.(gnames{g}).(dnames{d}).(fnames{f}); 
            incldat =  groupData.(gnames{g}).(dnames{d}).isPC;
            vdat.(nm) = tmpdat(incldat);
            datforstats{g} = tmpdat(incldat);
            cmat = [cmat; groupcolors.(gnames{g}).(dnames{d})];
        end
    end
%     vdat.xlabels = {'pre', 'post'};
    violinplot_half(vdat, [], 'ViolinColorMat', cmat, 'BoxWidth', 0.018, 'MedianSize', 20, 'ShowData', false);
    ylabel(yy{f-1})
    ylim(xxl{f-1})
    title([fnames{f}])
    cf_stats2txt2(datforstats, statfid, panelL{iPlot}, 'place cells', yUnits{f-1}, analysisnames{f-1}, gnames, tablefilename)

    
    
    iPlot = iPlot+1;
end




end

