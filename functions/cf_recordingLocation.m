%script_probePosition
%
% ALP 5/5/2022

clear;
%% experiment info
[params, dirs, indices, groups] = projectInfo('chronicflicker_annulartrack', ...
    'datethresh', 200800);

%%% common stuff
dayindex = indices.dayindex;
allindex = indices.allindex;
grp = groups.animals;
fDay = indices.flickerDay;
groupdayinds = horzcat(ismember(dayindex(:,1), grp{1}), ...
    ismember(dayindex(:,1), grp{2}));
recType = indices.recType;
brainReg = 'CA3';

g1before = groupdayinds(:,1) & fDay < 5;
g1after = groupdayinds(:,1) & fDay > 5;
g2before = groupdayinds(:,2) & fDay < 5;
g2after = groupdayinds(:,2) & fDay > 5;
concatinds_40randompre_40randompost = [g1before g2before g1after g2after]; %PRE THEN POST, 40HZ THEN RANDOM

%%% directories
analysisdir = [dirs.projectdir, 'probeplacement\'];
savefigdir = [analysisdir, 'fig_poweroverprobe\'];

%% define probe info
[x,y] = getchannelmap('Takahashi', 100); 

chmapfinal = [25,29,30,1,2,3,27,6,31,8,4,24,5,22,21,7,18,20,23,13,15,...
    9,17,11,19,16,32,26,28,14,12,10,33,39,37,51,53,55,52,50,56,48,54,46,...
    49,44,43,47,58,42,45,38,59,34,57,61,41,60,40,36,35,64,63,62];
channels = 0:63; 

linearizedch{1} = [4 15 2 22 9 23 18 1 19 11 3 10 0 17 13 7 16 29 6 27 30 8 28 20 25 5 14 24 21 31 26 12];
linearizedch{2} = [48 59 54 44 61 45 40 52 41 62 53 60 50 57 63 47 56 36 46 34 55 33 38 42 35 49 58 32 43 39 51 37];

%% plot to see how it looks 
plot(x, y, '.', 'MarkerEdgeColor', [0.85 0.85 0.85], 'MarkerSize', 20)

%% loop over days
for d = 1:size(dayindex,1)
    disp(['getting LFP power for ', num2str(dayindex(d,2))])
    anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(dayindex(d,1)), ...
        '_', num2str(dayindex(d,2)), '\', brainReg, '\']; 
    
    index = allindex(allindex(:,2) == dayindex(d,2),:); 
    lfppowerfile = [analysisdir, 'dat_poweroverprobe\' brainReg, '_lfppower_', num2str(dayindex(d,2)), '.mat']; 
    
    if ~exist(lfppowerfile, 'file')
         temp = getprobeplacementinfo(index, anprocesseddatadir, channels); 
         lfpPower.data = temp; 
         lfpPower.brainReg = brainReg; 
         lfpPower.index = index;
         lfpPower.channels = channels; 
         save(lfppowerfile, 'lfpPower')
    else
        load(lfppowerfile)
        temp = lfpPower.data;
    end
    
    power(d) = temp;
    clear temp
end

%% adjust low gamma power for bad channel
for d = 1:size(dayindex,1)
    tempgamma = power(d).lowgamma; 
    
    sortedLG = sort(tempgamma); 
    
    if sortedLG(end) > 1.5*sortedLG(end-1)
        tempgamma(tempgamma == max(tempgamma)) = NaN; 
    end
    power(d).lowgamma = tempgamma; 
end

%% plot!
nDays = sum(concatinds_40randompre_40randompost,1);
lfps = fieldnames(power(1));
hem = {'R', 'L'}; 
daytype = {'0/1', '9/10'}; 
daytype2 = {'pre', 'post'};
colors{1} = cbrewer('seq', 'Reds', length(channels));
colors{2} = cbrewer('seq', 'Blues', length(channels));
lineplotcoords = [x(chmapfinal == 7) x(chmapfinal == 26+32)];

gind = 1;
for d = 1:2
    figure(d)
    hold on
    set(gcf, 'color', 'w')
    set(gcf, 'Position', [126 285 1687 552])
    sgtitle(['ripple power - days ', daytype{d}, ' - ', brainReg])
    annotation('textbox', [0.05 0.8 0.1 0.1], 'String', '40Hz')
    annotation('textbox', [0.05 0.4 0.1 0.1], 'String', 'Random')
    
    figure(d+2)
    hold on
    set(gcf, 'color', 'w')
    set(gcf, 'Position', [126 285 1687 552])
    sgtitle(['low gamma - days ', daytype{d}, ' - ', brainReg])
    annotation('textbox', [0.05 0.8 0.1 0.1], 'String', '40Hz')
    annotation('textbox', [0.05 0.4 0.1 0.1], 'String', 'Random')
    
    iSub = 1; 
    for g = 1:2
        inclData = power(concatinds_40randompre_40randompost(:,gind));
        grouphemis = indices.hemisphere(concatinds_40randompre_40randompost(:,gind)); 
        if g == 1
            nsubplots = max(nDays(gind:gind+1)); 
        end
        
        for i = 1:length(inclData)
            zPower = structfun(@(x) (x-nanmean(x))/nanstd(x), inclData(i), 'UniformOutput', false);
            fracPower = structfun(@(x) x./max(x), inclData(i), 'UniformOutput', false);
            normPower = structfun(@(x) normalize(x, 'range'), inclData(i), 'UniformOutput', false);
            for f = 1:numel(lfps)
                if f == 1
                    iFig = d;
                else
                    iFig = d+2; 
                end

                %size = 40+20*zPower.(lfps{f});
%                size = 50.*normPower.(lfps{f});
                size = 30.*fracPower.(lfps{f}); 
                size = size(chmapfinal);
                size(size<0) = 1;
                size(size == 0) = 1;
                
                [~, iSort] = sort(size);
                plotcolor(iSort,:) = colors{f};
                
                [~, cR] = max(inclData(i).(lfps{f}));
                iCR = find(cR == chmapfinal);
                yCoordCR(f) = y(iCR); 
                
                plotcolor(iCR,:) = [0 153/255 51/255];
                %     plotcolor(iCR,:) = [0 204/255 0];
                
                figure(iFig)
                subplot(2,nsubplots,iSub); hold on
                scatter(x, y-yCoordCR(1), size, plotcolor, 'filled')
                plot(lineplotcoords, [200 200], 'k-')
                ylim([-200 200])
                xticks([]); yticks([]); axis off
                title(hem{grouphemis(i)})
                
            end
            iSub = iSub+1; 
        end
        gind = gind+1;
    end
    
    saveas(figure(d), [savefigdir, brainReg, '_', lfps{1}, '_', daytype2{d}], 'fig')
    saveas(figure(d), [savefigdir, brainReg, '_', lfps{1}, '_', daytype2{d}], 'png')
    
    saveas(figure(d+2), [savefigdir, brainReg, '_', lfps{2}, '_', daytype2{d}], 'fig')
    saveas(figure(d+2), [savefigdir, brainReg, '_', lfps{2}, '_', daytype2{d}], 'png')
    
end






%% plot
% zPower = structfun(@(x) x/max(x), power, 'UniformOutput', false);
% size = 100.*zPower.ripple;
figure
hold on
set(gcf, 'color', 'w')


lfps = fieldnames(power(1));
plottitles = {'ripple power', 'gamma power'};
CHcolor = [1 0 0; 0 0 1];

ploti = 1;
for d = 1:4
    zPower = structfun(@(x) (x-nanmean(x))/nanstd(x), power(d), 'UniformOutput', false);
    fracPower = structfun(@(x) x./max(x), power(d), 'UniformOutput', false);
    normPower = structfun(@(x) normalize(x, 'range'), power(d), 'UniformOutput', false);
    for f = 1:numel(lfps)
        %size = 40+20*zPower.(lfps{f});
        size = 50.*normPower.(lfps{f});
        size = size(chmapfinal);
        size(size<0) = 1;
        size(size == 0) = 1;
        
        [~, iSort] = sort(size);
        plotcolor(iSort,:) = colors{f};
        
        [~, cR] = max(power(d).(lfps{f}));
        iCR = find(cR == chmapfinal);
        
        %     colors = repmat([0.85 0.85 0.85], length(channels),1);
        plotcolor(iCR,:) = [0 153/255 51/255];
        %     plotcolor(iCR,:) = [0 204/255 0];
        %     plotcolor(iCR,:) = [204/255 0 204/255];
        
        
        subplot(2,4,ploti)
        hold on
        scatter(x, y, size, plotcolor, 'filled')
        plot(lineplotcoords, [200 200], 'k-')
        ylim([-20 200])
        xticks([])
        yticks([])
        axis off
        title(plottitles{f})
        
        ploti = ploti+4; 
    end
    ploti = ploti-7; 
end

figure
hold on
for d = 1:2
    fracPower = structfun(@(x) x./max(x), power(d), 'UniformOutput', false);
    
    [~, numSort] = sort(chmapfinal); 
    mat = fracPower.ripple; 
    mat = mat(chmapfinal); 
    
    mat = [mat' x y']; 
    mat(33:64, 2:3) = mat(33:64,2:3)+100; 
    
    [~, imaxval] = max(mat(:,1));
    
    maxcoord = mat(imaxval, 2:3); 
    
    for i = 1:length(channels)
%         tmpx = [maxcoord(1) mat(i,2)];
%         tmpy = [maxcoord(2) mat(i,3)];
        tmp1 = maxcoord;
        tmp2 = mat(i,2:3); 
        
        dist(i) = norm(tmp1 - tmp2);
        if tmp1(2) > tmp2(2)
            dist(i) = -dist(i); 
        end
    end
    
    pwrDist = [mat(:,1) dist']; 
    
    [~, iSortDist] = sort(pwrDist(:,2));
    pwrDist = pwrDist(iSortDist,:);
    
   
    plot(pwrDist(:,1), pwrDist(:,2));

end
% % 
% % 
% 




















