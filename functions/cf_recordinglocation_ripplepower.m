function cf_recordinglocation_ripplepower(dayindex, metadata, brainReg, figdir)
%cf_recordinglocation_ripple power
%
%   ALP 1/16/2023

%% helpful group vectors
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 

concatinds_40randompre_40randompost = [isGamma&isPre isRandom&isPre isGamma&isPost isRandom&isPost];

%% define probe information

[x,y] = cf_getchannelmap('Takahashi', 100); 

chmapfinal = [25,29,30,1,2,3,27,6,31,8,4,24,5,22,21,7,18,20,23,13,15,...
    9,17,11,19,16,32,26,28,14,12,10,33,39,37,51,53,55,52,50,56,48,54,46,...
    49,44,43,47,58,42,45,38,59,34,57,61,41,60,40,36,35,64,63,62];
channels = 0:63; 

linearizedch{1} = [4 15 2 22 9 23 18 1 19 11 3 10 0 17 13 7 16 29 6 27 30 8 28 20 25 5 14 24 21 31 26 12];
linearizedch{2} = [48 59 54 44 61 45 40 52 41 62 53 60 50 57 63 47 56 36 46 34 55 33 38 42 35 49 58 32 43 39 51 37];

%% loop over days to load the power structure
analysisdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\probeplacement\']; %previously saved

for d = 1:size(dayindex,1)
    lfppowerfile = [analysisdir, 'dat_poweroverprobe\' brainReg, '_lfppower_', num2str(dayindex(d,2)), '.mat'];
    load(lfppowerfile)
    temp = lfpPower.data;
    
    power(d) = temp;
    clear temp
end

%% plot!
nDays = sum(concatinds_40randompre_40randompost,1);
lfps = fieldnames(power(1));
hem = {'R', 'L'}; 
daytype = {'0/1', '9/10'}; 
daytype2 = {'pre', 'post'};
colors{1} = cbrewer('seq', 'Purples', length(channels));
colors{2} = cbrewer('seq', 'Blues', length(channels));
lineplotcoords = [x(chmapfinal == 7) x(chmapfinal == 26+32)];

gind = 3;
for d = 2
    figure
    hold on
    set(gcf, 'units', 'inches')
    set(gcf, 'color', 'w')
    set(gcf, 'Position', [4.5729 3.9271 7 2.5])
    sgtitle(['ripple power - days ', daytype{d}, ' - ', brainReg])
    annotation('textbox', [0.2 0.85 0.1 0.1], 'String', '40Hz')
    annotation('textbox', [0.2 0.45 0.1 0.1], 'String', 'Random')
    
    iSub = 1;
    for g = 1:2
        inclData = power(concatinds_40randompre_40randompost(:,gind));
        grouphemis = metadata.Hemisphere(concatinds_40randompre_40randompost(:,gind));
        if g == 1
            nsubplots = max(nDays(gind:gind+1));
        end
        anIDs = dayindex(concatinds_40randompre_40randompost(:,gind),1);
        
        for i = 1:length(inclData)
            zPower = structfun(@(x) (x-nanmean(x))/nanstd(x), inclData(i), 'UniformOutput', false);
            fracPower = structfun(@(x) x./max(x), inclData(i), 'UniformOutput', false);
            normPower = structfun(@(x) normalize(x, 'range'), inclData(i), 'UniformOutput', false);
            
            for f = 1
                
                %size = 40+20*zPower.(lfps{f});
                %                size = 50.*normPower.(lfps{f});
                dotsize = 5.*fracPower.(lfps{f});
                dotsize = dotsize(chmapfinal);
                dotsize(dotsize<0) = 1;
                dotsize(dotsize == 0) = 1;
                
                [~, iSort] = sort(dotsize);
                plotcolor(iSort,:) = colors{f};
                
                [~, cR] = max(inclData(i).(lfps{f}));
                iCR = find(cR == chmapfinal);
                yCoordCR(f) = y(iCR);
                
                plotcolor(iCR,:) = [1 0 0];
%                 plotcolor(iCR,:) = [0 153/255 51/255];
                %     plotcolor(iCR,:) = [0 204/255 0];
                
                subplot(2,nsubplots,iSub); hold on
                scatter(x, y-yCoordCR(1), dotsize, plotcolor, 'filled')
                plot(lineplotcoords, [200 200], 'k-')
                ylim([-200 200])
                xticks([]); yticks([]); axis off
                title([num2str(anIDs(i)), ' - ' grouphemis{i}])
                 
            end
            iSub = iSub+1;
        end
        gind = gind+1;
    end
end

makefigurepretty(gcf)
savefigALP(figdir, ['probeposition_', daytype2{d}, '_', brainReg], 'filetype', 'pdf')

end

