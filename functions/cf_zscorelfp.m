function cf_zscorelfp(dirs, brainReg, allindex, filename)
%cf_zscorelfp
%     zscore the lfp file that is input as filename
%       default is to use the best ripple channel for each day as the
%       loading channel, from the region specified in brainReg
%ALP 3/18/2023

ripplechdir = ['\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\', brainReg, '\']; 
dayindex = unique(allindex(:,1:2), 'rows');

for d = 1:size(dayindex,1)
    files = allindex(allindex(:,2) == dayindex(d,2),3);
    animal = dayindex(d,1);
    day = dayindex(d,2);
    
    %load best ripple channel 
    load([ripplechdir, 'bestChannel_', num2str(dayindex(d,2)), '.mat'])
    ripplechannel = bestRippleChan.all;
    
    tmpLFP = [];
    for f = 1:length(files)
        %%% load LFP
        dat = load([dirs.processeddatadir, 'A', num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '\', brainReg, '\' num2str(ripplechannel), '\EEG\', filename, num2str(files(f)), '.mat']);
        tmpLFP = [tmpLFP; dat.(filename){dayindex(d,1)}{dayindex(d,2)}{files(f)}.data(:,3)];
        
        clear dat
    end
    
    ZScore(d).animal = animal;
    ZScore(d).day = day;
    ZScore(d).files = files;
    ZScore(d).region = brainReg;
    ZScore(d).([filename,'_mean']) = mean(tmpLFP);
    ZScore(d).([filename, '_std']) = std(tmpLFP);
    ZScore(d).ripplechannel = ripplechannel;
     
    clear tmpLFP animal day ripplechannel
end

ZScore = struct2table(ZScore);

savedir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\experimentinfo\';
savefilename = ['zscore_', filename, '_alldays_', brainReg, '.mat'];

save([savedir, savefilename])

end

