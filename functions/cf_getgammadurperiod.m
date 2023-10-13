function gammapower_z = cf_getgammadurperiod(processeddatadir, index, windows, samprate, brainreg)
%cf_getgammadurperiod
%
%ALP 3/18/2023

%load ZScore information
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\experimentinfo\zscore_lowgamma_alldays_', brainreg, '.mat'])

dayind = find(ZScore.day == index(1,2));

%load gamma
load([processeddatadir, brainreg, '\' num2str(ZScore.ripplechannel(dayind)), '\EEG\lowgamma', num2str(index(1,3)), '.mat'])

%check sampling rates are the same
if ~(lowgamma{index(1,1)}{index(1,2)}{index(1,3)}.samprate == samprate)
    error('sampling rates are not the same')
end
dat = lowgamma{index(1,1)}{index(1,2)}{index(1,3)}.data(:,3);

% loop over periods
gammapower_z = NaN(size(windows,1), 1);
for i = 1:size(windows,1)
    if windows(i,3) > length(dat)
        continue
    end
    avg_lfp = mean(dat(windows(i,1):windows(i,3))); 
    gammapower_z(i) = (avg_lfp - ZScore.lowgamma_mean(dayind))/ZScore.lowgamma_std(dayind);
end

end