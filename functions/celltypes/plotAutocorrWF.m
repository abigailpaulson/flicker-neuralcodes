function plotAutocorrWF(dayindex, autocorr, autocorr_mean, wf, unit, savedfiguresdir, iden)
%this function plots the autocorrelogram and the corresponding waveform for
%individual clusters 
%SP 4.27.18
%updated ALP 1/16/2020

binsize = 5;
time = (binsize*(length(autocorr)-1))/2;
figure; clf;
subplot(1,2,1)
bar(-time:binsize:time, autocorr); hold on;
centerofmass = autocorr_mean*5;
plot(centerofmass,2,'*g')
title('Autocorrelogram with center of mass (green)')
xlim([-time time])
hold on;
subplot(1,2,2)
plot(wf)
title(['Waveform - cluster ' num2str(unit) ' ', iden, num2str(dayindex(1)) num2str(dayindex(2))]);
xlim([0 length(wf)])

%savefigures
datadir = [savedfiguresdir 'ExampleAutocorrelograms\' num2str(dayindex(2)), '\'];
savefigALP(datadir, 'meanACs_examples', 'png');
end