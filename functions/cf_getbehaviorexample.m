function cf_getbehaviorexample(dirs, params, allindex, metadata, fh, ax)
%cf_getbehaviorexample
%ALP 12/5/22

%%% example day
% exday = 200822;
% exfile = 2;
% extrials = 4:6;
% 
% exday = 200831;
% exfile = 1;
% extrials = 8:10;

exday = 210211;
exfile = 1;
extrials = 14:16;

trackInfo = getTrackInfo_cflicker('chronicflicker_annulartrack');
RZ = trackInfo.rewardZone;

tempindex = allindex(allindex(:,2) == exday,:);
index = tempindex(tempindex(:,3) == exfile,:);

%%% directories
anprocesseddatadir = [dirs.processeddatadir, 'A', num2str(index(1,1)), '_', ...
    num2str(index(1,2)), '\'];
positiondir = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\'];

%%% load position info and behavior info
load([positiondir, 'behaviorInfo_A', num2str(index(1,1)), '_', num2str(index(1,2)), '.mat'])
load([positiondir, 'positionInfo_A', num2str(index(1,1)), '_', num2str(index(1,2)), '.mat'])
posFiles = [positionInfo.index];
posFiles = posFiles(3:3:end); 
iP = find(posFiles == index(1,3)); 

%%% load spiking data. use old PC data for now but I should update...
[singleunits, clusterMetrics, spikeVect] = cf_getspikes(anprocesseddatadir, dirs, params, ...
    index, params.brainReg);
brainReg = {clusterMetrics.brainReg};


%%% load ephys data
bestchdir = '\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\ripples\bestRippleChan\';
load([bestchdir, 'CA1\bestChannel_', num2str(index(1,2)), '.mat'])
ca1ch = bestRippleChan.all; 
load([bestchdir, 'CA3\bestChannel_', num2str(index(1,2)), '.mat'])
ca3ch = bestRippleChan.all; 

ca1 = load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\eeg', num2str(index(3)), '.mat']);
ca3 = load([anprocesseddatadir, 'CA3\', num2str(ca3ch), '\eeg', num2str(index(3)), '.mat']);
load([anprocesseddatadir, 'CA1\', num2str(ca1ch), '\ripples', num2str(index(3)), '.mat'])
LFPData(1,:) = ca1.eeg{index(1)}{index(2)}{index(3)}.data;
LFPData(2,:) = ca3.eeg{index(1)}{index(2)}{index(3)}.data;
SWRData = ripples{index(1)}{index(2)}{index(3)};
LFPtime = 0:1/2000:length(LFPData(1,:))/2000;

%%% load rate maps for sorting the spiking
load(['\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\spatialmaps_all_500_full_', num2str(index(1,2)), '.mat'])
peakpos = ratemaps.peakpos; 
region = ones(1,length(brainReg));
region(strcmp(brainReg,'CA3')) = 2;

%%% behavior data
theta = behaviorInfo(index(1,3)).theta;
time = round(positionInfo(iP).timeSmooth,5); 
speed = smoothdata(positionInfo(iP).speedSmooth, 'gaussian', 1*(1/0.02));
lickTimes = round(behaviorInfo(index(1,3)).lickTimes,5);
lickPos = behaviorInfo(index(1,3)).lickPos;
rewTimes = round(behaviorInfo(index(1,3)).rewardTimes,5);
[~,rewInd] = ismember(rewTimes, time); 
isRew = ind2vect(rewInd, length(time), 30);
[~, lickInd] = ismember(lickTimes, time);
isLick = ind2vect(lickInd, length(time), 1); 
rewPos = behaviorInfo(index(1,3)).rewardPos_d2r;

%%% divide up by trials
startind = [1; find(diff(theta)<-150)+1];
endind = startind-1;
endind = [endind(2:end); length(theta)];

%-------------------------------------------------------------------------
iT = extrials(2); 
t_time = time(startind(iT):endind(iT));
t_speed = speed(startind(iT):endind(iT));
t_theta = theta(startind(iT):endind(iT));
LFPstartind = lookup2(t_time(1), LFPtime);
LFPendind = lookup2(t_time(end), LFPtime);
t_ca1_eeg = LFPData(1,LFPstartind:LFPendind);
t_ca3_eeg = LFPData(2,LFPstartind:LFPendind);
% t_ca1_spiking = ;
% t_ca3_spiking = ; 
t_rewards = isRew(startind(iT):endind(iT));
t_licks = isLick(startind(iT):endind(iT));
%t_ripples = SWRData.midtime(logical(isExcluded(SWRData.midind, [LFPstartind LFPendind])));
t_ripples = SWRData.midtime;

isRZ1 = t_theta >= RZ(1,1) & t_theta < RZ(1,2);
isRZ2 = t_theta >= RZ(2,1) & t_theta < RZ(2,2);
isRZ = isRZ1 | isRZ2; 
t_time_RZ = t_time(isRZ);


figure(fh)
axes(ax{1})
cla
hold on
rnum = 6;
lfpinds = t_ripples(rnum)*2000 - 0.150*2000;
lfpendinds = t_ripples(rnum)*2000 + 0.150*2000;
plot(LFPtime(lfpinds:lfpendinds), LFPData(2,lfpinds:lfpendinds), 'k-')
plot(LFPtime(lfpinds:lfpendinds), 750+LFPData(1,lfpinds:lfpendinds), 'k-')
plot([LFPtime(lfpinds) LFPtime(lfpinds+0.05*2000)], [-200 -200], 'r-')
plot([LFPtime(lfpinds)+0.01 LFPtime(lfpinds)+0.01], [1000 1200], 'r-')
text(LFPtime(lfpinds), -300, '50 ms')
text(LFPtime(lfpinds), 1200, '200 uV')
ylim([-500 1300])
xlim([LFPtime(lfpinds) LFPtime(lfpendinds)])
xticks([])
yticks([])
title( 'SWR')

figure(fh)
axes(ax{2})
cla
hold on
lfpinds = LFPstartind + 20*2000;
lfpendinds = lfpinds+1.2*2000;
plot(LFPtime(lfpinds:lfpendinds), LFPData(2,lfpinds:lfpendinds), 'k-')
plot(LFPtime(lfpinds:lfpendinds), 700+LFPData(1,lfpinds:lfpendinds), 'k-')
plot([LFPtime(lfpinds) LFPtime(lfpinds+0.5*2000)], [-200 -200], 'r-')
plot([LFPtime(lfpinds)+0.01 LFPtime(lfpinds)+0.01], [1000 1200], 'r-')
text(LFPtime(lfpinds), -300, '500 ms')
text(LFPtime(lfpinds), 1200, '200 uV')
xlim([LFPtime(lfpinds) LFPtime(lfpendinds)])
xticks([])
yticks([])
title('Theta')


figure(fh)
axes(ax{3})
cla
hold on
set(gca, 'Units', 'Inches');
plot(t_time, 3+isRZ, 'b-')
plot(t_time, (t_theta./max(t_theta))*0.9+2.4, 'k-')
plot(t_time, (t_speed./max(t_speed))*0.65+1.3, 'k-')
plot(t_time, t_licks*0.35+0.6, 'k-')
plot(t_time, t_rewards*0.35-0.45, 'k-')
plot(t_ripples, -0.5*ones(1,length(t_ripples)), 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 3)
%%patch([LFPtime(lfpinds) LFPtime(lfpinds) LFPtime(lfpendinds) LFPtime(lfpendinds)], [0 3 3 0], 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
% plot(LFPtime(LFPstartind:LFPendind), t_ca1_eeg./max(t_ca1_eeg))
yticks([0 1 2])
text(t_time(end-300), 0.1, 'rewards', 'FontName', 'Arial')
text(t_time(end-200), 1.2, 'licks', 'FontName', 'Arial')
text(t_time(end-1000), 2.2, 'speed', 'FontName', 'Arial')
text(t_time(200), 3.2, 'position', 'FontName', 'Arial')
%xlabel('time (s)')

unitID = 1:length(region);
rastercolor = cbrewer('seq', 'BuPu', length(unitID)+20);
iC = 20;
yplot = -1;
%for i = 1:2 %try sorting all together ALP 10/10/23
    %isReg = region == i;
    isReg = true(1,length(peakpos));
    include = isReg ;
%     & ratemaps.includeCell;
    regPeaks = peakpos(isReg);
    regCells = unitID(isReg);
    [~,iSort] = sort(regPeaks);
    sortCells = regCells(iSort);
    
    %initialize colors
%     rastercolor = cbrewer('seq', 'YlGnBu', length(sortCells)+10);
%     iC = 10;
    for c = fliplr(sortCells)
        cellspikes = singleunits(index(1,3)).data(c).spikeTimes;
        isTrialSpikes = isExcluded(singleunits(index(1,3)).data(c).spikeTimes, [t_time(1) t_time(end)]);
        trialspikes = cellspikes(logical(isTrialSpikes)); 
        
        plot(trialspikes, yplot*ones(1,length(trialspikes)), '.', 'Color', params.colors.regions.('CA3'), 'MarkerSize', 0.005)
        yplot = yplot-0.1;
        iC = iC+1;
    end
%end

plot([t_time(50) t_time(50+10*50)], [yplot-0.2 yplot-0.2], 'r-')
text(t_time(100), yplot-0.75, '10 s', 'FontName', 'Arial')

xlim([t_time(1) t_time(end)])
ylim([yplot-0.5 4])
xticks([])


end

