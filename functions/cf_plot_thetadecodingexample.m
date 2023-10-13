function cf_plot_thetadecodingexample(dirs, params, allindex, metadata, fh, ax)
%cf_plot_thetadecodingexample
%
%ALP 1/18/23

dayindex = unique(allindex(:,1:2), 'rows'); 
exampleday = [30 200917]; %from visually inspecting plots 

iDay = find(dayindex(:,2) == exampleday(2));

%%% directory
dir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'ALL_decoding_thetaSeq_180_PC_alltrials.mat';

%%% set edges for plotting
timeEdges = -0.2:0.02:0.2; %copied from params section of cf_decoding_thetaseq. ideally will resave structure with params files 
posEdges = -81:3:99; %same as above
adjustedEdges = -180/2:3:180/2; %180 for now


%%% load data
load([dir, filename])
exSeq = allData(iDay).mnDecodedSeq; 
exTheta = allData(iDay).mnThetaTrace; 
exTheta = exTheta./max(exTheta)*10; %normalize into 10 deg for plotting


iPanel = 1;
%%% figures
figure(fh)
axes(ax{iPanel})
cla
hold on
imagesc(timeEdges(1:end-1)+0.01, adjustedEdges(1:end-1)+1.5, exSeq) %shift cuz imagesc plots the center of the bin over the edge
plot(timeEdges(1:end-1), zeros(1,size(exSeq,2)), 'w--')
plot(zeros(1,size(exSeq,1)), adjustedEdges(1:end-1), 'w--')
plot(timeEdges(1:end-1), 30+exTheta, 'w-')
xlim([-0.08 0.08])
ylim([-45 45])
yticks([-45 0 45])
xticks([-0.08 0 0.08])
ylabel('relative decoded position (deg)')
xlabel('time from theta trough (s)')
c = colorbar; 
c.Ticks = [round(c.Limits(1),3) floor(c.Limits(2)*1000)/1000];
c.Label.String = 'Probability'; 
iPanel = iPanel+1;

figure(fh)
axes(ax{iPanel})
cla
hold on
imagesc(timeEdges(1:end-1)+0.01, adjustedEdges(1:end-1)+1.5, exSeq)
plot(timeEdges(1:end-1), zeros(1,size(exSeq,2)), 'w--')
plot(zeros(1,size(exSeq,1)), adjustedEdges(1:end-1), 'w--')
xlim([-0.08 0.08])
ylim([-45 45])
yticks([-45 0 45])
patch([-0.08 -0.08 0 0],[0 45 45 0], 'w')
patch([0 0 0.08 0.08],[-45 0 0 -45], 'w')
ylabel('relative decoded position (deg)')
xlabel('time from theta trough (s)')
xticks([-0.08 0 0.08])
c = colorbar;
c.Ticks = [];




end

