function cf_plot_prospectivecoding_byspeed_example(dirs, params, allindex, metadata, fh, ax,  statfid, panelL, tablefilename)
%Figure 5 example
%ALP 7/20/23
%% set up colors, parameters, etc
datadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_thetaseq\';
filename = 'figure5_data_230723.mat';
load([datadir, filename])
postR_group = grpstats(postRPCR, {'group', 'Ctrl_speed_subset'}, {'mean', 'sem'});
behavioredges = -81:3:99;
params.dec_edges = -81:9:99; %could also try 9

gnames = {'random', 'gamma'};
dnames = {'pre', 'post'};
snames = {'speed1', 'speed2'};
groupcolors = params.colors; 

%%% custom colors for these plots
GS1 = hex2rgb('80aec5'); %gamma speed subset 1
GS2 = hex2rgb('196598'); %gamma speed subset 2 
RS1 = hex2rgb('a2cb7d'); %random speed subset 1
RS2 = hex2rgb('2a8924'); %random speed subset 2

c.gamma.speed1 = GS1;
c.gamma.speed2 = GS2;
c.random.speed1 = RS1;
c.random.speed2 = RS2; 

timeEdges = -0.2:0.02:0.2; %copied from params section of cf_decoding_thetaseq. ideally will resave structure with params files
posEdges = -81:3:99; %same as above
adjustedEdges = -180/2:3:180/2; %180 for now

%% plot 
iPlot = 1;
exT = [476 460];

figure(fh)
hold on

%%% slow 
axes(ax{iPlot})
hold on
tmpSeq = tmpData.PostR_Var5{exT(iPlot)};
imagesc(timeEdges(1:end-1)+0.01, adjustedEdges(1:end-1), tmpSeq, [0.01 0.035])
plot(timeEdges(1:end-1), zeros(1,size(tmpSeq,2)), 'w--')
plot(zeros(1,size(tmpSeq,1)), adjustedEdges(1:end-1), 'w--')
xlim([-0.08 0.08])
ylim([-45 45])
yticks([-45 0 45])
xticks([-0.08 0 0.08])
xlabel('Time from theta trough (s)')
ylabel('Relative decoded position (deg)')
title ('SLOW')
c = colorbar; 
c.Ticks = [0.01 0.035];
c.Label.String = 'Probability'; 

iPlot = iPlot+1;
%%% fast
axes(ax{iPlot})
hold on
tmpSeq = tmpData.PostR_Var5{exT(iPlot)};
imagesc(timeEdges(1:end-1)+0.01, adjustedEdges(1:end-1), tmpSeq, [0.01 0.045])
plot(timeEdges(1:end-1), zeros(1,size(tmpSeq,2)), 'w--')
plot(zeros(1,size(tmpSeq,1)), adjustedEdges(1:end-1), 'w--')
xlim([-0.08 0.08])
ylim([-45 45])
yticks([-45 0 45])
xticks([-0.08 0 0.08])
xlabel('Time from theta trough (s)')
ylabel('Relative decoded position (deg)')
title ('Fast')
c = colorbar; 
c.Ticks = [0.01 0.045];
c.Label.String = 'Probability'; 


end

