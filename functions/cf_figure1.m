function cf_figure1(dirs, params, allindex, metadata)
%cf_figure1
%
%ALP 12/13/2022

%% helpful formatting stuff


%% subplot coordinates
f1 = figure('units','inch','position',[0,0,6.5,6]);
ax_a = makesubplotwithletter(f1,4,3, 1, 'A', 'nosubplot');
ax_b = makesubplotwithletter(f1,4,3, 2, 'B', 'nosubplot');
ax_c = makesubplotwithletter(f1,4,3, 3, 'C', 'nosubplot');

%% behavior example
ax_d = makesubplotwithletter(f1,8,3, 7, 'D', 'rowbuffer', 0.05, 'spanW', 0.85);
ax_d2 = makesubplotwithletter(f1,8,3, 10, 'D', 'rowbuffer', 0.05, 'spanW', 0.85);
ax_e = makesubplotwithletter(f1,4,3, 8, 'E', 'spanW', 2, 'spanH', 2);

cf_getbehaviorexample(dirs, params, allindex, metadata, f1, {ax_d, ax_d2, ax_e})

%% behavior, by trials
% ax_f = makesubplotwithletter(f1,4,4,9,'F');
% ax_g = makesubplotwithletter(f1,4,4,10,'G');
% ax_h = makesubplotwithletter(f1,4,4,11,'H');
% ax_i = makesubplotwithletter(f1,4,4,12,'I');
% 
% cf_plot_behaviorperformance_groups(dirs, params, allindex, metadata, f1, {ax_f, ax_g, ax_h, ax_i})

%% final touches
makefigurepretty(f1)

%% save PDF
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\'; 
filename = ['figure1_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')

end

