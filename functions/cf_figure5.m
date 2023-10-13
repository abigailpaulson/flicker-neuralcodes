function cf_figure5(dirs, params, allindex, metadata)
%Figure 5
%ALP 6/29/2023

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Figure5\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figure5';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure 5 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);
hold on

%% Panel A is a decoding example, I don't have that yet
ax_a1 = makesubplotwithletter(fh,3,4, 1, 'A', 'spanW', 1.3);
ax_a2 = makesubplotwithletter(fh,3,4, 3, 'A', 'spanW', 1.3);

cf_plot_prospectivecoding_byspeed_example(dirs, params, allindex, metadata, fh, {ax_a1, ax_a2}, statsFID, {'A', 'A'}, tablefilename)

%% Panel B, C, D are from the theta seq decoding figure
ax_b1 = makesubplotwithletter(fh,3,4, 5, 'B');
ax_b2 = makesubplotwithletter(fh,3,4, 6, 'B');
ax_c = makesubplotwithletter(fh,3,4,7, 'C', 'spanW', 2);
ax_d1 = makesubplotwithletter(fh,3,4, 9, 'D');
ax_d2 = makesubplotwithletter(fh,3,4, 10, 'D');
ax_e = makesubplotwithletter(fh,3,4, 11, 'E', 'spanW', 2);


cf_plot_prospectivecoding_byspeed(dirs, params, allindex, metadata, fh, {ax_b1, ax_b2, ax_c, ax_d1, ax_d2, ax_e}, statsFID, {'C1', 'C2', 'E1', 'E2'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figure5_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')


end

