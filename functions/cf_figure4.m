function cf_figure4(dirs, params, allindex, metadata)
%Figure 4
%ALP 3/23/2023

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Figure4\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figure4';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure 4 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);
hold on

%% Panel A, B, C are ripple decoding supplementary figures
ax_a = makesubplotwithletter(fh,3,3, 1, 'A');
ax_b = makesubplotwithletter(fh,3,3, 2, 'B');
ax_c = makesubplotwithletter(fh,3,3, 4, 'C');
ax_d = makesubplotwithletter(fh,3,3, 5, 'D');
ax_e = makesubplotwithletter(fh,3,3, 6, 'E');

%%% panel A make ripple example


cf_plot_rippleExample(dirs, params, allindex, metadata, fh, {ax_a}, statsFID, {'A'}, tablefilename)


cf_plot_rippledecoding_RZripples(dirs, params, allindex, metadata, fh, {ax_c, ax_d, ax_e}, statsFID, {'C', 'D', 'E'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figure4_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')



end

