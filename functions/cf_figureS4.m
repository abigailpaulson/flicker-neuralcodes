function cf_figureS4(dirs, params, allindex, metadata)
%cf_figureS4
%   supplementary figure 4
%ALP 3/21/2023

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Supplement4\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figureS4';
statsfilename = [figname, '_captionInfo_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Supplement 4 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);

%% Panel A, B, C are ripple decoding supplementary figures
ax_a =makesubplotwithletter(fh,3,3, 1, 'A');
ax_b = makesubplotwithletter(fh,3,3, 2, 'B');
ax_c = makesubplotwithletter(fh,3,3, 3, 'C');

cf_plot_rippledecoding_allripples(dirs, params, allindex, metadata, fh, {ax_a, ax_b, ax_c}, statsFID, {'A', 'B', 'C'}, tablefilename)

%% Panel D, E, and F are sharpwave ripple properties

ax_d = makesubplotwithletter(fh,3,3, 4, 'D');
ax_e = makesubplotwithletter(fh,3,3, 5, 'E');
ax_f = makesubplotwithletter(fh, 3,3, 6, 'F');

%RZ ripple properties 9/28/23
cf_plot_rippleproperties_RZ(dirs, params, allindex, metadata, fh, {ax_d, ax_e, ax_f}, statsFID, {'D', 'E', 'F'}, tablefilename)


%% Panel G, H, I, J are activation / coactivation properties
ax_g = makesubplotwithletter(fh,3,3, 7, 'G');
ax_h = makesubplotwithletter(fh,3,3, 8, 'H');
ax_i = makesubplotwithletter(fh,3,3, 9, 'I');
% ax_f = makesubplotwithletter(fh,4,5, 9, 'F');
% ax_g = makesubplotwithletter(fh,4,5, 10, 'G');

cf_plot_rippleactivation(dirs, params, allindex, metadata, fh, {ax_g, ax_h, ax_i}, statsFID, {'G', 'H', 'I'}, tablefilename)

%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save the figure as a pdf
filename = ['figureS4_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')






end

