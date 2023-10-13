function cf_figure3(dirs, params, allindex, metadata)
%chronic flicker prospective coding manuscript
%   figure 3
%ALP 1/18/2023

%% helpful stuff
dayindex = unique(allindex(:,1:2), 'rows');
figdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\manuscriptfigures\Figure3\'; 

%% create stats file, rewrite it if it already exists
%this is so that everytime the file is run, the stats are regenerated
figname = 'figure3';
statsfilename = [figname, '_stats_', num2str(datestr(now, 'yymmdd')), '.txt'];
statsFID = fopen([figdir, statsfilename], 'w');
fprintf(statsFID, '%s\r\n', ['Figure 3 stats -- ' num2str(datestr(now))]);

tablefilename = [figname, '_stats_'];
tablefilename = [figdir, tablefilename];

%% initialize figure
fh = figure('units','inch','position',[0,0,6.5,6]);

%% Panel A is a diagram
%done in illustrator 
% ax_a = makesubplotwithletter(fh,3,3, 1, 'A', 'nosubplot');

%% new Panel A is an example of spiking and theta
%ALP 4/18/2023
%   I think by the time i make the desired edits to panel B below there
%   will be room for the example. I'm not sure though.... 
ax_a = makesubplotwithletter(fh,4,3, 1, 'A');
cf_plot_thetaspikingexample(dirs, params, allindex, metadata, fh, {ax_a})

%% Panel B and C is an example session
ax_b = makesubplotwithletter(fh,4,3, 2, 'B');
ax_c = makesubplotwithletter(fh,4,3, 3, 'C'); %this is going to be the same as the example in b1 but with boxes over the edges to show the forward and behind info

cf_plot_thetadecodingexample(dirs, params, allindex, metadata, fh, {ax_b, ax_c})

%% Panel D and E - group theta decoding over position
%%% panel E - quantify prospective coding approaching and in 1/2 of the reward zone

ax_d = makesubplotwithletter(fh,4,4,5,'D', 'spanW', 1.6); %line plot through space
ax_e = makesubplotwithletter(fh,4,4,7,'E');

cf_plot_thetadecoding_results(dirs, params, allindex, metadata, fh, {ax_d, ax_e}, statsFID, {'D', 'E'}, tablefilename)


%% add in speed control 
ax_f = makesubplotwithletter(fh, 4,4,9, 'F', 'spanW', 1.2);

cf_plot_RRZprospectivecoding_byspeed(dirs, params, allindex, metadata, fh, {ax_f}, statsFID, {'F1', 'F2'}, tablefilename)


%% panel E - quantify goal representation, place cells etc
%want to show similar amounts of goal represenation between groups
ax_h = makesubplotwithletter(fh, 4,4, 11, 'H'); 

cf_plot_placecells_goalrepresentation(dirs, params, allindex, metadata, fh, {ax_h}, statsFID, {'H'}, tablefilename)



%% finish stuff up
%close stats file
fclose(statsFID);

%make figure pretty 
makefigurepretty(fh,1)

%% save
filename = ['figure3_', num2str(datestr(now,'yymmdd'))];
savefigALP(figdir, filename, 'filetype', 'pdf')



end

