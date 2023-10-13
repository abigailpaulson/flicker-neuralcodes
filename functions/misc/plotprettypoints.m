function plotprettypoints(g, b, data, varargin)
% plotprettypoints
%
%   g - current figure handle
%   b - bar handle OR vector of x data
%   data - data to plot
%   sz - size of data points, 25 if empty
%
% ALP 6/13/2022

%%% set defaults
sz = 25;
darkpts = 0;
scwidth = 0.8;

%%% variable inputs
for v = 1:2:length(varargin)
    if strcmp(varargin{v}, 'size')
        sz = varargin{v+1};
    elseif strcmp(varargin{v}, 'color')
        CData = varargin{v+1}; 
    elseif strcmp(varargin{v}, 'darken')
        darkpts = 1; 
    elseif strcmp(varargin{v}, 'scatterwidth')
        scwidth = varargin{v+1};
    end
end

%%% check if bar input
if isobject(b)
    width = 0.9*b.BarWidth;
    xdata = b.XData; 
    CData = b.CData;
else
    width = 0.9*scwidth; 
    xdata = b; 
end

for i = 1:length(xdata)
    %%% darken color
    c = CData(i,:);
    if darkpts
        c = c.*0.7;
    end
    
    %%% create jitter vector
    r = rand([1,length(data{i})]);
    j = r*(width)-(width/2); 
    j = j+xdata(i);
    
    yvals = data{i};
    xvals = j; 
    
    figure(g)
    hold on
    scatter(xvals, yvals, sz, c, 'filled')
end


end

