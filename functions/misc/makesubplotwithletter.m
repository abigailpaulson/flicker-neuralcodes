function ax = makesubplotwithletter(fh, nRows, nCols, nPanel, letter, varargin)
%makesubplotwithletter
% modified from getsubplotcoords
%
%ALP 12/14/2022

sw = 1; sh = 1; sp = 1;
%defaults
lbuffer = 0.07;
bbuffer = 0.05;
rbuffer = 0.02;
tbuffer = 0.03;
colbuffer = 0.08;
rowbuffer = 0.08; %changed from 0.08 to make room for titles
for v = 1:2:length(varargin)
    if strcmp(varargin{v}, 'spanW')
        sw = varargin{v+1};
    elseif strcmp(varargin{v}, 'spanH')
        sh = varargin{v+1};
    elseif strcmp(varargin{v}, 'nosubplot')
        sp = 0;
    elseif strcmp(varargin{v}, 'Lbuffer')
        lbuffer = varargin{v+1};
    elseif strcmp(varargin{v}, 'Bbuffer')
        bbuffer = varargin{v+1};
    elseif strcmp(varargin{v}, 'Rbuffer')
        rbuffer = varargin{v+1};
    elseif strcmp(varargin{v}, 'Tbuffer')
        tbuffer = varargin{v+1};
    elseif strcmp(varargin{v}, 'colbuffer')
        colbuffer = varargin{v+1};
    elseif strcmp(varargin{v}, 'rowbuffer')
        rowbuffer = varargin{v+1};
    end
end

sp = logical(sp);

% sidebuffer = [0.07 0.05 0.02 0.03]; %L side, bottom, R side, top
% midbuffer = [0.08 0.08]; %buffer between columns, buffer between rows
sidebuffer = [lbuffer bbuffer rbuffer tbuffer];
midbuffer = [colbuffer rowbuffer];

panelheight = ((1-sidebuffer(2)-sidebuffer(4))-(midbuffer(2)*(nRows-1)))/nRows;
panelwidth = (1-sidebuffer(1)-sidebuffer(3)-(midbuffer(1)*(nCols-1)))/nCols;

for i = 1:nRows
    for j= 1:nCols
        xx(i,j) = sidebuffer(1)+midbuffer(1)*(j-1)+panelwidth*(j-1);
        yy(i,j) = sidebuffer(2)+midbuffer(2)*(i-1)+panelheight*(i-1);
    end
end
yy = flipud(yy);
axx = xx-0.07;
ayy = yy+(panelheight*sh+midbuffer(2)*(sh-1))-0.07;

%%% figure out what panel we're on
inds = 1:nRows*nCols;
indMat = reshape(inds, [nCols, nRows])'; %have to invert bc reshape works columnwise

[iPanelX, iPanelY] = find(indMat == nPanel);

%%% annotate
annotation(fh, 'textbox', [axx(iPanelX, iPanelY), ayy(iPanelX, iPanelY), 0.1, 0.1], 'String', letter, 'LineStyle', 'none', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', 'Margin', 1)

if sp
    ax = subplot('position', [xx(iPanelX, iPanelY) yy(iPanelX, iPanelY), panelwidth*sw+midbuffer(1)*(sw-1), panelheight*sh+midbuffer(2)*(sh-1)]);
else
    ax = [];
end





end

