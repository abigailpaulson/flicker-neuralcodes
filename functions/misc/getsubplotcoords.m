function [xx, yy, panelheight, panelwidth, midbuffer, axx, ayy] = getsubplotcoords(nRows, nCols)
%getsubplotcoords
%
%ALP 12/6/2022

sidebuffer = [0.07 0.05 0.02 0.03]; %L side, bottom, R side, top
midbuffer = [0.1 0.08]; %buffer between columns, buffer between rows

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
ayy = yy+panelheight-0.07;

end

