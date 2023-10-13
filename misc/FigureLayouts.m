
sidebuffer = [0.1 0.05 0.05 0.03]; %L side, bottom, R side, top
midbuffer = [0.08 0.08]; %buffer between columns, buffer between rows

nRows = 3;
nCols = 3;

panelheight = ((1-sidebuffer(2)-sidebuffer(4))-(midbuffer(2)*(nRows-1)))/nRows; 
panelwidth = (1-sidebuffer(1)-sidebuffer(3)-(midbuffer(1)*(nCols-1)))/nCols;

for i = 1:nRows
    for j= 1:nCols
        xcoord(i,j) = sidebuffer(1)+midbuffer(1)*(j-1)+panelwidth*(j-1);
        ycoord(i,j) = sidebuffer(2)+midbuffer(2)*(i-1)+panelheight*(i-1);
    end
end

ycoord = flipud(ycoord);

a_xcoord = xcoord-0.07; 
a_ycoord = ycoord+panelheight-0.07;
pletters = 'ABCDEFGHIJKLMNOPQRS'; 


%3x3 figure layout

figure('units','inch','position',[0.1,0,6.5,9]);
hold on
count = 1;
for i = 1:nRows
    for j = 1:nCols
       annotation(gcf, 'textbox', [a_xcoord(i,j), a_ycoord(i,j), 0.1, 0.1], 'String', pletters(count), 'LineStyle', 'none', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold')
       subplot('position', [xcoord(i,j) ycoord(i,j) panelwidth, panelheight])
       hold on
       ylabel('text here')
       xlabel('text here')
       count = count+1;
    end
end
makefigurepretty(gcf)







figure('units','inch','position',[0.1,0,6.5,9]);
annotation(gcf, 'textbox', [0.04 0.24, 0.1, 0.1], 'String', 'G', 'LineStyle', 'none', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold')
subplot('position', [0.1 0.05 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.41 0.05 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.72 0.05 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')

hold on
subplot('position', [0.1 0.383 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.41 0.383 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.72 0.383 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')

hold on
subplot('position', [0.1 0.716 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.41 0.716 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
subplot('position', [0.72 0.716 0.2333 0.25])
hold on
ylabel('text here')
xlabel('text here')
makefigurepretty(gcf)

figure('units','inch','position',[0,0,6.5,9]);
hold on
subplot('position', [0.1 0.05 0.33 0.33])

