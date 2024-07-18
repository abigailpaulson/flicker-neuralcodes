%cf_AllFigures
%   this script runs all figures for Paulson et al. 2023
%ALP 12/13/2022

clear
%% load metadata
[params, dirs, metadata, allindex] = projectInfo2('chronicflicker_annulartrack', 'condition', 'prepostVR');

%% set run flags
run.figure1 = 0;
run.figure2 = 0;
run.figure3 = 0;
run.figure4 = 0;    
run.figure5 = 0;

run.figureS1 = 0;
run.figureS2 = 0;
run.figureS3 = 0;
run.figureS4 = 0;
run.figureS5 = 0;
run.figureS6 = 0;
run.figureS7 = 0;

run.revieweronly1 = 1;

%% figure 1
if run.figure1
    cf_figure1(dirs, params, allindex, metadata)
end

if run.figureS1
    cf_figureS1(dirs, params, allindex, metadata)
end

if run.figure3
    cf_figure3(dirs, params, allindex, metadata)
end

if run.figureS3 %now supp 5
    cf_figureS3(dirs, params, allindex, metadata)
end

if run.figure4
    cf_figure4(dirs, params, allindex,metadata)
end

if run.figureS4 % now supp 6
    cf_figureS4(dirs, params, allindex, metadata)
end
    
if run.figure5
    cf_figure5(dirs, params, allindex, metadata)
end

if run.figureS5
    cf_figureS5(dirs, params, allindex, metadata)
end

if run.figureS6 % now supplementary figure 7
    cf_figureS6(dirs, params, allindex, metadata)
end

if run.figureS7
    cf_figureS7(dirs, params, allindex, metadata)
end

if run.revieweronly1
    cf_revieweronly1(dirs, params, allindex, metadata)
end