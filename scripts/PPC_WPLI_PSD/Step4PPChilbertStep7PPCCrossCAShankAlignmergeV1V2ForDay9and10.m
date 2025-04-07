clear all
LoadData='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step1\';
load([LoadData 'StimInfoList.mat'],'StimFileList');

SavePathLoad='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step4\';


SavePathRaw=['\\ad.gatech.edu\bme\labs\singer\Abby\code\flicker-neuralcodes\results\PPC\Manuscript\'];
if ~exist(SavePathRaw, 'dir')
mkdir(SavePathRaw)
end

SavePathRaw=[SavePathRaw '6bandCACrossShankColMergeV1V2Day910Notch40\'];
if ~exist(SavePathRaw, 'dir')
mkdir(SavePathRaw)
end

foldernm = '6bandCACrossShankColMergeV1V2Day910Notch40'; 

%%%%%%%%%%Integrate VR results

%%%%%%%%%%%%%Day 0 VR, and Day 1 VR before Flicker
% load([SavePathLoad '13bandVRhilbertPPCresultV2.mat']);

load([SavePathLoad '6bandVRhilbertPPCresultsCrossCAVRShankColPyrChAveSh12Notch40.mat']);
NeuroTemp1=NeuronData;
PPCFiletemp1=PPCvrFile;
TSppcVRNumTemp1=TSppcVRNumFile;
NonNanVRNumTemp1=NonNanVRNumFile;
%%%%%%%%%%%%%Day 0 VR, and Day 1 VR before Flicker


%%%%%%%%%%%%%Day 9 and 10, using all VR
load([SavePathLoad '6bandVRhilbertPPCresultsCrossCAVRShankColPyrChMergeV1V2AveSh12Notch40.mat']);
NeuroTemp2=NeuronData;
PPCFiletemp2=PPCvrFile;
TSppcVRNumTemp2=TSppcVRNumFile;
NonNanVRNumTemp2=NonNanVRNumFile;
%%%%%%%%%%%%%Day 9 and 10, using all VR

%% run after stuff is loaded
savefiguredir = ['']; 


NonNanVRNum=[];
TSVRNum=[];
PPCvr=[];
CellAllProp.CellType=[];
CellAllProp.SubjFile=[];
CellAllProp.brainReg=[];
CellAllProp.TSIDtotal=[];
CellAllProp.TSIDtotal=[];
CellAllProp.Subj=[];

FrePlot=mean(passband);
% FreName={'Delta' 'Theta' 'Beta','G1','G2','G3-40Hz','G4','G5','G6','G7','G8','G9','G10'};
% FreName={'Delta' 'Theta' 'Beta','G1','G2','G3','G4'};
FreName={'Delta' 'Theta' 'Beta','Sgamma','Mgamma','Fgamma'};

FreName2={};
for iFF=1:length(FreName)
    FreName2{iFF}=[num2str(passband(1,iFF)) '-' num2str(passband(2,iFF)) ' Hz'] ;
end

% FI=[2 3 4 5 6 7 8];
FI=1:length(FreName);

for i=1:length(NeuronData)
%            if i==40
%               continue;
%            end


if StimFileList(i).StimDay(1)==0||StimFileList(i).StimDay(1)==1
    if i==40
        Iv=isnan(NeuroTemp1(i).CellType);
        if ~isempty(Iv)
            i
            NeuroTemp1(i).CellType(Iv)=[];
            NeuroTemp1(i).brainReg(Iv)=[];
            NeuroTemp1(i).TSIDtotal(Iv)=[];
            
            PPCFiletemp1{i}(Iv,:,:,:)=[];
            TSppcVRNumTemp1{i}(Iv,:,:)=[];
            NonNanVRNumTemp1{i}(Iv,:,:)=[];
        end
    end
    
    if isempty(PPCFiletemp1{i})
        continue
    end
    
    PPCvr=cat(1,PPCvr,squeeze(PPCFiletemp1{i}(:,:,1,:)));
    NonNanVRNum=cat(1,NonNanVRNum,squeeze(NonNanVRNumTemp1{i}(:,1,:)));
    TSVRNum=cat(1,TSVRNum,squeeze(TSppcVRNumTemp1{i}(:,1,:)));
    
    
elseif StimFileList(i).StimDay(1)==9||StimFileList(i).StimDay(1)==10
    if i==40
        Iv=isnan(NeuroTemp2(i).CellType);
        if ~isempty(Iv)
            i
            NeuroTemp2(i).CellType(Iv)=[];
            NeuroTemp2(i).brainReg(Iv)=[];
            NeuroTemp2(i).TSIDtotal(Iv)=[];
            PPCFiletemp2{i}(Iv,:,:)=[];
            TSppcVRNumTemp2{i}(Iv,:)=[];
            NonNanVRNumTemp2{i}(Iv,:)=[];
        end
    end
    if isempty(PPCFiletemp2{i})
        continue
    end
    
    
    PPCvr=cat(1,PPCvr,PPCFiletemp2{i});
    NonNanVRNum=cat(1,NonNanVRNum,NonNanVRNumTemp2{i});
    TSVRNum=cat(1,TSVRNum,TSppcVRNumTemp2{i});
    
    
    
else
    
end



CellAllProp.CellType=[CellAllProp.CellType;NeuroTemp2(i).CellType];
CellAllProp.brainReg=[CellAllProp.brainReg;NeuroTemp2(i).brainReg];
CellAllProp.TSIDtotal=[CellAllProp.TSIDtotal;NeuroTemp2(i).TSIDtotal];
CellAllProp.SubjFile=[CellAllProp.SubjFile;repmat(i,length(NeuroTemp2(i).TSIDtotal),1)];
CellAllProp.Subj=[CellAllProp.Subj;repmat(StimFileList(i).Subj,length(NeuroTemp2(i).brainReg),1)];

end
%%%%%%%%%%Integrate VR results


% StimFileList




StimDay=[];
Subj=[];
FlickerT=[];
for i=1:length(StimFileList)
    Subj(i)=StimFileList(i).Subj;
    StimDay(i)=StimFileList(i).StimDay(1);
    if isempty(setdiff(StimFileList(i).FlickerType,0))
        FlickerT(i)=0;
    else
        FlickerT(i)=setdiff(StimFileList(i).FlickerType,0);
    end
end
StimDayAll=StimDay;

Temp1=[];
Temp2=[];
Temp3=[];
for i=1:length(StimFileList)
    n=sum(CellAllProp.SubjFile==i);
    Temp1=[Temp1;repmat(StimDay(i),n,1)];
    Temp2=[Temp2;repmat(FlickerT(i),n,1)];
    Temp3=[Temp3;repmat(StimFileList(i).Group,n,1)];
end
%
Temp4=Temp1;
Temp4(Temp1==10)=4;
Temp4(Temp1==9)=3;
Temp4(Temp1==1)=2;
Temp4(Temp1==0)=1;
CellAllProp.StimDay=Temp1;
CellAllProp.Flicker=Temp2;
CellAllProp.FlickerG=Temp3;
CellAllProp.StimDSeq=Temp4;
CellAllProp.CellID=[1:length(CellAllProp.CellType)]';








DaySeqGroup=[1 2;3 4];
DayGroup=[0 1;9 10];
%
DayName={'Pre' 'Post'};
StimGroup=[1 4];
StimName={'random','40Hz'};
CellGroup=[1 2];
CellName={'Pyr','Int'};
RegionGroup=[1 2];
RegionName={'CA1','CA3'};
%



%%%%%%%%%%%Set up threshold to include cells for analysis
SpikeThAll=[50 50 100 100 200 200];
TrialThAll=[5 10 5 10 5 10];

%%%%%%%%%%Set up threshold to include cells for analysis, for manuscript
SpikeThAll=[100];
TrialThAll=[5];
% FI=2;

YlimCellTheta(1,:)=[-0.001 0.008];
YlimCellTheta(2,:)=[-0.001 0.01];

YtickCellTheta(1,:)=[0:0.004:0.008];
YtickCellTheta(2,:)=[0:0.015:0.03];

YlimCellTheta(1,:)=[-0.001 0.008];
YlimCellTheta(2,:)=[-0.001 0.01];


%%%%%%%%%%%Violin Plot, rank test and t test
for jTh=1:length(SpikeThAll)
    
    SpikeTh=SpikeThAll(jTh);
    TrialTh=TrialThAll(jTh);
    SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
    if ~exist(SavePathRaw0, 'dir')
    mkdir(SavePathRaw0);
    end
    %%%%%%%%%%%Set up threshold to include cells for analysis
    
    nBoot=10000;
    
    % YlimCell=[-0.00001 0.0003;-0.001 0.008;-0.00001 0.0003;-0.00001 0.0003;];
    % YlimCellBand=[-0.00001 0.0005;-0.001 0.008;repmat([-0.00001 0.0003],3,1);repmat([-0.00001 0.00005],3,1);repmat([-0.00001 0.00002],5,1);];
    YlimCellBand=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);];  %%CA1
    
    
    YlimCellRegion{1}=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);];  %%CA1
    YlimCellRegion{1}(4:end,:)=YlimCellRegion{1}(4:end,:)/2;
    YlimCellRegion{2}=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);];  %%CA3
    
    
    
    YtickCellBand=YlimCellBand;
    
    % YlimCellTheta(2,:)=[-0.001 0.01];
    %
    % YtickCellTheta(1,:)=[0:0.004:0.008];
    % YtickCellTheta(2,:)=[0:0.015:0.03];
    %
    %
    %
    % YlimCellGamma(1,:)=[-0.00001 0.0002];
    % YlimCellGamma(2,:)=[-0.00001 0.0005];
    % YtickCellGamma(1,:)=[-0.00001:0.00025:0.0005];
    % YtickCellGamma(2,:)=[-0.00001:0.001:0.002];
    
    
    PlotColor2{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0
    PlotColor2{1,2}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 1
    PlotColor2{1,3}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9
    PlotColor2{1,4}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 10
    PlotColor2{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
    PlotColor2{2,2}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
    PlotColor2{2,3}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
    PlotColor2{2,4}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
    
    PlotColor4=[0.698 0.8745 0.5412;0.6510 0.8078 0.8902;0.2 0.6275 0.1725;0.1216 0.4706 0.7059]; %%%%%%%40Hz Day 1
    
    
    % close all
    GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
    
    P.xLeft=0.03;        %%%%%%Left Margin
    P.xRight=0.02;       %%%%%%Right Margin
    P.yTop=0.08;         %%%%%%Top Margin
    P.yBottom=0.08;      %%%%%%Bottom Margin
    P.xInt=0.05;         %%%%%%Width-interval between subplots
    P.yInt=0.08;         %%%%%%Height-interval between subplots
    
    
    
    clear UnPairedDayFlicker
    for iGroup=1 %iGroup = 3 is correct trial, 1 is all trial %using all trials ALP 5/2/24

        if iGroup==2
            continue;
        end
        
        SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
        if ~exist(SavePath, 'dir')
        mkdir(SavePath);
        end
        
        
        statsfilename = ['PPC_', foldernm, '_stats_', num2str(datestr(now, 'yymmdd')), '.txt'];
        statsFID = fopen([SavePath, statsfilename], 'w');
        fprintf(statsFID, '%s\r\n', ['PPC stats -- ' num2str(datestr(now))]);
        
        
        % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
        % % % length(ValidPair)
        % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
        % % % length(ValidPair)
        % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
        % % % length(ValidPair)
        
        
        ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
        length(ValidD0)
        ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
        length(ValidD0)
        ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
        length(ValidD0)
        
        
        GroupPair.CorrName='fdr';
        GroupPair.Test='Ranktest';
        GroupPair.Q=0.05;
        GroupPair.Pair=[1 1 2 3;2 3 4 4];
        
        %    GroupPair.Pair=[1 2;3 4];
        
        GroupPair.SignY=0.01;
        GroupPair.Plot=1;
        GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
        GroupPair.SamplePlot=0; %%%%%%%%%Plot Individual Sample Point
        GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
        GroupPair.LimY=[0 GroupPair.SignY*1.2];
        GroupPair.Marker={'o'};
        GroupPair.ViolinLR=[0 1 0 1];
        
        % PlotColor2=[0.6 0.6 0.6;0.1 0.1 0.1];
        PlotColor{1}=[0.6 0.6 0.6;0.8 0.1 0.2];
        PlotColor{2}=[0.6 0.6 0.6;0.2 0.8 0.1];
        
        
        % % Param3=Param2;
        % % Param3.TimeRepeatAnova=1;
        % % Param3.GroupRepeatAnova=0;
        % % Param3.RepeatAnova=1;
        % % Param3.Paired=0;
        % %
        % %
        VRname={'PreFlicker' 'PostFlicker'};
        SaveTemp=[SavePath 'VR-40vsRandom\'];
        if ~exist(SaveTemp, 'dir')
            mkdir(SaveTemp)
        end
        
        for iCell=1:length(CellGroup)
            dataforstats = []; %ALP 9/15/23
            
            for iFF=1:length(FI)
                jFF=FI(iFF);
                SaveTempF=[SaveTemp];
                if ~exist(SaveTempF, 'dir')
                mkdir(SaveTempF);
                end
                
                % %     YlimCell=repmat(YlimCellBand(jFF,:),2,1);
                % %     YtickCell=repmat(YtickCellBand(jFF,:),2,1);
                % %     if jFF<=2
                % %         YlimCell=YlimCellTheta;
                % %         YtickCell=YtickCellTheta;
                % %     else
                % %         YlimCell=YlimCellGamma;
                % %         YtickCell=YtickCellGamma;
                % %
                % %     end
                
                for iReg=1:length(RegionGroup)
                    YlimCell=repmat(YlimCellRegion{iReg}(jFF,:),2,1);
                    %YtickCell=YlimCell;
                    YtickCell = [YlimCell(:,1) zeros(length(YlimCell(:,1)), 1) YlimCell(:,2)];
                    
                    %                     subplotLU(length(RegionGroup),length(FI),iReg,iFF,P)
                    GroupPair.SignY=YlimCell(iCell,2);
                    GroupPair.LimY=[0 YlimCell(iCell,2)];
                    clear DataPlot xSubjID DataPlotName;
                    DataPlot={}; xSubjID={}; DataPlotName={}; ViolinPlotName = {}; xCellID = {}; xAnID = {};%ALP added 5/9/24
                    for iDay=1:size(DayGroup,1)
                        
                        for iS=1:length(StimGroup)
                            Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
                            Index2=find(CellAllProp.FlickerG==StimGroup(iS));
                            Index3=find(CellAllProp.CellType==CellGroup(iCell));
                            
                            
                            IndexNeed=intersect(intersect(Index1,Index2),Index3);
                            
                            Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
                            IndexNeed=intersect(Index4,IndexNeed);
                            
                            CellNeed=intersect(IndexNeed,ValidD0);
                            
                            if isempty(CellNeed)
                                continue;
                            end
                            DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
                            xSubjID{end+1}=CellAllProp.SubjFile(CellNeed); %this is the day index (1:64)
                            xAnID{end+1}=CellAllProp.Subj(CellNeed); %ALP 5/9/24 this is the animal ID
                            xCellID{end+1} = CellAllProp.CellID(CellNeed); % ALP 9/20/23
                            DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
                            ViolinPlotName{end+1} = [DayName{iDay}, '_', StimName{iS}];

                            PlotColor3(iS,:)=[PlotColor2{iS,iDay}(2,:)];
                            
                        end
                        % %
                    end
                    
                    %%% get data for stats ALP 9/15/23
                    for sg = 3:4 %stats group
                        tmp = []; tmpgroup = []; tmpreg = []; tmpfreq = [];
                        tmp = [DataPlot{sg} xAnID{sg} xCellID{sg}]; %edited from xSubjID to xAnID ALP 5/22/24
                        tmpgroup = sg.*ones(size(tmp,1),1);
                        tmpreg = iReg.*ones(size(tmp,1),1);
                        tmpfreq = iFF.*ones(size(tmp,1),1);
                        stat_append = [tmp tmpgroup tmpreg tmpfreq];
                        dataforstats = [dataforstats; stat_append];
                    end
                    %data for stats
                    %data animalID groupID cellID regionID frequencyID
                    %%%
                    
                    figure('Position', [273 608 172 191])
                    hold on
                    vdat = [];
                    % %     subplotLU(1,length(DayGroup),iCell,iDay);
                    GroupPair.Marker={'o','^','o','^'};
                    
                    PathSave=[SaveTempF RegionName{iReg} CellName{iCell} FreName{jFF} '.txt'];
                    Datatype=0;
                    GroupPair.SamplePairedPlot=0;
                    x=[1 1.4 3 3.4];
                    
                    %%% ALP plot
                    for iv = 3:4
                        vdat.(ViolinPlotName{iv}) = DataPlot{iv};
                    end
                    violinplot_half(vdat, [], 'ViolinColorMat', PlotColor4(3:4,:), 'BoxWidth', 0.018, 'MedianSize', 25, 'ShowData', false);
                    set(gca, 'xtick', [], 'ylim',YlimCell(iCell,:),'ytick',YtickCell(iCell,:),'xticklabel',[]);
                    ylabel([RegionName{iReg} ' PPC']);
                    title([CellName{iCell} ' ', ViolinPlotName{3}, ' ', num2str(length(DataPlot{3})), ' ', ViolinPlotName{4}, ' ', num2str(length(DataPlot{4}))]);
                    
                    txtentry = [RegionName{iReg}, ' ', CellName{iCell}, ' ' FreName2{jFF}];
                    tblname = [RegionName{iReg}, '_', CellName{iCell}, '_' FreName2{jFF}];
                    cf_stats2txt2(DataPlot(3:4), statsFID, txtentry, 'cells', 'PPC', 'PPC', StimName, [SaveTemp, tblname])
                    xlabel(FreName2{jFF});
                    
                    figname = [RegionName{iReg}, '_', CellName{iCell}, '_', FreName{jFF}];
                    savefigALP(SaveTemp, figname)
                    
                    
                    %%% per animal plots
                    vdat = [];
                    stimname = {'', '', 'random', 'gamma'};
                    figure('Units', 'inches', 'Position', [0 0 6.5 1.8])
                    hold on
                    cmat = [];
                    iPlot = 1;
                    for sg = 3:4
                        animals = unique(xAnID{sg});
                        nUnitsTemp = tabulate(xAnID{sg});
                        tmpAn = ismember(nUnitsTemp(:,1), animals);
                        nUnits = nUnitsTemp(tmpAn,2);
                        for an = 1:length(animals)
                            isDat = xAnID{sg} == animals(an);
                            vdat.(['an', num2str(animals(an))]) = DataPlot{sg}(isDat);
                            cmat{iPlot} = PlotColor4(sg,:);
                            iPlot = iPlot+1;
                        end
                       disp(['VR ' stimname{sg}, ' ', CellName{iCell} ' ', FreName{jFF}, ' ', RegionName{iReg}, ' cell range: ', num2str(min(nUnits)), ' - ', num2str(max(nUnits))])

                    end
                    v = violinplot(vdat, [], 'BoxWidth', 0.02, 'ShowData', false);
                    for vv = 1:length(v)
                        v(vv).ViolinColor = cmat{vv};
                        v(vv).BoxColor = cmat{vv};
                    end
                    ylabel([RegionName{iReg}, ' PPC'])
                    title([CellName{iCell} ' ', FreName{jFF}, ' ', RegionName{iReg}]);
                    figname = ['PerAnimal', RegionName{iReg}, '_', CellName{iCell}, '_', FreName{jFF}];
                    makefigurepretty(gcf, 1)
                    %savefigALP(SaveTemp, figname, 'ftype', 'pdf')
                    
                    %%%% end ALP
                    
                    %PathSave=[SaveTempF RegionName{iReg} CellName{iCell} FreName{jFF} '.txt'];
                    %UnPairedDayFlicker{iGroup,iCell,iFF}=ErrorViolinHalf(x,DataPlot,PlotColor4,Datatype,PathSave,GroupPair,[1 2 3 4]);
                    
                    %     VRUnPaired40VsRandom{iCell,iDay,iFF}=ErrorBoxPlotLU(1:2,DataPlot,PlotColor3,[],GroupPair,[1 2]);
                    %set(gca,'xlim',[0 4],'xtick',x,'ylim',YlimCell(iCell,:),'ytick',YtickCell(iCell,:),'xticklabel',[]);
                    %     if iFF==1
                    %         ylabel([RegionName{iReg} ' PPC']);
                    %     end
                    %%%%PosXY(1,:): x coordinate;PosXY(2,:): y coordinate; coordinate of plot
                    %%%%PosXY(3,:): x coordinate;PosXY(4,:): y coordinate; coordinate of labels
                    
                    %         if iFF==2
                    %
                    %     yTemp=YtickCell(iCell,2);
                    %     tickStep=mean(diff(YtickCell(iCell,:)));
                    %     yTemp=[yTemp yTemp-tickStep*0.05];
                    %     tempXY=[repmat(x(1),1,2);yTemp;repmat(x(1)+0.2,1,2);yTemp];
                    %     LuLegend(tempXY,0,DataPlotName(1:2),PlotColor4(1:2,:),6);
                    %
                    % %     yTemp=[yTemp yTemp-tickStep*0.1];
                    %     tempXY=[repmat(x(3),1,2);yTemp;repmat(x(3)+0.2,1,2);yTemp];
                    %     LuLegend(tempXY,0,DataPlotName(3:4),PlotColor4(3:4,:),6);
                    %         end
                    %                 if iReg==2
                    %                    xlabel(FreName2{jFF});
                    %                 end
                    % %     LuFontStandard;
                    %
                    %     end
                    
                    
                end
                %     papersizePX=[0 0 6*length(FI) 6*length(RegionGroup)];
                %     set(gcf, 'PaperUnits', 'centimeters');
                %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
                %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'pdf');
                %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'png');
                %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker.eps'],'epsc');
                %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'fig');
                %     close all
                
            end
            % % save([SaveTemp 'statis.mat'],'VRUnPaired40VsRandom');
            % % close all
            
            %data for stats
            %data animalID groupID regionID frequencyID
            PPCTable = array2table(dataforstats, 'VariableNames', {'ppc', 'anID', 'cellID', 'groupID', 'regID', 'freqID'});
            tablefilename = ['\\ad.gatech.edu\bme\labs\singer\Abby\code\flicker-neuralcodes\results\LMM_R\', 'TableData_PPC_VR_CrossReg_', CellName{iCell}, '_', GroupName{iGroup} '.txt']; %added groupname ALP 5/2/24 
            writetable(PPCTable, tablefilename)
            
            clear dataforstats
        end
        
        fclose(statsFID);
    end
end


%%% ALP 9/15/2023
%%% get data for LME in R
% want region = 2 %ca3 spiking to ca1 LFP
% want freq = [2 4 5]
% need to add cell ID
% make giant matrix and save



%%% 
    
    
    
    
    %%%%%%%%%%%%Mixed-Linear Model for above analysis
%     Rfolder='Y:\singer\LuZhang\Project2ChronicF\Analysis\RFolder\';
    
    
    %%%%%%%%Delete all previous LME results in .txt
    for jTh=1:length(SpikeThAll)
        % for jTh=3:3
        
        SpikeTh=SpikeThAll(jTh);
        TrialTh=TrialThAll(jTh);
        SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
        % mkdir(SavePathRaw0);
        %%%%%%%%%%%Set up threshold to include cells for analysis
        
        close all
        GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
        
        clear UnPairedDayFlicker
        for iGroup=1:3
            if iGroup==2
                continue;
            end
            %
            %     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
            %     mkdir(SavePath);
            
            
            
            
            % %
            VRname={'PreFlicker' 'PostFlicker'};
            SaveTemp=[SavePath 'VR-40vsRandom\'];
            % mkdir(SaveTemp)
            
            for iCell=1:length(CellGroup)
                
                %     delete([SaveTemp 'LME*.txt']);
                
            end
            
        end
        
    end
    %%%%%%%%Delete all previous LME results in .txt
    
    
    %%%%%%Animal as random effect, and plot each animal's data in R
    for jTh=1:length(SpikeThAll)
        % for jTh=3:3
        
        SpikeTh=SpikeThAll(jTh);
        TrialTh=TrialThAll(jTh);
        SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
        % mkdir(SavePathRaw0);
        %%%%%%%%%%%Set up threshold to include cells for analysis
        
        close all
        GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
        
        clear UnPairedDayFlicker
        for iGroup=1:3
            if iGroup==2
                continue;
            end
            
            SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
            %     mkdir(SavePath);
            
            
            
            % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
            % % % length(ValidPair)
            % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
            % % % length(ValidPair)
            % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
            % % % length(ValidPair)
            
            
            ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
            length(ValidD0)
            ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
            length(ValidD0)
            ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
            length(ValidD0)
            
            
            % %
            VRname={'PreFlicker' 'PostFlicker'};
            SaveTemp=[SavePath 'VR-40vsRandom\'];
            % mkdir(SaveTemp)
            
            for iCell=1:length(CellGroup)
                
                tempFile=[Rfolder 'tempLME.txt'];
                %     fh=fopen(tempFile,'a');
                
                for iFF=1:length(FI)
                    jFF=FI(iFF);
                    SaveTempF=[SaveTemp];
                    %     mkdir(SaveTempF);
                    % %     if jFF<=2
                    % %         YlimCell=YlimCellTheta;
                    % %         YtickCell=YtickCellTheta;
                    % %     else
                    % %         YlimCell=YlimCellGamma;
                    % %         YtickCell=YtickCellGamma;
                    % %
                    % %     end
                    
                    for iReg=1:length(RegionGroup)
                        clear DataPlot xSubjID  xSubj DataPlotName FlickerID DayID StimDayID CellID StimDSeq;
                        DataPlot={}; xSubjID={};xSubj={}; DataPlotName={}; FlickerID={}; DayID={}; StimDayID={};CellID={};StimDSeq={};
                        
                        for iDay=1:size(DayGroup,1)
                            
                            for iS=1:length(StimGroup)
                                Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
                                Index2=find(CellAllProp.FlickerG==StimGroup(iS));
                                Index3=find(CellAllProp.CellType==CellGroup(iCell));
                                
                                
                                IndexNeed=intersect(intersect(Index1,Index2),Index3);
                                
                                Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
                                IndexNeed=intersect(Index4,IndexNeed);
                                
                                CellNeed=intersect(IndexNeed,ValidD0);
                                
                                if isempty(CellNeed)
                                    continue;
                                end
                                DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
                                xSubjID{end+1}=CellAllProp.SubjFile(CellNeed);
                                xSubj{end+1}=CellAllProp.Subj(CellNeed);
                                FlickerID{end+1}=repmat(iS-1,length(CellNeed),1);
                                DayID{end+1}=repmat(iDay-1,length(CellNeed),1);
                                CellID{end+1}=CellAllProp.CellID(CellNeed);
                                StimDayID{end+1}=CellAllProp.StimDay(CellNeed);
                                StimDSeq{end+1}=CellAllProp.StimDSeq(CellNeed);
                                
                                DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
                                
                                
                            end
                            % %
                        end
                        
                        DataLME=[];
                        Cov1=[];
                        Cov2=[];
                        Cov3=[];
                        Cov4=[];
                        Cov5=[];
                        Cov6=[];
                        Cov7=[];
                        
                        for iData=1:length(DataPlot)
                            DataLME=[DataLME;DataPlot{iData}(:)];
                            Cov1=[Cov1;xSubjID{iData}(:)];
                            Cov2=[Cov2;xSubj{iData}(:)];
                            Cov3=[Cov3;FlickerID{iData}(:)];
                            Cov4=[Cov4;DayID{iData}(:)];
                            Cov5=[Cov5;StimDayID{iData}(:)];
                            Cov6=[Cov6;CellID{iData}(:)];
                            Cov7=[Cov7;StimDSeq{iData}(:)];
                            
                        end
                        tbl = array2table([DataLME(:),Cov1(:),Cov2(:),Cov3(:) Cov4(:) Cov5(:) Cov6(:) Cov7(:)],'VariableNames',{'PPC','Session','Animal','Flicker','Day','StimDay','Cell','StimDSeq'});
                        %   lme = fitlme(tbl,'Firing ~ Speed * Spa + (Speed|CellI)');
                        %   lme0 = fitlme(tbl,'Firing ~ Speed + (1|CellSubj)','FitMethod','ML');
                        %   lme1 = fitlme(tbl,'Firing ~ Speed + Spa + (1|CellSubj)','FitMethod','ML');
                        %   lme2 = fitlme(tbl,'Firing ~ Speed * Spa -Spa+ (1|CellSubj)','FitMethod','ML');
                        %   lme3 = fitlme(tbl,'Firing ~ Speed*Spa + (1|CellSubj)','FitMethod','ML');
                        
                        %   results = compare(lme0,lme1);   %%%%%%%speed affect vs speed+spatial effect
                        %   results = compare(lme1,lme2);   %%%%%%%speed affect vs speed*spatial effect
                        %   results = compare(lme2,lme3);   %%%%%%%speed affect vs speed*spatial effect
                        
                        %     writetable(tbl,[Rfolder 'tempData.csv'])
                        %     fprintf(fh,'\n%s\n',[RegionName{iReg} ' ' FreName2{jFF}]);
                        
                        %     RunRcode([Rfolder 'RplotAnimalDay.R'],'C:\Program Files\R\R-4.0.2\bin')
                        % % % %     RunRcode([Rfolder 'RplotAnimalDay.R'])
                        %     delete([Rfolder 'tempData.csv']);
                        
                        
                        ResultFigureFile=[SaveTemp RegionName{iReg} CellName{iCell} FreName{jFF} 'AnimalDay.png'];
                        %delete(ResultFigureFile);
                        tempFigure=[Rfolder 'PPCAnimalDay.png'];
                        %     movefile(tempFigure,ResultFigureFile);
                        
                    end
                    
                    
                end
                %     fclose(fh);
                ResultFile=[SaveTemp 'LME' CellName{iCell} 'PPCvsFlickerAndDay.txt'];
                % %     delete([SaveTemp 'LME*.txt']);
                %     movefile(tempFile,ResultFile);
                
            end
            % % save([SaveTemp 'statis.mat'],'VRUnPaired40VsRandom');
            % % close all
            
        end
        
    end
    %%%%%%Animal as random effect, and plot each animal's data in R
    
    
    %%%%%%Mixed ANOVA
    for jTh=1:length(SpikeThAll)
        % for jTh=3:3
        
        SpikeTh=SpikeThAll(jTh);
        TrialTh=TrialThAll(jTh);
        SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
        % mkdir(SavePathRaw0);
        %%%%%%%%%%%Set up threshold to include cells for analysis
        
        close all
        GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
        
        clear UnPairedDayFlicker
        for iGroup=1:3
            if iGroup==2
                continue;
            end
            
            SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
            %     mkdir(SavePath);
            
            
            
            % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
            % % % length(ValidPair)
            % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
            % % % length(ValidPair)
            % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
            % % % length(ValidPair)
            
            
            ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
            length(ValidD0)
            ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
            length(ValidD0)
            ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
            length(ValidD0)
            
            
            % %
            VRname={'PreFlicker' 'PostFlicker'};
            SaveTemp=[SavePath 'VR-40vsRandom\'];
            % mkdir(SaveTemp)
            SaveTemp=[SaveTemp 'MixedAnova\'];
            % mkdir(SaveTemp)
            
            for iCell=1:length(CellGroup)
                
                
                for iFF=1:length(FI)
                    jFF=FI(iFF);
                    SaveTempF=[SaveTemp];
                    %     mkdir(SaveTempF);
                    % %     if jFF<=2
                    % %         YlimCell=YlimCellTheta;
                    % %         YtickCell=YtickCellTheta;
                    % %     else
                    % %         YlimCell=YlimCellGamma;
                    % %         YtickCell=YtickCellGamma;
                    % %
                    % %     end
                    
                    for iReg=1:length(RegionGroup)
                        clear DataPlot xSubjID  xSubj DataPlotName FlickerID DayID StimDayID CellID StimDSeq;
                        DataPlot={}; xSubjID={};xSubj={}; DataPlotName={}; FlickerID={}; DayID={}; StimDayID={};CellID={};StimDSeq={};
                        
                        for iDay=1:size(DayGroup,1)
                            
                            for iS=1:length(StimGroup)
                                Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
                                Index2=find(CellAllProp.FlickerG==StimGroup(iS));
                                Index3=find(CellAllProp.CellType==CellGroup(iCell));
                                
                                
                                IndexNeed=intersect(intersect(Index1,Index2),Index3);
                                
                                Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
                                IndexNeed=intersect(Index4,IndexNeed);
                                
                                CellNeed=intersect(IndexNeed,ValidD0);
                                
                                if isempty(CellNeed)
                                    continue;
                                end
                                DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
                                xSubjID{end+1}=CellAllProp.SubjFile(CellNeed);
                                xSubj{end+1}=CellAllProp.Subj(CellNeed);
                                FlickerID{end+1}=repmat(iS-1,length(CellNeed),1);
                                DayID{end+1}=repmat(iDay-1,length(CellNeed),1);
                                CellID{end+1}=CellAllProp.CellID(CellNeed);
                                StimDayID{end+1}=CellAllProp.StimDay(CellNeed);
                                StimDSeq{end+1}=CellAllProp.StimDSeq(CellNeed);
                                
                                DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
                                
                                
                            end
                            % %
                        end
                        
                        DataLME=[];
                        Cov1=[];
                        Cov2=[];
                        Cov3=[];
                        Cov4=[];
                        Cov5=[];
                        Cov6=[];
                        Cov7=[];
                        
                        for iData=1:length(DataPlot)
                            DataLME=[DataLME;DataPlot{iData}(:)];
                            Cov1=[Cov1;xSubjID{iData}(:)];
                            Cov2=[Cov2;xSubj{iData}(:)];
                            Cov3=[Cov3;FlickerID{iData}(:)];
                            Cov4=[Cov4;DayID{iData}(:)];
                            Cov5=[Cov5;StimDayID{iData}(:)];
                            Cov6=[Cov6;CellID{iData}(:)];
                            Cov7=[Cov7;StimDSeq{iData}(:)];
                            
                        end
                        
                        %     MixAnova2F1W1B_withR(DataLME,[Cov4(:) Cov3(:) Cov2(:)],[SaveTempF RegionName{iReg} CellName{iCell} FreName2{jFF}]);
                    end
                    
                    
                end
                
            end
            % % save([SaveTemp 'statis.mat'],'VRUnPaired40VsRandom');
            % % close all
            
        end
        
    end
    %%%%%%Mixed ANOVA
    
    
    
    
    
    
    %%%%%%
    SpikeTh=100;
    TrialTh=5;
    iGroup=3;    %%%%%%%%Only Correct trial
    StimDay=StimDayAll;
    
    
    ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
    length(ValidD0)
    ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
    length(ValidD0)
    ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
    length(ValidD0)
    
    
    
    clear TableW;
    TableW{1,1}='Animals';   %%1
    TableW{1,end+1}='Session';
    TableW{1,end+1}='Flicker'; %%2
    TableW{1,end+1}='FlickerDay';
    TableW{1,end+1}='PrePostFlicker';
    
    TableW{1,end+1}='CellNum';
    TableW{1,end+1}='CA1PyrNum';
    TableW{1,end+1}='CA1IntNum';
    TableW{1,end+1}='CA3PyrNum';
    TableW{1,end+1}='CA3IntNum';
    
    TableW{1,end+1}='CA1PyrNumPPC';
    TableW{1,end+1}='CA1IntNumPPC';
    TableW{1,end+1}='CA3PyrNumPPC';
    TableW{1,end+1}='CA3IntNumPPC';
    
    
    
    
    for i=1:length(StimFileList)
        TableW{end+1,1}=StimFileList(i).Subj;
        TableW{end,2}=StimFileList(i).Session;
        if StimFileList(i).Group==1
            TableW{end,3}='Random';
        elseif StimFileList(i).Group==4
            TableW{end,3}='40Hz';
            
        else
            
        end
        TableW{end,4}=StimDay(i);
        if StimDay(i)==0|StimDay(i)==1
            TableW{end,5}='PreFlicker';
        elseif StimDay(i)==9|StimDay(i)==10
            TableW{end,5}='PostFlicker';
            
        else
            
        end
        
        
        
        TableW{end,6}=sum(CellAllProp.SubjFile==i);
        
        
        IP1=find(CellAllProp.CellType==1&CellAllProp.brainReg==1&CellAllProp.SubjFile==i);
        TableW{end,7}=length(IP1);
        IP2=find(CellAllProp.CellType==2&CellAllProp.brainReg==1&CellAllProp.SubjFile==i);
        TableW{end,8}=length(IP2);
        IP3=find(CellAllProp.CellType==1&CellAllProp.brainReg==2&CellAllProp.SubjFile==i);
        TableW{end,9}=length(IP3);
        IP4=find(CellAllProp.CellType==2&CellAllProp.brainReg==2&CellAllProp.SubjFile==i);
        TableW{end,10}=length(IP4);
        
        % %     valid=find(sum(sum(NonNanTrialNumSpa>=5,3),2)>=18&CellAllProp.FiringR>=1&ValidCellI==1);
        
        
        TableW{end,11}=length(intersect(IP1,ValidD0));
        TableW{end,12}=length(intersect(IP2,ValidD0));
        TableW{end,13}=length(intersect(IP3,ValidD0));
        TableW{end,14}=length(intersect(IP4,ValidD0));
        
    end
    
    IP1=find(CellAllProp.CellType==1&CellAllProp.brainReg==1&CellAllProp.FlickerG==4);
    IP1=intersect(IP1,ValidD0);
    tempIP=find(CellAllProp.StimDay==0|CellAllProp.StimDay==1);
    IP1=intersect(IP1,tempIP);
    length(IP1)
    % xlswrite([SavePath 'AnimalCellSummary.xlsx'],TableW,'AllData','A1:M64');
    % writetable(cell2table(TableW),[SavePath 'AnimalCellSummary'],'FileType','spreadsheet')
    
    clear TableW1;
    TableW1=TableW(2:end,:);
    TempW1=cell2table(TableW1,'VariableNames',TableW(1,:));
    % writetable(TempW1,[SavePathRaw0 'AnimalCellCrossCA'],'FileType','spreadsheet')
    
    StimGroup=[];
    for i=1:length(StimFileList)
        StimGroup(i)=StimFileList(i).Group;
    end
    
    
    clear TableW1;
    TableW1=TableW(2:end,:);
    % TableW1(InvalidFile,:)=[];
    % GenoFileW1(InvalidFile)=[];
    % xlswrite([SavePath 'AnimalCellSummary.xls'],[TableW(1,:);TableW1],'IncludedAnimal');
    % TableTemp=table2array(cell2table(TableW1(:,[1 size(TableW1,2)-6:size(TableW1,2)])));
    % AniID=unique(TableTemp(:,1));
    TableTemp=table2array(cell2table(TableW1(:,[1 6:size(TableW1,2)])));
    AniID=unique(TableTemp(:,1));
    
    
    TableW2={};
    for i=1:length(AniID)
        I1=find(TableTemp(:,1)==AniID(i));
        TableW2(end+1,1:2)=TableW1(I1(1),[1 3]);
        I2=intersect(I1,find(StimDay==0|StimDay==1));
        TableW2{end,3}='PreFlicker';
        if ~isempty(I2)
            TableW2(end,4:12)=table2cell(array2table(sum(TableTemp(I2,2:end),1)));
        end
        
        I2=intersect(I1,find(StimDay==9|StimDay==10));
        TableW2{end+1,3}='PostFlicker';
        TableW2(end,1:2)=TableW2(end-1,1:2);
        if ~isempty(I2)
            TableW2(end,4:12)=table2cell(array2table(sum(TableTemp(I2,2:end),1)));
        end
    end
    
    FlickerTemp=[];
    PrePostTemp=[];
    
    for i=1:size(TableW2,1)
        if strfind(TableW2{i,2},'40Hz')
            FlickerTemp(i)=4;
        elseif strfind(TableW2{i,2},'Random')
            FlickerTemp(i)=1;
        else
            
        end
        
        if strfind(TableW2{i,3},'Pre')
            PrePostTemp(i)=1;
        elseif strfind(TableW2{i,3},'Post')
            PrePostTemp(i)=2;
        else
            
        end
    end
    
    FType=[4 1];
    
    TableW4={};
    for i=1:length(FType)
        TempW2=table2array(cell2table(TableW2(1:end,4:end)));
        I1=find(FlickerTemp==FType(i));
        TableW3=TableW2(I1,:);
        TableW3{end+1,3}='SumPre.';
        TableW3(end,4:end)=table2cell(array2table(sum(TempW2(intersect(I1,find(PrePostTemp==1)),:))));
        TableW3{end+1,3}='SumPost';
        TableW3(end,4:end)=table2cell(array2table(sum(TempW2(intersect(I1,find(PrePostTemp==2)),:))));
        TableW4=[TableW4;TableW3];
    end
    
    TempW1=cell2table(TableW4,'VariableNames',TableW(1,[1 3 5 6:end]));
    % writetable(TempW1,[SavePathRaw0 'AnimalCellCrossCAorganize'],'FileType','spreadsheet')
    
    
    %
    % for jTh=1:length(SpikeThAll)
    %
    % SpikeTh=SpikeThAll(jTh);
    % TrialTh=TrialThAll(jTh);
    % SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
    % mkdir(SavePathRaw0);
    % %%%%%%%%%%%Set up threshold to include cells for analysis
    %
    % nBoot=10000;
    %
    % % YlimCell=[-0.00001 0.0003;-0.001 0.008;-0.00001 0.0003;-0.00001 0.0003;];
    % % YlimCellBand=[-0.00001 0.0005;-0.001 0.008;repmat([-0.00001 0.0003],3,1);repmat([-0.00001 0.00005],3,1);repmat([-0.00001 0.00002],5,1);];
    % YlimCellBand=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);];
    %
    % YtickCellBand=YlimCellBand;
    %
    % % YlimCellTheta(2,:)=[-0.001 0.01];
    % %
    % % YtickCellTheta(1,:)=[0:0.004:0.008];
    % % YtickCellTheta(2,:)=[0:0.015:0.03];
    % %
    % %
    % %
    % % YlimCellGamma(1,:)=[-0.00001 0.0002];
    % % YlimCellGamma(2,:)=[-0.00001 0.0005];
    % % YtickCellGamma(1,:)=[-0.00001:0.00025:0.0005];
    % % YtickCellGamma(2,:)=[-0.00001:0.001:0.002];
    %
    %
    % PlotColor2{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0
    % PlotColor2{1,2}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 1
    % PlotColor2{1,3}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9
    % PlotColor2{1,4}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 10
    % PlotColor2{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
    % PlotColor2{2,2}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
    % PlotColor2{2,3}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
    % PlotColor2{2,4}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
    %
    % PlotColor4=[0.698 0.8745 0.5412;0.6510 0.8078 0.8902;0.2 0.6275 0.1725;0.1216 0.4706 0.7059]; %%%%%%%40Hz Day 1
    %
    %
    % close all
    % GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
    %
    %    P.xLeft=0.03;        %%%%%%Left Margin
    %    P.xRight=0.02;       %%%%%%Right Margin
    %    P.yTop=0.02;         %%%%%%Top Margin
    %    P.yBottom=0.08;      %%%%%%Bottom Margin
    %    P.xInt=0.02;         %%%%%%Width-interval between subplots
    %    P.yInt=0.04;         %%%%%%Height-interval between subplots
    %
    %
    %
    % clear UnPairedDayFlicker
    % for iGroup=1:3
    %     if iGroup==2
    %        continue;
    %     end
    %
    %     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
    %     mkdir(SavePath);
    %
    %     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
    %     mkdir(SavePath);
    %
    %
    % % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
    % % % % length(ValidPair)
    % % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
    % % % % length(ValidPair)
    % % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
    % % % % length(ValidPair)
    %
    %
    % ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
    % length(ValidD0)
    % ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
    % length(ValidD0)
    % ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
    % length(ValidD0)
    %
    %
    %    GroupPair.CorrName='fdr';
    %    GroupPair.Test='Ranktest';
    %    GroupPair.Q=0.05;
    %    GroupPair.Pair=[1 1 2 3;2 3 4 4];
    %
    % %    GroupPair.Pair=[1 2;3 4];
    %
    %    GroupPair.SignY=0.01;
    %    GroupPair.Plot=1;
    %    GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
    %    GroupPair.SamplePlot=0; %%%%%%%%%Plot Individual Sample Point
    %    GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
    %    GroupPair.LimY=[0 GroupPair.SignY*1.2];
    %    GroupPair.Marker={'o'};
    %    GroupPair.ViolinLR=[0 1 0 1];
    %
    % % PlotColor2=[0.6 0.6 0.6;0.1 0.1 0.1];
    % PlotColor{1}=[0.6 0.6 0.6;0.8 0.1 0.2];
    % PlotColor{2}=[0.6 0.6 0.6;0.2 0.8 0.1];
    %
    %
    %
    %
    %
    % % % Param3=Param2;
    % % % Param3.TimeRepeatAnova=1;
    % % % Param3.GroupRepeatAnova=0;
    % % % Param3.RepeatAnova=1;
    % % % Param3.Paired=0;
    % % %
    % % %
    % VRname={'PreFlicker' 'PostFlicker'};
    % SaveTemp=[SavePath 'VR-40vsRandom\'];
    % mkdir(SaveTemp)
    %
    % for iCell=1:length(CellGroup)
    %
    %     tempFile=[Rfolder 'tempLME.txt'];
    %     fh=fopen(tempFile,'a');
    %
    % for iFF=1:length(FI)
    %     jFF=FI(iFF);
    %     SaveTempF=[SaveTemp];
    %     mkdir(SaveTempF);
    %
    %
    %     for iReg=1:length(RegionGroup)
    %     clear DataPlot xSubjID  xSubj DataPlotName FlickerID DayID;
%      DataPlot={}; xSubjID={};xSubj={}; DataPlotName={}; FlickerID={}; DayID={};
% for iDay=1:size(DayGroup,1)
%         
%     for iS=1:length(StimGroup) 
%         Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
%         Index2=find(CellAllProp.FlickerG==StimGroup(iS));
%         Index3=find(CellAllProp.CellType==CellGroup(iCell));
%         
%       
%         IndexNeed=intersect(intersect(Index1,Index2),Index3);
%         
%         Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
%         IndexNeed=intersect(Index4,IndexNeed);
% 
%         CellNeed=intersect(IndexNeed,ValidD0);
%            
%         if isempty(CellNeed)
%            continue; 
%         end
%         DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
%         xSubjID{end+1}=CellAllProp.SubjFile(CellNeed);
%         xSubj{end+1}=CellAllProp.Subj(CellNeed);
%         FlickerID{end+1}=repmat(iS-1,length(CellNeed),1);
%         DayID{end+1}=repmat(iDay-1,length(CellNeed),1);
%        
%         DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
% 
%    
%     end
% % %     
% end
% 
% DataLME=[];
% Cov1=[];
% Cov2=[];
% Cov3=[];
% Cov4=[];
% Cov5=[];
% 
% for iData=1:length(DataPlot)
%     DataLME=[DataLME;DataPlot{iData}(:)];
%     Cov1=[Cov1;xSubjID{iData}(:)];
%     Cov2=[Cov2;xSubj{iData}(:)];
%     Cov3=[Cov3;FlickerID{iData}(:)];
%     Cov4=[Cov4;DayID{iData}(:)];
% end
%     tbl = array2table([DataLME(:),Cov1(:),Cov2(:),Cov3(:) Cov4(:)],'VariableNames',{'PPC','Session','Animal','Flicker','Day'});
% %     lme = fitlme(tbl,'Firing ~ Speed * Spa + (Speed|CellI)');
% 
% %     lme0 = fitlme(tbl,'Firing ~ Speed + (1|CellSubj)','FitMethod','ML');
% %     lme1 = fitlme(tbl,'Firing ~ Speed + Spa + (1|CellSubj)','FitMethod','ML');
% %     lme2 = fitlme(tbl,'Firing ~ Speed * Spa -Spa+ (1|CellSubj)','FitMethod','ML');
% %     lme3 = fitlme(tbl,'Firing ~ Speed*Spa + (1|CellSubj)','FitMethod','ML');
% 
% %     results = compare(lme0,lme1);   %%%%%%%speed affect vs speed+spatial effect
% %     results = compare(lme1,lme2);   %%%%%%%speed affect vs speed*spatial effect
% %     results = compare(lme2,lme3);   %%%%%%%speed affect vs speed*spatial effect
%     
%     fprintf(fh,'\n%s\n',[RegionName{iReg} ' ' FreName2{jFF}]);
%     writetable(tbl,[Rfolder 'tempData.csv'])
% %     RunRcode([Rfolder 'LME1.R'],'C:\Program Files\R\R-4.0.2\bin')
%     RunRcode([Rfolder 'LME1FlickerAndDay.R'])
%     delete([Rfolder 'tempData.csv']);
% 
% 
%     end
% 
% 
% end
%     fclose(fh);
%     ResultFile=[SaveTemp 'LME' CellName{iCell} 'PPCvsFlickerAndDay.txt'];
%     delete(ResultFile);
%     movefile(tempFile,ResultFile);
% 
% end
% 
% end
% 
% end
% 
% 
% for jTh=1:length(SpikeThAll)
% 
% SpikeTh=SpikeThAll(jTh);
% TrialTh=TrialThAll(jTh);
% SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
% mkdir(SavePathRaw0);
% %%%%%%%%%%%Set up threshold to include cells for analysis
% 
% nBoot=10000;
% 
% % YlimCell=[-0.00001 0.0003;-0.001 0.008;-0.00001 0.0003;-0.00001 0.0003;];
% % YlimCellBand=[-0.00001 0.0005;-0.001 0.008;repmat([-0.00001 0.0003],3,1);repmat([-0.00001 0.00005],3,1);repmat([-0.00001 0.00002],5,1);];
% YlimCellBand=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);];
% 
% YtickCellBand=YlimCellBand;
% 
% % YlimCellTheta(2,:)=[-0.001 0.01];
% % 
% % YtickCellTheta(1,:)=[0:0.004:0.008];
% % YtickCellTheta(2,:)=[0:0.015:0.03];
% % 
% % 
% % 
% % YlimCellGamma(1,:)=[-0.00001 0.0002];
% % YlimCellGamma(2,:)=[-0.00001 0.0005];
% % YtickCellGamma(1,:)=[-0.00001:0.00025:0.0005];
% % YtickCellGamma(2,:)=[-0.00001:0.001:0.002];
% 
% 
% PlotColor2{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0
% PlotColor2{1,2}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 1
% PlotColor2{1,3}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9
% PlotColor2{1,4}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 10
% PlotColor2{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
% PlotColor2{2,2}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
% PlotColor2{2,3}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
% PlotColor2{2,4}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
% 
% PlotColor4=[0.698 0.8745 0.5412;0.6510 0.8078 0.8902;0.2 0.6275 0.1725;0.1216 0.4706 0.7059]; %%%%%%%40Hz Day 1
% 
% 
% close all
% GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
% 
%    P.xLeft=0.03;        %%%%%%Left Margin
%    P.xRight=0.02;       %%%%%%Right Margin
%    P.yTop=0.02;         %%%%%%Top Margin
%    P.yBottom=0.08;      %%%%%%Bottom Margin
%    P.xInt=0.02;         %%%%%%Width-interval between subplots
%    P.yInt=0.04;         %%%%%%Height-interval between subplots
% 
% 
% 
% clear UnPairedDayFlicker
% for iGroup=1:3
%     if iGroup==2
%        continue; 
%     end
%     
%     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
%     mkdir(SavePath);
%     
%     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
%     mkdir(SavePath);
% 
% 
% % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
% % % % length(ValidPair)
% % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
% % % % length(ValidPair)
% % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
% % % % length(ValidPair)
% 
% 
% ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
% length(ValidD0)
% ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
% length(ValidD0)
% ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
% length(ValidD0)
% 
% 
%    GroupPair.CorrName='fdr';
%    GroupPair.Test='Ranktest';
%    GroupPair.Q=0.05;
%    GroupPair.Pair=[1 1 2 3;2 3 4 4];
% 
% %    GroupPair.Pair=[1 2;3 4];
% 
%    GroupPair.SignY=0.01;
%    GroupPair.Plot=1;
%    GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
%    GroupPair.SamplePlot=0; %%%%%%%%%Plot Individual Sample Point
%    GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
%    GroupPair.LimY=[0 GroupPair.SignY*1.2];
%    GroupPair.Marker={'o'};
%    GroupPair.ViolinLR=[0 1 0 1];
%    
% % PlotColor2=[0.6 0.6 0.6;0.1 0.1 0.1];
% PlotColor{1}=[0.6 0.6 0.6;0.8 0.1 0.2];
% PlotColor{2}=[0.6 0.6 0.6;0.2 0.8 0.1];
% 
% 
% 
% 
% 
% % % Param3=Param2;
% % % Param3.TimeRepeatAnova=1;
% % % Param3.GroupRepeatAnova=0;
% % % Param3.RepeatAnova=1;
% % % Param3.Paired=0;
% % % 
% % % 
% VRname={'PreFlicker' 'PostFlicker'};
% SaveTemp=[SavePath 'VR-40vsRandom\'];
% mkdir(SaveTemp)
% 
% for iCell=1:length(CellGroup)
%     
%     tempFile=[Rfolder 'tempLME.txt'];
%     fh=fopen(tempFile,'a');
% 
% for iFF=1:length(FI)
%     jFF=FI(iFF);
%     SaveTempF=[SaveTemp];
%     mkdir(SaveTempF);
%     
% 
%     for iReg=1:length(RegionGroup)
%     clear DataPlot xSubjID  xSubj DataPlotName FlickerID DayID;
%      DataPlot={}; xSubjID={};xSubj={}; DataPlotName={}; FlickerID={}; DayID={};
% for iDay=1:size(DayGroup,1)
%         
%     for iS=1:length(StimGroup) 
%         Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
%         Index2=find(CellAllProp.FlickerG==StimGroup(iS));
%         Index3=find(CellAllProp.CellType==CellGroup(iCell));
%         
%       
%         IndexNeed=intersect(intersect(Index1,Index2),Index3);
%         
%         Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
%         IndexNeed=intersect(Index4,IndexNeed);
% 
%         CellNeed=intersect(IndexNeed,ValidD0);
%            
%         if isempty(CellNeed)
%            continue; 
%         end
%         DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
%         xSubjID{end+1}=CellAllProp.SubjFile(CellNeed);
%         xSubj{end+1}=CellAllProp.Subj(CellNeed);
%         FlickerID{end+1}=repmat(iS-1,length(CellNeed),1);
%         DayID{end+1}=repmat(iDay-1,length(CellNeed),1);
%        
%         DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
% 
%    
%     end
% % %     
% end
% 
% DataLME=[];
% Cov1=[];
% Cov2=[];
% Cov3=[];
% Cov4=[];
% Cov5=[];
% 
% for iData=1:length(DataPlot)
%     DataLME=[DataLME;DataPlot{iData}(:)];
%     Cov1=[Cov1;xSubjID{iData}(:)];
%     Cov2=[Cov2;xSubj{iData}(:)];
%     Cov3=[Cov3;FlickerID{iData}(:)];
%     Cov4=[Cov4;DayID{iData}(:)];
% end
%     tbl = array2table([DataLME(:),Cov1(:),Cov2(:),Cov3(:) Cov4(:)],'VariableNames',{'PPC','Session','Animal','Flicker','Day'});
% %     lme = fitlme(tbl,'Firing ~ Speed * Spa + (Speed|CellI)');
% 
% %     lme0 = fitlme(tbl,'Firing ~ Speed + (1|CellSubj)','FitMethod','ML');
% %     lme1 = fitlme(tbl,'Firing ~ Speed + Spa + (1|CellSubj)','FitMethod','ML');
% %     lme2 = fitlme(tbl,'Firing ~ Speed * Spa -Spa+ (1|CellSubj)','FitMethod','ML');
% %     lme3 = fitlme(tbl,'Firing ~ Speed*Spa + (1|CellSubj)','FitMethod','ML');
% 
% %     results = compare(lme0,lme1);   %%%%%%%speed affect vs speed+spatial effect
% %     results = compare(lme1,lme2);   %%%%%%%speed affect vs speed*spatial effect
% %     results = compare(lme2,lme3);   %%%%%%%speed affect vs speed*spatial effect
%     
%     fprintf(fh,'\n%s\n',[RegionName{iReg} ' ' FreName2{jFF}]);
%     writetable(tbl,[Rfolder 'tempData.csv'])
% %     RunRcode([Rfolder 'LME1.R'],'C:\Program Files\R\R-4.0.2\bin')
%     RunRcode([Rfolder 'LME2FlickerAndDay.R'])
%     delete([Rfolder 'tempData.csv']);
% 
% 
%     end
% 
% 
% end
%     fclose(fh);
%     ResultFile=[SaveTemp 'LME' CellName{iCell} 'PPCvsFlickerAndDaySession.txt'];
%     delete(ResultFile);
%     movefile(tempFile,ResultFile);
% 
% end
% 
% end
% 
% end
% 
% 
% % %%%%%%%%%%%SingleSubj Plot, Hierarchy Bootstrap
% % for jTh=1:length(SpikeThAll)
% % 
% % SpikeTh=SpikeThAll(jTh);
% % TrialTh=TrialThAll(jTh);
% % SavePathRaw0=[SavePathRaw 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
% % mkdir(SavePathRaw0);
% % %%%%%%%%%%%Set up threshold to include cells for analysis
% % 
% % nBoot=10000;
% % 
% % % YlimCell=[-0.00001 0.0003;-0.001 0.008;-0.00001 0.0003;-0.00001 0.0003;];
% % % YlimCellBand=[-0.00001 0.0005;-0.001 0.008;repmat([-0.00001 0.0003],3,1);repmat([-0.00001 0.00005],3,1);repmat([-0.00001 0.00002],5,1);];
% % YlimCellBand=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0002],3,1);]*0.6;
% % 
% % YtickCellBand=YlimCellBand;
% % 
% % % YlimCellTheta(2,:)=[-0.001 0.01];
% % % 
% % % YtickCellTheta(1,:)=[0:0.004:0.008];
% % % YtickCellTheta(2,:)=[0:0.015:0.03];
% % % 
% % % 
% % % 
% % % YlimCellGamma(1,:)=[-0.00001 0.0002];
% % % YlimCellGamma(2,:)=[-0.00001 0.0005];
% % % YtickCellGamma(1,:)=[-0.00001:0.00025:0.0005];
% % % YtickCellGamma(2,:)=[-0.00001:0.001:0.002];
% % 
% % 
% % PlotColor2{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0
% % PlotColor2{1,2}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 1
% % PlotColor2{1,3}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9
% % PlotColor2{1,4}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 10
% % PlotColor2{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
% % PlotColor2{2,2}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
% % PlotColor2{2,3}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
% % PlotColor2{2,4}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
% % 
% % PlotColor4=[0.698 0.8745 0.5412;0.6510 0.8078 0.8902;0.2 0.6275 0.1725;0.1216 0.4706 0.7059]; %%%%%%%40Hz Day 1
% % 
% % 
% % close all
% % GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};
% % 
% %    P.xLeft=0.03;        %%%%%%Left Margin
% %    P.xRight=0.02;       %%%%%%Right Margin
% %    P.yTop=0.02;         %%%%%%Top Margin
% %    P.yBottom=0.08;      %%%%%%Bottom Margin
% %    P.xInt=0.02;         %%%%%%Width-interval between subplots
% %    P.yInt=0.04;         %%%%%%Height-interval between subplots
% % 
% % 
% % 
% % clear UnPairedDayFlicker
% % for iGroup=1:3
% %     if iGroup==2
% %        continue; 
% %     end
% %     
% %     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
% %     mkdir(SavePath);
% %     
% %     SavePath=[SavePathRaw0 GroupName{iGroup} '\'];
% %     mkdir(SavePath);
% % 
% % 
% % % % % ValidPair=find(sum(isnan(squeeze(PPCvr(:,:,iGroup))),2)==0);   %%%%%%Only correct trial
% % % % % length(ValidPair)
% % % % % ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,iGroup)>=SpikeTh,2)==1));
% % % % % length(ValidPair)
% % % % % ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,iGroup)>=TrialTh,2)==1));
% % % % % length(ValidPair)
% % 
% % 
% % ValidD0=find(sum(isnan(squeeze(PPCvr(:,1,iGroup))),2)==0);   %%%%%%Only correct trial
% % length(ValidD0)
% % ValidD0=intersect(ValidD0,find(TSVRNum(:,iGroup)>=SpikeTh));
% % length(ValidD0)
% % ValidD0=intersect(ValidD0,find(NonNanVRNum(:,iGroup)>=TrialTh));
% % length(ValidD0)
% % 
% % 
% %    GroupPair.CorrName='fdr';
% %    GroupPair.Test='Ranktest';
% %    GroupPair.Q=0.05;
% %    GroupPair.Pair=[1 1 2 3;2 3 4 4];
% % 
% % %    GroupPair.Pair=[1 2;3 4];
% % 
% %    GroupPair.SignY=0.01;
% %    GroupPair.Plot=1;
% %    GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
% %    GroupPair.SamplePlot=0; %%%%%%%%%Plot Individual Sample Point
% %    GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
% %    GroupPair.LimY=[0 GroupPair.SignY*1.2];
% %    GroupPair.Marker={'o'};
% %    GroupPair.ViolinLR=[0 1 0 1];
% %    
% % % PlotColor2=[0.6 0.6 0.6;0.1 0.1 0.1];
% % PlotColor{1}=[0.6 0.6 0.6;0.8 0.1 0.2];
% % PlotColor{2}=[0.6 0.6 0.6;0.2 0.8 0.1];
% % 
% % 
% % 
% % 
% % 
% % % % Param3=Param2;
% % % % Param3.TimeRepeatAnova=1;
% % % % Param3.GroupRepeatAnova=0;
% % % % Param3.RepeatAnova=1;
% % % % Param3.Paired=0;
% % % % 
% % % % 
% % VRname={'PreFlicker' 'PostFlicker'};
% % SaveTemp=[SavePath 'VR-40vsRandom\'];
% % mkdir(SaveTemp)
% % 
% % for iCell=1:length(CellGroup)
% %     figure;
% % for iFF=1:length(FI)
% %     jFF=FI(iFF);
% %     SaveTempF=[SaveTemp];
% %     mkdir(SaveTempF);
% %     
% %     YlimCell=repmat(YlimCellBand(jFF,:),2,1);
% %     YtickCell=repmat(YtickCellBand(jFF,:),2,1);
% % % %     if jFF<=2
% % % %         YlimCell=YlimCellTheta;
% % % %         YtickCell=YtickCellTheta;
% % % %     else
% % % %         YlimCell=YlimCellGamma;
% % % %         YtickCell=YtickCellGamma;
% % % % 
% % % %     end
% 
%     for iReg=1:length(RegionGroup)
%      subplotLU(length(RegionGroup),length(FI),iReg,iFF,P)
%     GroupPair.SignY=YlimCell(iCell,2);
%     GroupPair.LimY=[0 YlimCell(iCell,2)];
%     clear DataPlot xSubjID DataPlotName;
%      DataPlot={}; xSubjID={}; DataPlotName={};
% for iDay=1:size(DayGroup,1)
%         
%     for iS=1:length(StimGroup) 
%         Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
%         Index2=find(CellAllProp.FlickerG==StimGroup(iS));
%         Index3=find(CellAllProp.CellType==CellGroup(iCell));
%         
%       
%         IndexNeed=intersect(intersect(Index1,Index2),Index3);
%         
%         Index4=find(CellAllProp.brainReg==RegionGroup(iReg));
%         IndexNeed=intersect(Index4,IndexNeed);
% 
%         CellNeed=intersect(IndexNeed,ValidD0);
%            
%         if isempty(CellNeed)
%            continue; 
%         end
%         DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,iGroup));
%         xSubjID{end+1}=CellAllProp.Subj(CellNeed);
%         DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];
%         
%         PlotColor3(iS,:)=[PlotColor2{iS,iDay}(2,:)];
% 
%    
%     end
% % %     
% end
% 
% % %     subplotLU(1,length(DayGroup),iCell,iDay);
%     GroupPair.Marker={'o','^','o','^'};
% 
%     PathSave=[SaveTempF RegionName{iReg} CellName{iCell} FreName{jFF} '.txt'];
%     Datatype=0;
%     GroupPair.SamplePairedPlot=0;
%     x=[1 1.4 3 3.4];
% 
% 
%     PathSave=[SaveTempF RegionName{iReg} CellName{iCell} FreName{jFF} 'HY.txt'];
%     VRUnPaired40VsRandom{iCell,iDay,iFF}=ErrorBarHierarchy(x,DataPlot,xSubjID,PlotColor4,nBoot,Datatype,PathSave,GroupPair,[1 2 3 4]);
% %     VRUnPaired40VsRandom{iCell,iDay,iFF}=ErrorBoxPlotLU(1:2,DataPlot,PlotColor3,[],GroupPair,[1 2]);
%     set(gca,'xlim',[0 4],'xtick',x,'ylim',YlimCell(iCell,:),'ytick',YtickCell(iCell,:),'xticklabel',[]);
%     if iFF==1
%     ylabel([RegionName{iReg} ' PPC']);
%     end
%     %%%%PosXY(1,:): x coordinate;PosXY(2,:): y coordinate; coordinate of plot
%     %%%%PosXY(3,:): x coordinate;PosXY(4,:): y coordinate; coordinate of labels
%     
%         if iFF==1&&iReg==1
% 
%     yTemp=YtickCell(iCell,2);
%     tickStep=mean(diff(YtickCell(iCell,:)));
%     yTemp=[yTemp yTemp-tickStep*0.1];
%     tempXY=[repmat(x(2)+0.1,1,2);yTemp;repmat(x(2)+0.3,1,2);yTemp];
%     LuLegend(tempXY,0,DataPlotName(1:2),PlotColor4(1:2,:));
% 
%     tempXY=[repmat(x(4)+0.1,1,2);yTemp;repmat(x(4)+0.3,1,2);yTemp];
%     LuLegend(tempXY,0,DataPlotName(3:4),PlotColor4(3:4,:));
%         end
%                 if iReg==2
%                    xlabel(FreName2{jFF});
%                 end
% %     LuFontStandard;
% 
%     end
% 
% 
% end
%     papersizePX=[0 0 6*length(FI) 6*length(RegionGroup)];
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlickerHY'],'tiff');
% %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker.eps'],'epsc'); 
% %     saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'fig'); 
%     close all
% 
% end
% % % save([SaveTemp 'statis.mat'],'VRUnPaired40VsRandom');
% % % close all
% 
% end
% 
% end
