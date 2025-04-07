clear all
LoadData='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step1\';
load([LoadData 'StimInfoList.mat'],'StimFileList');
SavePathLoad='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step4\';
load([SavePathLoad '6bandVRhilbertPPCresultsLocalCAVRShankColPyrChNotch40.mat'],'PPCvrFile','NonNanVRNumFile','TSppcVRNumFile');
load([SavePathLoad '6bandflickerhilbertPPCresultsLocalCAShankColPyrChNotch40.mat']);
% % SavePathRaw=[SavePathLoad '13bandHirearchy\'];
% % mkdir(SavePathRaw)
%%%%%%%%%%Integrate VR results
SavePathRaw=[SavePathLoad 'Manuscript\'];
mkdir(SavePathRaw)
SavePathRaw=[SavePathRaw '6bandCALocalShankColVRFlickerNotch40\'];
mkdir(SavePathRaw)
%%%%%%%%%%Integrate VR results




clear CellAllProp
NonNanVRNum=[];
TSVRNum=[];
PPCvr=[];
CellAllProp.CellType=[];
CellAllProp.SubjFile=[];
CellAllProp.brainReg=[];
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
PPCvr=[];
NeuroTemp=NeuronData;
PPCFiletemp=PPCvrFile;
TSppcVRNumTemp=TSppcVRNumFile;
NonNanVRNumTemp=NonNanVRNumFile;



NonNanVRNum=[];
TSVRNum=[];

% % FI=2;
for i=1:length(NeuronData)

    if i==40
    Iv=isnan(NeuronData(i).CellType);
    if ~isempty(Iv)
        i
    PPCFiletemp{i}(Iv,:,:,:)=[];
    TSppcVRNumTemp{i}(Iv,:,:)=[];
    NonNanVRNumTemp{i}(Iv,:,:)=[];

    NeuroTemp(i).CellType(Iv)=[];
    NeuroTemp(i).brainReg(Iv)=[];
    NeuroTemp(i).TSIDtotal(Iv)=[];

    end



    end
    PPCvr=cat(1,PPCvr,PPCFiletemp{i});
    NonNanVRNum=cat(1,NonNanVRNum,NonNanVRNumTemp{i});
    TSVRNum=cat(1,TSVRNum,TSppcVRNumTemp{i});

    
    CellAllProp.CellType=[CellAllProp.CellType;NeuroTemp(i).CellType];
    CellAllProp.brainReg=[CellAllProp.brainReg;NeuroTemp(i).brainReg];
    CellAllProp.TSIDtotal=[CellAllProp.TSIDtotal;NeuroTemp(i).TSIDtotal];
    CellAllProp.SubjFile=[CellAllProp.SubjFile;repmat(i,length(NeuroTemp(i).TSIDtotal),1)];
    CellAllProp.Subj=[CellAllProp.Subj;repmat(StimFileList(i).Subj,length(NeuroTemp(i).TSIDtotal),1)];


end
%%%%%%%%%%Integrate VR results

%%%%%%%%%%Integrate Flicker results
PPCflicker=[];
NeuroTemp=NeuronData;
PPCFiletemp=PPCflickerFile;
TSppcflickerNumTemp=TSppcflickerNumFile;
NonNanflickerNumTemp=NonNanflickerNumFile;



NonNanflickerNum=[];
TSflickerNum=[];
for i=1:length(NeuronData)

    if i==40
    Iv=isnan(NeuronData(i).CellType);
    if ~isempty(Iv)
        i
    NeuroTemp(i).CellType(Iv)=[];
    PPCFiletemp{i}(Iv,:,:,:)=[];
    TSppcflickerNumTemp{i}(Iv,:,:)=[];
    NonNanflickerNumTemp{i}(Iv,:,:)=[];
    end
    end
    PPCflicker=cat(1,PPCflicker,PPCFiletemp{i});
    NonNanflickerNum=cat(1,NonNanflickerNum,NonNanflickerNumTemp{i});
    TSflickerNum=cat(1,TSflickerNum,TSppcflickerNumTemp{i});

    
end
%%%%%%%%%%Integrate Flicker results
%%%%%%%%%%Integrate VR results


% 
% 
% 
% 
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

% 
% 
% 
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


%%%%%%%%%%%currently using 100 spikes, 5 trials threshold 12/09/2022 For %%%%%%%%%%%Abby's manuscript
SpikeThAll=[100];
TrialThAll=[5];




Rfolder='Y:\singer\LuZhang\Project2ChronicF\Analysis\RFolder\';
%%%%%%%%%%%%LME:Only Flicker period, random vs 40, pre vs post, Sepreate Cell,Region
for jTh=1:length(SpikeThAll)

SpikeTh=SpikeThAll(jTh);
TrialTh=TrialThAll(jTh);
SavePathRaw00=[SavePathRaw '2DayMergeFlickerCrossDay\'];
mkdir(SavePathRaw00)
SavePathRaw0=[SavePathRaw00 'FlickerVR' num2str(SpikeTh) 'Spike' num2str(TrialTh) 'Trial\'];
mkdir(SavePathRaw0);
%%%%%%%%%%%Set up threshold to include cells for analysis

nBoot=10000;

YlimCellTheta(1,:)=[-0.001 0.008];
YlimCellTheta(2,:)=[-0.001 0.01];

YtickCellTheta(1,:)=[0:0.004:0.008];
YtickCellTheta(2,:)=[0:0.015:0.03];



YlimCellGamma(1,:)=[-0.00001 0.0002];
YlimCellGamma(2,:)=[-0.00001 0.0005];
YtickCellGamma(1,:)=[-0.00001:0.00025:0.0005];
YtickCellGamma(2,:)=[-0.00001:0.001:0.002];

YlimCellBand=[-0.00001 0.0008;-0.001 0.008;repmat([-0.00001 0.0008],2,1);repmat([-0.00001 0.0003],3,1);];
YtickCellBand=YlimCellBand;

PlotColor2{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0
PlotColor2{1,2}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 1
PlotColor2{1,3}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9
PlotColor2{1,4}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 10
PlotColor2{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
PlotColor2{2,2}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0
PlotColor2{2,3}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1
PlotColor2{2,4}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 1

PlotColor3{1,1}=[0.6 0.6 0.6;0.698 0.8745 0.5412;0.1 0.1 0.1]; %%%%%%%StimRandom Day 0 and 1
PlotColor3{1,2}=[0.6 0.6 0.6;0.2 0.6275 0.1725;0.1 0.1 0.1]; %%%%%%%StimRandom Day 9 and 10
PlotColor3{2,1}=[0.6 0.6 0.6;0.6510 0.8078 0.8902;0.1 0.1 0.1]; %%%%%%%40Hz Day 0 and 1
PlotColor3{2,2}=[0.6 0.6 0.6;0.1216 0.4706 0.7059;0.1 0.1 0.1]; %%%%%%%40Hz Day 9 and 10

PlotColor4=[0.698 0.8745 0.5412;0.6510 0.8078 0.8902;0.2 0.6275 0.1725;0.1216 0.4706 0.7059]; %%%%%%%40Hz Day 1


close all
GroupName={'AllTrial','IncorrectTrial','CorrectTrial'};

   P.xLeft=0.03;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.08;         %%%%%%Top Margin
   P.yBottom=0.08;      %%%%%%Bottom Margin
   P.xInt=0.05;         %%%%%%Width-interval between subplots
   P.yInt=0.08;         %%%%%%Height-interval between subplots



clear UnPairedDayFlicker
    
SavePath=[SavePathRaw0];
mkdir(SavePath);


% ValidPair=find(sum(sum(isnan(squeeze(PPCvr(:,:,:,iGroup))),3),2)==0);   %%%%%%Only correct trial
% length(ValidPair)
% ValidPair=intersect(ValidPair,find(sum(TSVRNum(:,:,iGroup)>=SpikeTh,2)==2));
% length(ValidPair)
% ValidPair=intersect(ValidPair,find(sum(NonNanVRNum(:,:,iGroup)>=TrialTh,2)==2));
% length(ValidPair)
% 
% 
% ValidD0=find(sum(isnan(squeeze(PPCvr(:,:,1,iGroup))),2)==0);   %%%%%%Only correct trial
% length(ValidD0)
% ValidD0=intersect(ValidD0,find(TSVRNum(:,1,iGroup)>=SpikeTh));
% length(ValidD0)
% ValidD0=intersect(ValidD0,find(NonNanVRNum(:,1,iGroup)>=TrialTh));
% length(ValidD0)


ValidFlicker=find(sum(isnan(PPCflicker),2)==0);   %%%%%%Only correct trial
length(ValidFlicker)
ValidFlicker=intersect(ValidFlicker,find(TSflickerNum>=SpikeTh));
length(ValidFlicker)
ValidFlicker=intersect(ValidFlicker,find(NonNanflickerNum>=TrialTh));
length(ValidFlicker)




   
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




SaveTemp=SavePath;
mkdir(SaveTemp);

SaveTempANOVA=[SaveTemp 'MixedAnova\'];
mkdir(SaveTempANOVA)





% % Param3=Param2;
% % Param3.TimeRepeatAnova=1;
% % Param3.GroupRepeatAnova=0;
% % Param3.RepeatAnova=1;
% % Param3.Paired=0;
% % 
% % 
for iCell=1:length(CellGroup)
    tempFile=[Rfolder 'tempLME.txt'];
    fh=fopen(tempFile,'a');

    figure;
for iFF=1:length(FI)
    jFF=FI(iFF);
    SaveTempF=[SaveTemp];
    mkdir(SaveTempF);
% %     if jFF<=2
% %         YlimCell=YlimCellTheta;
% %         YtickCell=YtickCellTheta;
% %     else
% %         YlimCell=YlimCellGamma;
% %         YtickCell=YtickCellGamma;
% % 
% %     end
    YlimCell=repmat(YlimCellBand(jFF,:),2,1);
    YtickCell=repmat(YtickCellBand(jFF,:),2,1);

    for iReg=1:length(RegionGroup)
     subplotLU(length(RegionGroup),length(FI),iReg,iFF,P)
    GroupPair.SignY=YlimCell(iCell,2)*0.9;
    GroupPair.LimY=[0 YlimCell(iCell,2)];
%     clear DataPlot xSubjID DataPlotName FlickerID DayID x colorX;
    clear DataPlot xSubjID  xSubj DataPlotName FlickerID DayID StimDayID CellID StimDSeq x colorX;
     DataPlot={}; xSubjID={};xSubj={}; DataPlotName={}; FlickerID={}; DayID={}; StimDayID={};CellID={};StimDSeq={};xSubj={};StimDay={};colorX=[];

     x=[];x0=-1;Groupx=[];CountX=1;GroupPair.Pair=[];
for iDay=1:length(DayGroup)
    for iS=1:length(StimGroup) 
        Index1=find(CellAllProp.StimDay==DayGroup(iDay,1)|CellAllProp.StimDay==DayGroup(iDay,2));
        Index2=find(CellAllProp.FlickerG==StimGroup(iS));
        Index3=find(CellAllProp.CellType==CellGroup(iCell));
        Index4=find(CellAllProp.brainReg==RegionGroup(iReg));

        IndexNeed=intersect(intersect(Index1,Index2),Index3);
        IndexNeed=intersect(IndexNeed,Index4);
     
%         CellNeed=intersect(IndexNeed,ValidPair);
        CellNeed=intersect(IndexNeed,ValidFlicker);
   
     if isempty(CellNeed)
        continue; 
     end
%    clear DataPlot xSubjID
%     DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,1,iGroup));
%     xSubjID{1}=CellAllProp.SubjFile(CellNeed);
    
    DataPlot{end+1}=squeeze(PPCflicker(CellNeed,jFF));
%     xSubjID{2}=CellAllProp.SubjFile(CellNeed);
    xSubjID{end+1}=CellAllProp.SubjFile(CellNeed);
    xSubj{end+1}=CellAllProp.Subj(CellNeed);
    FlickerID{end+1}=repmat(iS-1,length(CellNeed),1);
    DayID{end+1}=repmat(iDay-1,length(CellNeed),1);
    CellID{end+1}=CellAllProp.CellID(CellNeed);
    StimDayID{end+1}=CellAllProp.StimDay(CellNeed);
    StimDSeq{end+1}=CellAllProp.StimDSeq(CellNeed);
        
    DataPlotName{end+1}=[DayName{iDay} ' ' StimName{iS} ', n =' num2str(length(DataPlot{end}))];

%     DataPlot{end+1}=squeeze(PPCvr(CellNeed,jFF,2,iGroup));
%     xSubjID{3}=CellAllProp.SubjFile(CellNeed);
% % %     subplotLU(length(StimGroup),length(DayGroup),iS,iDay);
    colorX=[colorX;PlotColor3{iS,iDay}];
    Groupx=[Groupx zeros(1,3)+CountX];
%     tempPair=[1 1 2;2 3 3];
%     GroupPair.Pair=[GroupPair.Pair tempPair+(CountX-1)*3];
    GroupPair.Pair=[1 1 2 3;2 3 4 4];

    CountX=CountX+1;

    
    Datatype=1;
    x=[x;[1:3]'+x0+1];
    x0=x(end);
%     PathSave=[SaveTempF 'bandPPC' StimName{iS} CellName{iCell} DayName{iDay} '.txt'];
%     WithinDayComparison{iS,iCell,iDay,iFF}=ErrorBarHierarchy(1:3,DataPlot,xSubjID,PlotColor2{iS,iDay},nBoot,Datatype,PathSave,GroupPair,[1 1 1]);
%     WithinDayComparison{iS,iCell,iDay,iFF}=ErrorBoxPlotLU(1:3,DataPlot,PlotColor2{iS,iDay},[],GroupPair,[1 1 1]);  
    end
    
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
    MixAnova2F1W1B_withR(DataLME,[Cov4(:) Cov3(:) Cov2(:)],[SaveTempANOVA RegionName{iReg} CellName{iCell} FreName2{jFF}]);


    tbl = array2table([DataLME(:),Cov1(:),Cov2(:),Cov3(:) Cov4(:) Cov5(:) Cov6(:) Cov7(:)],'VariableNames',{'PPC','Session','Animal','Flicker','Day','StimDay','Cell','StimDSeq'});
    writetable(tbl,[Rfolder 'tempData.csv'])
    fprintf(fh,'\n%s\n',[RegionName{iReg} ' ' FreName2{jFF}]);

    RunRcode([Rfolder 'RplotAnimalDayFlicker.R'],'C:\Program Files\R\R-4.0.2\bin')
% % %     RunRcode([Rfolder 'RplotAnimalDay.R'])
    delete([Rfolder 'tempData.csv']);
    
    
    ResultFigureFile=[SaveTemp RegionName{iReg} CellName{iCell} FreName{jFF} 'AnimalDay.png'];
    %delete(ResultFigureFile);
    tempFigure=[Rfolder 'PPCAnimalDay.png'];
    movefile(tempFigure,ResultFigureFile);
    PathSave=[SaveTempF 'bandFlicker' RegionName{iReg} CellName{iCell} FreName{jFF} '.txt'];
%     WithinDayComparison{iS,iCell,iDay,iFF}=ErrorBarHierarchy(1:3,DataPlot,xSubjID,PlotColor2{iS,iDay},nBoot,Datatype,PathSave,GroupPair,[1 1 1]);
%     WithinDayComparison{iReg,iFF}=ErrorBoxPlotLU(x,DataPlot,colorX,PathSave,GroupPair,Groupx);  
%     set(gca,'xlim',[0 x(end)+1],'xtick',x,'ylim',YlimCell(iCell,:),'ytick',YtickCell(iCell,:),'xticklabel',[]);

% %     GroupPair.Marker={'o','^','o','^'};
% %     PathSave=[SaveTempF RegionName{iReg} CellName{iCell} FreName{jFF} '.txt'];
% %     Datatype=0;
% %     GroupPair.SamplePairedPlot=0;
% % %     VRUnPaired40VsRandom{iCell,iDay,iFF}=ErrorBarHierarchy(1:2,DataPlot,xSubjID,PlotColor3,nBoot,Datatype,PathSave,GroupPair,[1 2]);
% % %     VRUnPaired40VsRandom{iCell,iDay,iFF}=ErrorBoxPlotLU(1:2,DataPlot,PlotColor3,[],GroupPair,[1 2]);
% %     x=[1 1.2 3 3.2];
    x=[1 1.4 3 3.4];

    UnPairedDayFlicker{iCell,iReg,iFF}=ErrorViolinHalf(x,DataPlot,PlotColor4,0,PathSave,GroupPair,[1 2 3 4]);
    set(gca,'xlim',[0 4],'xtick',x,'ylim',YlimCell(iCell,:),'ytick',YtickCell(iCell,:),'xticklabel',[]);
    if iFF==1
    ylabel([RegionName{iReg} ' PPC']);
    end
    %%%%PosXY(1,:): x coordinate;PosXY(2,:): y coordinate; coordinate of plot
    %%%%PosXY(3,:): x coordinate;PosXY(4,:): y coordinate; coordinate of labels
    
        if iFF==2

    yTemp=YtickCell(iCell,2);
    tickStep=mean(diff(YtickCell(iCell,:)));
    yTemp=[yTemp yTemp-tickStep*0.05];
    tempXY=[repmat(x(1),1,2);yTemp;repmat(x(1)+0.2,1,2);yTemp];
    LuLegend(tempXY,0,DataPlotName(1:2),PlotColor4(1:2,:),6);

%     yTemp=[yTemp yTemp-tickStep*0.1];
    tempXY=[repmat(x(3),1,2);yTemp;repmat(x(3)+0.2,1,2);yTemp];
    LuLegend(tempXY,0,DataPlotName(3:4),PlotColor4(3:4,:),6);
%     LuFontStandard;
        end
         if iReg==2
            xlabel(FreName2{jFF});
         end
%     LuFontStandard;

    end


end
    papersizePX=[0 0 6*length(FI) 6*length(RegionGroup)];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'pdf');
    saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'png');
    saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker.eps'],'epsc'); 
    saveas(gcf,[SaveTempF CellName{iCell} '40RandomPrePostFlicker'],'fig'); 
    close all

    fclose(fh);
    ResultFile=[SaveTemp 'LME' CellName{iCell} 'PPCvsFlickerAndDay.txt'];
% %     delete([SaveTemp 'LME*.txt']);
    movefile(tempFile,ResultFile);

end
    save([SaveTemp 'statis.mat'],'UnPairedDayFlicker');
% % close all


end
%%%%%%%%%%%%LME:Only Flicker period,random vs 40, pre vs post, Sepreate Cell,Region









%%%%%%
SpikeTh=100;
TrialTh=5;
iGroup=3;    %%%%%%%%Only Correct trial
StimDay=StimDayAll;



ValidFlicker=find(sum(isnan(PPCflicker),2)==0);   %%%%%%Only correct trial
length(ValidFlicker)
ValidFlicker=intersect(ValidFlicker,find(TSflickerNum>=SpikeTh));
length(ValidFlicker)
ValidFlicker=intersect(ValidFlicker,find(NonNanflickerNum>=TrialTh));
length(ValidFlicker)




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
    

    TableW{end,11}=length(intersect(IP1,ValidFlicker));
    TableW{end,12}=length(intersect(IP2,ValidFlicker));
    TableW{end,13}=length(intersect(IP3,ValidFlicker));
    TableW{end,14}=length(intersect(IP4,ValidFlicker));

end

IP1=find(CellAllProp.CellType==1&CellAllProp.brainReg==1&CellAllProp.FlickerG==4);
IP1=intersect(IP1,ValidFlicker);
tempIP=find(CellAllProp.StimDay==0|CellAllProp.StimDay==1);
IP1=intersect(IP1,tempIP);
length(IP1)
% xlswrite([SavePath 'AnimalCellSummary.xlsx'],TableW,'AllData','A1:M64');
% writetable(cell2table(TableW),[SavePath 'AnimalCellSummary'],'FileType','spreadsheet')

clear TableW1;
TableW1=TableW(2:end,:);
TempW1=cell2table(TableW1,'VariableNames',TableW(1,:));
writetable(TempW1,[SavePathRaw0 'AnimalCellLocalCA'],'FileType','spreadsheet')

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
writetable(TempW1,[SavePathRaw0 'AnimalCellLocalCAorganize'],'FileType','spreadsheet')


