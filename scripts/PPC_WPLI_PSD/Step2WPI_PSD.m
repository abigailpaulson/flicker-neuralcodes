%% CA3-CA1 LFP WPLI(Weighted phase lag index) and PSD during VR period

%% load data and set parameters
clear all
LoadData='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step1\';
load([LoadData 'StimInfoList.mat'],'StimFileList');
%%%%%%%Channels chosen by best Col among 5 Cols of one Shank
% load([LoadData 'DetermineCA1BoardPyrRunList.mat'],'RunList');
load('Y:\singer\LuZhang\Project2ChronicF\Results\step1\ThetaAlign\ShankColAlignCA1PyrOriCA3PyrCh.mat','ThetaList');
RunList=ThetaList; clear ThetaList;



SavePath='Y:\singer\LuZhang\Project2ChronicF\Results\step2\';
mkdir(SavePath);
SavePath=[SavePath 'NoSpeedThwpiPSD_PreGoal\'];
mkdir(SavePath);
SavePathSubj=[SavePath 'Subj\'];
mkdir(SavePathSubj);

%%  This section loads the calculated results, the saved results was previoiusly done at line 558
% load([SavePath 'ShColAlignWPLIFile.mat'],'PSDca1All','PSDca3All','CohAll','WPLIAll','CohGroup','psdParameter','SavePath','FLimBand','statis','CompareN','ComparePair');
%
%%
% mkdir(SavePath)f
ChList=0:63;   %%%%ChID for data storing
% 
% 
samprate=2000;% 
DownR=4;
samprateD=samprate/DownR;

psdParameter.Fs=samprateD;
psdParameter.window=512;
psdParameter.noverlap=0;
psdParameter.nfft=512;
nFre=psdParameter.nfft/2+1;

% WaveParam.Fband=[1;150];
% WaveParam.MorseParam=[3 9 4];
% WaveParam.samprate=samprateD;
% WaveParam.DownSample=1;
% WaveParam.NSW=3;
% WaveParam.NTW=11;
% WaveParam.flag_SMOOTH=1;


% % Notch60 = designfilt('bandstopiir','FilterOrder',2, ...
% %                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
% %                'DesignMethod','butter','SampleRate',samprate);
           
% Notch60 = designfilt('bandstopiir','FilterOrder',2, ...
%                'HalfPowerFrequency1',59.5,'HalfPowerFrequency2',60.5, ...
%                'SampleRate',samprate);
% Notch60 = designfilt('bandstopiir','FilterOrder',20, ...
%                'StopbandFrequency1',59.5,'StopbandFrequency2',60.5, ...
%                'SampleRate',samprate,'StopbandAttenuation',0.1);
       
           
% fvtool(Notch60)


%% Calculate CA3-CA1 cross-spectrum, and CA3/CA1 power spectrum

%%%%%%%%%%%%%%consider all trial
for i=1:length(RunList)
    SaveTemp=[SavePathSubj 'Subj' num2str(i)];
    mkdir(SaveTemp);

        tic
        
        clear temp2CA1 temp2CA3
           temp1CA1=[RunList(i).PathFile 'CA1\'];
           ChPyrCA1=ChList(RunList(i).CA1PyrChIndex);
           for iCh=1:length(ChPyrCA1)
           temp2CA1{iCh}=[temp1CA1 num2str(ChPyrCA1(iCh)) '\'];    %%%%%%%%CA1 Shank iCh
           end

           temp1CA3=[RunList(i).PathFile 'CA3\'];
           ChPyrCA3=ChList(RunList(i).CA3PyrChIndex);
           for iCh=1:length(ChPyrCA3)
           temp2CA3{iCh}=[temp1CA3 num2str(ChPyrCA3(iCh)) '\'];    %%%%%%%%CA3 Shank iCh
           end

           LFPChAll{1}=cell(1,length(ChPyrCA1));
           LFPChAll{2}=cell(1,length(ChPyrCA3));
        

              for iReg=1:length(LFPChAll) 
                 for iSh=1:length(LFPChAll{iReg})
                     LFPDataTrial{iReg}{iSh}=struct('Data',[],'Time',[]);
                     LFPDataTrial{iReg}{iSh}(1)=[];
                 end
                 MarkTrialI=0;
              end
           
   temp3CA1=[temp1CA1 'sorted\'];
   temp3CA3=[temp1CA3 'sorted\'];

   if     isempty(dir(temp3CA1))&&isempty(dir(temp3CA3))
          continue
   end

    if StimFileList(i).StimDay(1)==0
       VRtaskFileID=find(StimFileList(i).RType==0);  %%%%%For day 0 with no flickers, using all VR files.
    elseif StimFileList(i).StimDay(1)==1             %%%%%For day 1, VR files before flicker period is used.
       VRtaskFileID=find(StimFileList(i).RType==1);
    elseif StimFileList(i).StimDay(1)==9||StimFileList(i).StimDay(1)==10      %%%%For day 9 and 10, Merge VR from pre and post flicker period.
       VRtaskFileID=find(StimFileList(i).RType==1|StimFileList(i).RType==3);
    else
       
    end
    FlickerFileID=find(StimFileList(i).RType==2);


      for jFile=1:min(length(VRtaskFileID),length(StimFileList(i).Trial))
                 iFile=VRtaskFileID(jFile);



                CountF=0;    %%%%%%%whether count the file or not?0 not count yet;1 counted already
                TotalZone=[StimFileList(i).Trial(iFile).Zone1Type StimFileList(i).Trial(iFile).Zone2Type];
%                 TotalZoneTrialTs=[StimFileList(i).Trial(iFile).Zone1TrialTs StimFileList(i).Trial(iFile).Zone2TrialTs];
                TotalZoneTrialTs=[StimFileList(i).Trial(iFile).Zone1PreTs StimFileList(i).Trial(iFile).Zone2PreTs];
 

                if isempty(TotalZoneTrialTs)
                   continue
                end
                  [~,sortI]=sort(TotalZoneTrialTs(1,:));
                  TotalZoneTrialTs=TotalZoneTrialTs(:,sortI);
                  TotalZone=TotalZone(sortI);
                  clear sortI

                  temptype1=find(TotalZone==0);  %%%%%%%%0 no reward.
                  temptype2=find(TotalZone==1);  %%%%%%%%1 get reward;
                  
                 
%                 PeriodGroup{1}=MergePeriod(TotalZoneTrialTs);    %%0 no reward
%                 PeriodGroup{2}=MergePeriod(TotalZoneTrialTs(:,temptype1));    %%0 no reward
%                 PeriodGroup{3}=MergePeriod(TotalZoneTrialTs(:,temptype2));    %%1 get reward
%                 PeriodTask=MergePeriod(TotalZoneTrialTs(:,temptype2));
                  PeriodTask=TotalZoneTrialTs;  %%All trial considered
                  
                  PeriodTask(:,diff(PeriodTask)<2)=[];
                  
                  if isempty(PeriodTask)
                     continue
                  end
            
                  PosData=StimFileList(i).Trial(iFile).PosData;
                  if isempty(PosData)
                     continue; 
                  end
                   
                  clear LFPTempCA1
                 for iSh=1:length(ChPyrCA1)
                     if length(dir(temp2CA1{iSh}))>=1
                           load([temp2CA1{iSh} 'eeg' num2str(RunList(i).File(iFile))]);
                           LFPTempCA1{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);

                     end
                 end
                 
                 clear LFPTempCA3
                 for iSh=1:length(ChPyrCA3)
                     if length(dir(temp2CA3{iSh}))>=1
                         load([temp2CA3{iSh} 'eeg' num2str(RunList(i).File(iFile))]);
                         LFPTempCA3{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);

                     end
                 end
               

                
                 LFPTempCA{1}=LFPTempCA1;
                 LFPTempCA{2}=LFPTempCA3;
%                  clear LFPTempCA1 LFPTempCA3;
                 Time=[0:length(LFPTempCA1{1})-1]/samprateD;
                     for iTrial=1:size(PeriodTask,2)
                           tempP1=PeriodTask;
                           PeriodUsedI=round(tempP1*samprateD);
                           PeriodUsedI(PeriodUsedI<1)=1;
                           PeriodUsedI(PeriodUsedI>length(LFPTempCA1{1}))=length(LFPTempCA1{1});
                     end

                 for iReg=1:length(LFPTempCA) 
                   for iSh=1:length(LFPTempCA{iReg})
                     for iTrial=1:size(PeriodTask,2)
                          LFPDataTrial{iReg}{iSh}(MarkTrialI+iTrial,1).Data=LFPTempCA{iReg}{iSh}(PeriodUsedI(1,iTrial):PeriodUsedI(2,iTrial));
                          LFPDataTrial{iReg}{iSh}(MarkTrialI+iTrial,1).Time=Time(PeriodUsedI(1,iTrial):PeriodUsedI(2,iTrial));          
                     end
                   end
                 end
                MarkTrialI=MarkTrialI+size(PeriodTask,2);
        end
        clear TempTrial1 TempTrial2
            if isempty(LFPDataTrial{1}{iSh})
               continue
            end
            if isempty(LFPDataTrial{2}{iSh})
               continue
            end

        for iSh=1:length(LFPTempCA{1})
            TempTrial1=LFPDataTrial{1}{iSh};
           for jSh=1:length(LFPTempCA{2})
               TempTrial2=LFPDataTrial{2}{jSh};
                [TrialSpec.Sxy,TrialSpec.Sxx,TrialSpec.Syy,TrialSpec.w,TrialSpec.options,ValidIndex]=crossspec_NonEqualTriL(TempTrial1,TempTrial2,psdParameter);          
                CohGroup{i}{iSh,jSh}=crossspec_NonEqualTriL_TrialIndex(TrialSpec,ValidIndex);
% %                 a=crossspec_NonEqualTriL_TrialIndex(TrialSpec,ValidIndex);
% %                 figure;
% %                 plot(a.Fre,smooth(abs(a.wpli1),21))
% %                 a=crossspec_NonEqualTriL_TrialIndex(TrialSpec,ValidIndex);
% %                 figure;
% %                 plot(a.Fre,abs(a.wpli2))
% %                 figure;
% %                 plot(a.Fre,abs(a.wpli3))
% % 
%                  save([SaveTemp 'CA1Sh' num2str(iSh) 'CA3' num2str(jSh) '.mat'],'TrialSpec','ValidIndex','psdParameter');


%                  [TrialSpec.Sxy,TrialSpec.Sxx,TrialSpec.Syy,TrialSpec.Fre,ValidIndex]=crossspecMorse_NonEqualTriL(TempTrial1,TempTrial2,samprateD,WaveParam);


%                  save([SaveTemp 'Ch' num2str(iCh) 'Ch' num2str(jCh) '.mat'],'TrialSpec','ValidIndex','psdParameter');
%                  
%                 a=crossspecMorse_NonEqualTriL_TrialIndex(TrialSpec,ValidIndex);
%                 figure;
%                 plot(a.Fre,abs(a.wpli1))
%                 figure;
%                 plot(a.Fre,abs(a.wpli2))
%                 figure;
%                 plot(a.Fre,abs(a.wpli3))


           end
        end
     toc   


end

%% Calculate and Integrate WPLI and psd

clear PSDca1All PSDca3All  CohAll WPLIAll
for i=1:length(CohGroup)
    PxxCA1=[];
    for iSh=1:size(CohGroup{i},1)
        PxxCA1(iSh,:)=CohGroup{i}{iSh,1}.Pxx;
    end
    if size(PxxCA1,1)==1
       PSDca1All(i,:)=PxxCA1;
    else
       PSDca1All(i,:)=nanmean(PxxCA1);
    end

    PxxCA3=[];
    for iSh=1:size(CohGroup{i},2)
        PxxCA3(iSh,:)=CohGroup{i}{1,iSh}.Pyy;
    end
    if size(PxxCA3,1)==1
       PSDca3All(i,:)=PxxCA3;
    else
       PSDca3All(i,:)=nanmean(PxxCA3);
    end

    Cohtemp=[];
    WPLItemp=[];
    for iSh=1:size(CohGroup{i},1)
         for jSh=1:size(CohGroup{i},2)
             Cohtemp(end+1,:)=CohGroup{i}{iSh,jSh}.Cxy;
             WPLItemp(end+1,:)=CohGroup{i}{iSh,jSh}.wpli1;
         end
    end
    if size(Cohtemp,1)==1
    CohAll(i,:)=Cohtemp;
    WPLIAll(i,:)=WPLItemp;
    else
    CohAll(i,:)=nanmean(Cohtemp);
    WPLIAll(i,:)=nanmean(WPLItemp);
    end

end



i=1;
iSh=1;
Fre=CohGroup{i}{iSh,1}.Fre;
FLimBand=[1 160];    
FI=find(Fre>=FLimBand(1)&Fre<=FLimBand(2));
FrePlot=Fre(FI);


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

[SubjCandi,I]=unique(Subj);
SubjCandi-Subj(I);
% % FlickerCandi=FlickerT(I);



Invalid=[];
FlickerTG=[];

for iSubj=1:length(SubjCandi)
    I=find(Subj==SubjCandi(iSubj));
% %     iSubj
% %     StimDay(I)
    if ~isempty(setdiff([0 1 9 10],StimDay(I)))
       Invalid=[Invalid;iSubj];
    end
    FlickerTG(iSubj)=FlickerT(max(I));
end
% 
StimG{1}=setdiff(find(FlickerTG==1),Invalid);
StimG{2}=setdiff(find(FlickerTG==4),Invalid);

StimType=[];
StimEL=[];
for i=1:length(Subj)
    StimType(i)=FlickerTG(find(SubjCandi==Subj(i)));
    if ~isempty(intersect(StimDay(i),[0 1]))
       StimEL(i)=0;   %%%%%Pre chronic flicker
    else
       StimEL(i)=1;   %%%%%Post chronic flicker
    end
        
end
StimGroup=[1 4]; %%random or 40Hz
ELGroup=[0 1]; %%Pre or post flicker


Param2.PlotType=2;
Param2.SigPlot='Anova';
Param2.SigPlot='Ttest';
Param2.CorrName='fdr';   %%%methold for multi-compairson
Param2.Q=0.1;
% Param2.Ytick=[0 0.002 0.004];
Param2.LegendShow=0;
Param2.Legend=[];
Param2.TimeRepeatAnova=1;
Param2.GroupRepeatAnova=1;
Param2.RepeatAnova=1;
Param2.TimeCol=Fre;
Param2.Paired=1;
Param2.BinName='Fre';
Param2.Bin=Fre;
Param2.TimeComparison=0;
Param2.statisP=1;
Param2.Ytick=[0:0.1:0.2];

ylimC=[repmat([-0.0005 0.004],3,1);-0.0005 0.006];
ylimC=[repmat([-0.0005 0.003],3,1);-0.0005 0.006];
ylimC=[0.1 0.8;0.1 0.8];
ylimC=[-10 -1;-10 -1];

ylimCIn=[0.3 0.6;0.3 0.6];

PlotColor2=[0.698 0.8745 0.5412;0.2 0.6275 0.1725;0.6510 0.8078 0.8902;0.1216 0.4706 0.7059];




P.xLeft=0.01;
P.xRight=0.01;
P.yTop=0.01;
P.yBottom=0.01;
P.xInt=0.02;
P.yInt=0.04;

SessionI=[1 2 3];    %%%%%%Wrong trial and Correct trial
%%recently used in Paper
Param2.PlotType=2;
Param2.SigPlot='Anova';
% Param2.SigPlot='Ttest';
Param2.CorrName='fdr';   %%%methold for multi-compairson
Param2.Q=0.1;
% Param2.Ytick=[0 0.002 0.004];
Param2.LegendShow=0;
Param2.Legend=[];
Param2.TimeRepeatAnova=1;
Param2.GroupRepeatAnova=1;
Param2.RepeatAnova=1;
Param2.TimeCol=FrePlot;
Param2.Paired=1;
Param2.BinName='Fre';
Param2.Bin=FrePlot;
Param2.TimeComparison=0;
Param2.statisP=1;
Param2.Ytick=0:0.2:0.8;

mkdir([SavePath 'statis\']);
% delete([SavePath 'statis\*.txt']);


ylimC=[-10 -1;-10 -1];


Param3=Param2;


% %     Param3.ANOVAstats=statisAll(ip,1);
    Param3.TimeRepeatAnova=1;
    Param3.GroupRepeatAnova=0;
    Param3.RepeatAnova=1;
%     Param3.ANOVAstats=statisAll(ip,1);
    Param3.Paired=0;
    Param3.PathSave=[];
   
    Param3.PlotType=6;
    Param3.Crit_p=0.05;
%   Param3.ANOVAstats=genoWR(iC,iG);
    Param3.statisP=1;
%     Param3=rmfield(Param3,'ANOVAstats');
    Param3.Marker={'none','none'};

    ComparePair=[1 3 1 2 ;2 4 3 4];
    CompareN={'RandomPrePost' '40PrePost' 'PreFlicker' 'PostFlicker' };
    
   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.06;      %%%%%%Bottom Margin
   P.xInt=0.04;         %%%%%%Width-interval between subplots
   P.yInt=0.02;         %%%%%%Height-interval between subplots
    
    
  

PSDregion{1}=PSDca1All;
PSDregion{2}=PSDca3All;

RegionName{1}='CA1';
RegionName{2}='CA3';

ParamCoh=Param3;
for iFband=1:size(FLimBand,1)   
    figure;
for iReg=1:length(PSDregion)
     

FLim=FLimBand(iFband,:);
FTick=0:50:200;
Fre=CohGroup{1}{1}.Fre;
FI=find(Fre>=FLim(1)&Fre<=FLim(2));
FrePlot=Fre(FI);


PlotColor2=[0.698 0.8745 0.5412;0.2 0.6275 0.1725;0.6510 0.8078 0.8902;0.1216 0.4706 0.7059];

PSDnorm=abs(PSDregion{iReg}(:,FI));
A=nansum(PSDnorm(:,:),2);
A=repmat(A,1,length(FrePlot),1,1); 
% A=repmat(A,1,length(PSDTaskFile(1).Fre),1,1); 

PSDnorm=PSDnorm./A;
ParamCoh.Bin=FrePlot;
ParamCoh.TimeCol=FrePlot;

clear DataPlot
DataPlot={};
for iStim=1:length(StimGroup)
    for iPre=1:length(ELGroup)
        NeedI=find(StimEL==ELGroup(iPre)&StimType==StimGroup(iStim));
        DataPlot{end+1}=(abs(PSDnorm(NeedI,:)));
    end
end

    for ip=1:length(CompareN)
        subplotLU(2,length(CompareN),iReg,ip,P);
        IndexCom=ComparePair(:,ip);
        ParamCoh.PathSave=[SavePath 'statis\PSD' CompareN{iReg} CompareN{ip} 'FreNorm' num2str(FLim(1)) '-' num2str(FLim(2))];
        [~,statis{iFband,ip,iReg}]=RateHist_GroupPlot(FrePlot,DataPlot(IndexCom),PlotColor2(IndexCom,:),ParamCoh);
        set(gca,'yscale','log','ylim',[0.00001 1])

    end
 papersizePX=[0 0 26 12];


    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    saveas(gcf,[SavePath 'PSD' 'FreNorm' num2str(FLim(1)) '-' num2str(FLim(2))],'png'); 

end
end  



COHtype{1}=CohAll;
COHtype{2}=WPLIAll;
COHName{1}='Coh';
COHName{2}='WPLI';
ParamCoh=Param3;
ParamCoh.Ytick=[0:0.25:1];
ParamCoh.PlotType=3;

figure;
for iType=1:length(COHtype)
     

FLim=FLimBand(iFband,:);
FTick=0:50:200;
Fre=CohGroup{1}{1}.Fre;
FI=find(Fre>=FLim(1)&Fre<=FLim(2));
FrePlot=Fre(FI);


PlotColor2=[0.698 0.8745 0.5412;0.2 0.6275 0.1725;0.6510 0.8078 0.8902;0.1216 0.4706 0.7059];

COHPlot=abs(COHtype{iType}(:,FI));
ParamCoh.Bin=FrePlot;
ParamCoh.TimeCol=FrePlot;

clear DataPlot
DataPlot={};
for iStim=1:length(StimGroup)
    for iPre=1:length(ELGroup)
        NeedI=find(StimEL==ELGroup(iPre)&StimType==StimGroup(iStim));
        DataPlot{end+1}=(abs(COHPlot(NeedI,:)));
        DataPlot{end}(sum(isnan(DataPlot{end}),2)>0,:)=[];
    end
end

    for ip=1:length(CompareN)
        subplotLU(2,length(CompareN),iType,ip,P);
        IndexCom=ComparePair(:,ip);
        ParamCoh.PathSave=[SavePath 'statis\' COHName{iType} CompareN{ip}];
        [~,statisCOH{ip,iType}]=RateHist_GroupPlot(FrePlot,DataPlot(IndexCom),PlotColor2(IndexCom,:),ParamCoh);
        set(gca,'ylim',[0.0 1.1])

    end
    papersizePX=[0 0 26 6];

    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%     saveas(gcf,[SavePath 'PSD' 'FreNorm' num2str(FLim(1)) '-' num2str(FLim(2))],'png'); 

end
papersizePX=[0 0 26 12];

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

saveas(gcf,[SavePath 'COHandWPLI'],'png'); 

save([SavePath 'NFFT' num2str(psdParameter.nfft) 'Win' num2str(psdParameter.window) 'ShColAlignWPLIFile.mat'],'PSDca1All','PSDca3All','CohAll','WPLIAll','CohGroup','psdParameter','SavePath','FLimBand','statis','CompareN','ComparePair');




save([SavePath 'NFFT' num2str(psdParameter.nfft) 'Win' num2str(psdParameter.window) 'SubResultsVR.mat'],'COHtype','COHName','PSDregion','RegionName','Fre','CompareN','ComparePair','ParamCoh','PlotColor2');



%% Mixed ANOVA
% iFband=1;
% FLim=FLimBand(iFband,:);
% DataAnova=WPLIAll;
% FI=find(Fre>=FLim(1)&Fre<=FLim(2));
% DataAnova=WPLIAll(:,FI);
% FrePlot=Fre(FI);
% nFre=length(FrePlot);
% gSubj=repmat(Subj(:),1,nFre);
% gFre=repmat(FrePlot(:)',length(Subj),1);
% gEL=repmat(StimEL(:),1,nFre);
% gStim=repmat(StimType(:),1,nFre);
% gSubj(isnan(DataAnova))=[];
% gFre(isnan(DataAnova))=[];
% gStim(isnan(DataAnova))=[];
% gEL(isnan(DataAnova))=[];
% DataAnova(isnan(DataAnova))=[];
% ResultFile=MixAnova3F2W1B_withR(DataAnova(:),[gFre(:) gEL(:) gStim(:) gSubj(:)],[SavePath 'statis\WPLI']);

%% 
% FOI=[6 10 25 60;10 20 55 100];
% FreName={'Theta','Beta','Sgamma','Mgamma'}
% for i=1:size(FreName,2)
%     tempFI=find(Fre>=FOI(1,i)&Fre<=FOI(2,i));
%     WPLIband(:,i)=nanmean(WPLIAll(:,tempFI),2);
% end
% 
% figure;
%      
% iType=2;
% FLim=FLimBand(iFband,:);
% FTick=0:50:200;
% Fre=CohGroup{1}{1}.Fre;
% FI=find(Fre>=FLim(1)&Fre<=FLim(2));
% FrePlot=Fre(FI);
% 
% 
% PlotColor2=[0.698 0.8745 0.5412;0.2 0.6275 0.1725;0.6510 0.8078 0.8902;0.1216 0.4706 0.7059];
% 
% 
% ParamBand=Param3;
% 
% COHPlot=WPLIband;
% ParamBand.Bin=1:4;
% ParamBand.TimeCol=1:4;
% 
% clear DataPlot
% DataPlot={};
% for iStim=1:length(StimGroup)
%     for iPre=1:length(ELGroup)
%         NeedI=find(StimEL==ELGroup(iPre)&StimType==StimGroup(iStim));
%         DataPlot{end+1}=(abs(COHPlot(NeedI,:)));
%         DataPlot{end}(sum(isnan(DataPlot{end}),2)>0,:)=[];
%     end
% end
% 
% ip=4;iType=2;
% 
% % subplotLU(2,length(CompareN),iType,ip,P);
% figure;
%         IndexCom=ComparePair(:,ip);
%         ParamBand.PathSave=[SavePath 'statis\WPLIband' COHName{iType} CompareN{ip}];
%         [~,statisBand{ip,iType}]=RateHist_GroupPlot(1:4,DataPlot(IndexCom),PlotColor2(IndexCom,:),ParamCoh);
%         set(gca,'ylim',[0.0 1.1])
% 
%     papersizePX=[0 0 26 6];
% 
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     saveas(gcf,[SavePath 'PSD' 'FreNorm' num2str(FLim(1)) '-' num2str(FLim(2))],'png'); 
% 
% 
% 
% 
% 
% A=struct2table(StimFileList)
% B=array2table(WPLIband,'VariableNames',FreName);
% C=[A(:,1:2) B];
% C.Session=num2str(round(C.Session));
% 
% writetable(C,[SavePath 'wpliVR.csv'])
