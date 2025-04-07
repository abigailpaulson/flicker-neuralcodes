clear all

KeyWord1='NFFT1024';
KeyWord2='Win1024';

load(['Y:\singer\LuZhang\Project2ChronicF\Results\step2\NoSpeedThwpiPSD_PreGoal\' KeyWord1 KeyWord2 'SubResultsVR.mat'],'COHtype','PSDregion');
Data(1).wpli=COHtype{2}    %%%WPLI
Data(1).PSDregion=PSDregion;
Data(1).Name='VR';
clear COHtype PSDregion
load(['Y:\singer\LuZhang\Project2ChronicF\Results\step2\NoSpeedThwpiPSD_Flicker\' KeyWord1 KeyWord2 'SubResultsflicker.mat']);
Data(2).wpli=COHtype{2}    %%%WPLI
Data(2).PSDregion=PSDregion;
Data(2).Name='Flicker';
SavePath='Y:\singer\LuZhang\Project2ChronicF\Results\step2\Manuscript\';
SavePath=['Y:\singer\LuZhang\Project2ChronicF\Results\step2\Manuscript\' KeyWord1 KeyWord2 '\'];
mkdir(SavePath)


iFband=1;
FLimBand=[1 120];

%%
pwpli=[];
ppsd=[];
pwpliRank=[];
ppsdRank=[];
for iData=1:length(Data)
     

    FLim=FLimBand(iFband,:);
    FTick=0:50:200;
    FI=find(Fre>=FLim(1)&Fre<=FLim(2));
    FrePlot=Fre(FI);

COHPlot=abs(Data(iData).wpli(:,FI));
ParamCoh.Bin=FrePlot;
ParamCoh.TimeCol=FrePlot;

clear DataPlot
    for iStim=1:length(StimGroup)
        NeedI=find(StimEL==1&StimType==StimGroup(iStim));
        DataPlot{iStim}=(abs(COHPlot(NeedI,:)));
        DataPlot{iStim}(sum(isnan(DataPlot{iStim}),2)>0,:)=[];
    end

    ip=4;
    ParamCoh.PathSave=[SavePath 'wpli' Data(iData).Name CompareN{ip}];
    WPLIstats{iData}=RateHist_Stats(DataPlot,ParamCoh);
    pwpli=[pwpli;WPLIstats{iData}.ttestP.Ppair(:)];
    pwpliRank=[pwpliRank;WPLIstats{iData}.ttestP.PpairRank(:)];
    for iReg=1:2
clear DataPlot
       PSDnorm=Data(iData).PSDregion{iReg};
       PSDnorm=abs(PSDnorm(:,FI));
        A=nansum(PSDnorm(:,:),2);
        A=repmat(A,1,length(FrePlot),1,1); 
% A=repmat(A,1,length(PSDTaskFile(1).Fre),1,1); 

       PSDnorm=PSDnorm./A;
       ParamCoh.PathSave=[SavePath 'psd' Data(iData).Name CompareN{ip} RegionName{iReg}];

      for iStim=1:length(StimGroup)
        NeedI=find(StimEL==1&StimType==StimGroup(iStim));
        DataPlot{iStim}=(abs(PSDnorm(NeedI,:)));
        DataPlot{iStim}(sum(isnan(DataPlot{iStim}),2)>0,:)=[];
      end
        PSDstats{iData,iReg}=RateHist_Stats(DataPlot,ParamCoh);
        ppsd=[ppsd;PSDstats{iData,iReg}.ttestP.Ppair(:)];
        ppsdRank=[ppsdRank;PSDstats{iData}.ttestP.PpairRank(:)];
    end
    ip=4;



end



%% 
Q=0.1
[h, wpli_crit_p, adj_p]=fdr_bh(pwpli,Q,'pdep','yes');
wpli_crit_p
[h, wpli_crit_p, adj_p]=fdr_bh([pwpli;ppsd(:)],Q,'pdep','yes');
wpli_crit_p

[h, wpli_crit_p, adj_p]=fdr_bh(pwpliRank,Q,'pdep','yes');
wpli_crit_p
[h, wpli_crit_p, adj_p]=fdr_bh([pwpliRank;ppsdRank(:)],Q,'pdep','yes');
wpli_crit_p


P.xLeft=0.15;
P.xRight=0.02;
P.yTop=0.02;
P.yBottom=0.15;
P.xInt=0.1;
P.yInt=0.05;
PlotColor2=[0.698 0.8745 0.5412;0.2 0.6275 0.1725;0.6510 0.8078 0.8902;0.1216 0.4706 0.7059];

ParamCoh.SigPlot='Ttest'
figure;
for iData=1:length(Data)
     

    FLim=FLimBand(iFband,:);
    FTick=0:50:200;
    FI=find(Fre>=FLim(1)&Fre<=FLim(2));
    FrePlot=Fre(FI);

COHPlot=abs(Data(iData).wpli(:,FI));
ParamCoh.Bin=FrePlot;
ParamCoh.TimeCol=FrePlot;

clear DataPlot
    for iStim=1:length(StimGroup)
        NeedI=find(StimEL==1&StimType==StimGroup(iStim));
        DataPlot{iStim}=(abs(COHPlot(NeedI,:)));
        DataPlot{iStim}(sum(isnan(DataPlot{iStim}),2)>0,:)=[];
    end

    ip=4;
%     ParamCoh.PathSave=[SavePath 'wpli' Data(iData).Name CompareN{ip}];
    subplotLU(2,3,iData,1,P);
    IndexCom=ComparePair(:,ip);
    ParamTemp=ParamCoh;
    ParamTemp.PathSave=[];
    ParamTemp.Crit_p=0.05;
    ParamTemp.Ytick=[0:0.2:1];

    ParamTemp.ANOVAstats=WPLIstats{iData};
    RateHist_GroupPlot(FrePlot,DataPlot,PlotColor2(IndexCom,:),ParamTemp);
    set(gca,'xtick',[0:30:120],'ytick',0:0.25:1,'xlim',[0 120])
    if iData==1
    set(gca,'xticklabel',[]);
end

    for iReg=1:2
clear DataPlot
       PSDnorm=Data(iData).PSDregion{iReg};
       PSDnorm=abs(PSDnorm(:,FI));
        A=nansum(PSDnorm(:,:),2);
        A=repmat(A,1,length(FrePlot),1,1); 
% A=repmat(A,1,length(PSDTaskFile(1).Fre),1,1); 

       PSDnorm=PSDnorm./A;

      for iStim=1:length(StimGroup)
        NeedI=find(StimEL==1&StimType==StimGroup(iStim));
        DataPlot{iStim}=(abs(PSDnorm(NeedI,:)));
        DataPlot{iStim}(sum(isnan(DataPlot{iStim}),2)>0,:)=[];
      end
            ParamTemp=ParamCoh;
            ParamTemp.PathSave=[];
            ParamTemp.Ytick=[0:0.2:1];
            ParamTemp.Crit_p=0.05;
            ParamTemp.ANOVAstats=PSDstats{iData,iReg};
            subplotLU(2,3,iData,1+iReg,P);

            RateHist_GroupPlot(FrePlot,DataPlot,PlotColor2(IndexCom,:),ParamTemp);
            set(gca,'yscale','log','ylim',[0.0001 1],'ytick',[0.0001 0.001 0.01 0.1 1],'xtick',[0:30:120],'xlim',[0 120])
if iData==1
    set(gca,'xticklabel',[]);
end
    end


end
LuFontStandard

papersizePX=[0 0 18 12];

set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

saveas(gcf,[SavePath 'WPLIpsd'],'png'); 
saveas(gcf,[SavePath 'WPLIpsd'],'svg'); 


length(FrePlot)
sum(WPLIstats{1}.ttestP.Ppair(:)'<0.05&FrePlot(:)'>=25&FrePlot(:)'<=55)
sum(WPLIstats{1}.ttestP.Ppair(:)'<0.05)

sum(WPLIstats{2}.ttestP.Ppair(:)'<0.05&FrePlot(:)'>=25&FrePlot(:)'<=55)
sum(WPLIstats{2}.ttestP.Ppair(:)'<0.05)



%% Mixed ANOVA
gSubj=[];gFre=[];gStim=[];gSubj=[];gEL=[];gVRflicker=[];DataAnova=[];
iFband=1;
FLim=FLimBand(iFband,:);



Rfolder='Y:\singer\LuZhang\Project2ChronicF\Analysis\RFolder\';

tempFile=[Rfolder 'tempLME.txt'];
fh=fopen(tempFile,'a');


for iData=1:2
WPLIAll=Data(iData).wpli;
DataAnovaTemp=WPLIAll;
FI=find(Fre>=FLim(1)&Fre<=FLim(2));
DataAnovaTemp=WPLIAll(:,FI);
FrePlot=Fre(FI);
nFre=length(FrePlot);
gSubjTemp=repmat(Subj(:),1,nFre);
gFreTemp=repmat(FrePlot(:)',length(Subj),1);
gELTemp=repmat(StimEL(:),1,nFre);
gStimTemp=repmat(StimType(:),1,nFre);
gSubjTemp(isnan(DataAnovaTemp))=[];
gFreTemp(isnan(DataAnovaTemp))=[];
gStimTemp(isnan(DataAnovaTemp))=[];
gELTemp(isnan(DataAnovaTemp))=[];
DataAnovaTemp(isnan(DataAnovaTemp))=[];

    tbl = array2table([DataAnovaTemp(:),gSubjTemp(:),gStimTemp(:),gFreTemp(:) gELTemp(:)],'VariableNames',{'WPLI','Animal','Flicker','Fre','Day'});
%   lme = fitlme(tbl,'Firing ~ Speed * Spa + (Speed|CellI)');
%   lme0 = fitlme(tbl,'Firing ~ Speed + (1|CellSubj)','FitMethod','ML');
%   lme1 = fitlme(tbl,'Firing ~ Speed + Spa + (1|CellSubj)','FitMethod','ML');
%   lme2 = fitlme(tbl,'Firing ~ Speed * Spa -Spa+ (1|CellSubj)','FitMethod','ML');
%   lme3 = fitlme(tbl,'Firing ~ Speed*Spa + (1|CellSubj)','FitMethod','ML');

%   results = compare(lme0,lme1);   %%%%%%%speed affect vs speed+spatial effect
%   results = compare(lme1,lme2);   %%%%%%%speed affect vs speed*spatial effect
%   results = compare(lme2,lme3);   %%%%%%%speed affect vs speed*spatial effect

    writetable(tbl,[Rfolder 'tempData.csv'])
    fprintf(fh,'\n%s\n',[Data(iData).Name]);

    RunRcode([Rfolder 'LME8wpliSepDay.R'],'C:\Program Files\R\R-4.0.2\bin');
    
    delete([Rfolder 'tempData.csv']);
    



DataAnova=[DataAnova(:);DataAnovaTemp(:)];
gSubj=[gSubj(:);gSubjTemp(:)];
gFre=[gFre(:);gFreTemp(:)];
gStim=[gStim(:);gStimTemp(:)];
gEL=[gEL(:);gELTemp(:)];
gVRflicker=[gVRflicker(:);zeros(length(DataAnovaTemp),1)+iData];
end
fclose(fh)

    ResultFile=[SavePath 'LMEwpli_FreUnFactorsBothFixRandom.txt']
    %delete(ResultFigureFile);
    movefile(tempFile,ResultFile);


PostFlickerI=find(gEL==1);
ResultFile=MixAnova3F2W1B_withR(DataAnova(PostFlickerI),[gFre(PostFlickerI) gVRflicker(PostFlickerI) gStim(PostFlickerI) gSubj(PostFlickerI)],[SavePath 'WPLI']);




