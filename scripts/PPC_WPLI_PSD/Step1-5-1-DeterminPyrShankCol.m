%% Determine the best LFP channels, mainly for CA1 region
% Fit ripple power by gaussian curve for channels on the same column of a shank
% Fit thetaphase offset fiiting with sigmond function for channels on the same column of a shank

clear all
load(['Y:\singer\LuZhang\Project2ChronicF\Results\step1\DeterminePyrRunList.mat'],'RunList');
RippleList=RunList;
% load(['Y:\singer\LuZhang\Project2ChronicF\Results\step1\ThetaAlign\ShankAlign.mat']);
load(['Y:\singer\LuZhang\Project2ChronicF\Results\step1\AllData.mat']);

ThetaList=ThetaAlignList;
for i=1:length(ThetaList)
    temp1=ThetaList(i).PathFile; 
    s1=findstr(temp1,'singer');
    ThetaList(i).PathFile=['Y:\' temp1(s1:end)];
end

ThetaList=rmfield(ThetaList,'CA1ThChIndex');
ThetaList=rmfield(ThetaList,'CA1PyrChIndex');
% ThetaList=rmfield(ThetaList,'PyrThChIndex');



SaveData='Y:\singer\LuZhang\Project2ChronicF\Results\step1\ThetaAlign\';

ColorCol=colormap("jet");
ColorCol=ColorCol(1:4:64,:);


samprate=2000;
MWidth=samprate./[12];      %%%%%%%%[12]Hz range between peaks;

ComColor= [10,103,155;176,34,242;78,211,34;191,126,0;178,51,51;51,204,178]/255;
PlotTimeIn=[19;22]*samprate;
[~,ChsortI]=sort(ProbeXY(:,2),'descend');
% [~,ChsortI]=sort(ProbeXY(:,2),'ascend');

[~,Edges]=discretize([0:280],28);
[C,Edges]=discretize(ProbeXY(ChsortI,1),Edges);
GroupX=unique(C);

for iX=1:length(GroupX)
   ChsortII{iX}=ChsortI(find(C==GroupX(iX)));
end

   ShankID=[1 2];
   ShankCh=[ones(32,1);ones(32,1)+1];


[~,EdgesSh]=discretize(ShankID,2);

[C,EdgesSh]=discretize(ShankCh(ChsortI,1),EdgesSh);
GroupXSh=unique(C);

for iX=1:length(GroupXSh)
   ChsortSh{iX}=ChsortI(find(C==GroupXSh(iX)));
end
   
   
    
   P.xLeft=0.15;
   P.xRight=0.02;
   P.yTop=0.02;
   P.yBottom=0.08;
   P.xInt=0.02;
   P.yInt=0.02;
%    P.Clim=[0.0000 0.0015];
%    P.Xlim=[-pi pi];
%    P.Xtick=[-pi 0 pi];
%    P.Xticklabel={'-\pi','0','\pi'};
%    P.Ylim=[20 180];
%    P.Ytick=[20 100 180 ];
%    P.Yticklabel={'20' '100' '180'};
%    P.Ylabel='Frequency Hz';
%    P.Xlabel='Theta phase rad';


   PosMatX=[0.9;0.9];
   PosMatY=[0.55;0.35];

   PlotCol=[-0.1 0 0.1];

   P.xLeft=0.01;
   P.xRight=0.01;
   P.xInt=0.002;
   P.yTop=0.01;
   P.yBottom=0.01;
   P.yInt=0.01;

   RippleVecFit=zeros(5,length(ThetaList),length(ChsortII))+nan;
   ThetaVecFit=zeros(5,length(ThetaList),length(ChsortII))+nan;
   RippleAcc=zeros(length(ThetaList),length(ChsortII))+nan;
   ThetaAcc=zeros(length(ThetaList),length(ChsortII))+nan;

   xf=-10:200;
   Rippleyf=zeros(length(xf),length(ThetaList),length(ChsortII))+nan;
   Thetayf=zeros(length(xf),length(ThetaList),length(ChsortII))+nan;

%%%%%%%%%%%
% figure;
for i=1:length(ThetaList)
%         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
%         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','ticmes new roman','horizontalalignment','left','verticalalignment','top');
%       ThetaList(i).MaxLag=[];
       if isempty(ThetaList(i).CA1CCAll)
          continue
       end

for iCol=1:length(ChsortII)
%    subplotLU(length(ChsortII),length(ThetaList),iCol,i,P);
% %    set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4],'ylim',[0 0.006]);

   if isempty(ThetaList(i).CA1CCAll{iCol})
      continue;
   end

%     plot(ProbeXY(ChsortII{iCol},2),RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500,'r.');
    [RippleVecFit(:,iCol,i),ChRipple(i,iCol),xf,Rippleyf(:,i,iCol)]=ripplePowerVsDepth(ChsortII{iCol},RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500,ProbeXY(ChsortII{iCol},2));
    RippleAcc(i,iCol)=RippleVecFit(1,iCol,i)/std(RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500);
    
% %     hold on;
% %     plot(xf,Rippleyf(:,i,iCol),'r-');
% %      text(xf(20),Rippleyf(20),['rmse' num2str(RippleAcc(i,iCol))],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')

%     text(xf(20),yf(20),['rmse' num2str(RippleVecFit(1,iCol,i))],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
    
% %    figure;
%     hold on;
    [~,ThetaList(i).MaxLag(ChsortII{iCol})]=max(ThetaList(i).CA1CCAll{iCol});
    if std(ThetaList(i).MaxLag(ChsortII{iCol}))<0.00000000001
       continue;
        
    end

    
%     plot(ProbeXY(ChsortII{iCol},2),ThetaList(i).MaxLag(ChsortII{iCol})*bin_width,'b.');
    [ThetaVecFit(:,iCol,i),ChTheta(i,iCol),xf,Thetayf(:,i,iCol)]=thetaLagVsDepth(ChsortII{iCol},ThetaList(i).MaxLag(ChsortII{iCol})*bin_width,ProbeXY(ChsortII{iCol},2));
    hold on;
%     plot(xf,Thetayf(:,i,iCol),'b-');
     ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
%         text(xf(20),Thetayf(20),['rmse' num2str(ThetaAcc(i,iCol))],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')

%         text(xf(20),yf(20),['rmse' num2str(ThetaVecFit(1,iCol,i))],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
    

    
% %    if iCol==1
% %       text(0,max(RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500),['Data' num2str(i)],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center');
% %    end


end

%     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);

end
% papersizePX=[0 0 3*length(ThetaList) 3*10];
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

figure;
subplot(1,2,1)
hist(RippleAcc(:),25)
subplot(1,2,2)
hist(ThetaAcc(:),25)

Acc=RippleAcc+ThetaAcc;
rTh=0.5;
tTh=0.5;
%%%%%%%%%%%
figure;
for i=1:length(ThetaList)
%         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
%         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','ticmes new roman','horizontalalignment','left','verticalalignment','top');
      ThetaList(i).MaxLag=[];
       if isempty(ThetaList(i).CA1CCAll)
          continue
       end

for iCol=1:length(ChsortII)
   subplotLU(length(ChsortII),length(ThetaList),iCol,i,P);
% %    set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4],'ylim',[0 0.006]);

   if isempty(ThetaList(i).CA1CCAll{iCol})
      continue;
   end

    plot(ProbeXY(ChsortII{iCol},2),RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500,'r.');   
    hold on;
    if RippleAcc(i,iCol)<=rTh
    plot(xf,Rippleyf(:,i,iCol),'r-','linewidth',1);
     text(xf(20),Rippleyf(20),['rmse' showNum(RippleAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
    else
     plot(xf,Rippleyf(:,i,iCol),'r-','linewidth',0.5);
     text(xf(20),Rippleyf(20),['rmse' showNum(RippleAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','k')
     
    end
    
%     text(xf(20),yf(20),['rmse' num2str(RippleVecFit(1,iCol,i))],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
    
% %    figure;
    hold on;
    
     [~,ThetaList(i).MaxLag(ChsortII{iCol})]=max(ThetaList(i).CA1CCAll{iCol});

    if std(ThetaList(i).MaxLag(ChsortII{iCol}))<0.00000000001
       continue;
        
    end

    
    plot(ProbeXY(ChsortII{iCol},2),ThetaList(i).MaxLag(ChsortII{iCol})*bin_width,'b.');
    hold on;
    if ThetaAcc(i,iCol)<tTh
    plot(xf,Thetayf(:,i,iCol),'b-','linewidth',2);
     ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
        text(xf(20),Thetayf(20),['rmse' showNum(ThetaAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
    else
        plot(xf,Thetayf(:,i,iCol),'b-','linewidth',0.5);
     ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
        text(xf(20),Thetayf(20),['rmse' showNum(ThetaAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','k')
    
    end
   if iCol==1
      text(0,max(RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500),['Data' num2str(i)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center');
   end


end

%     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);

end
papersizePX=[0 0 3*length(ThetaList) 3*10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SaveData 'ThetaRippleShankColFitAcc'],'tiff'); 


% % [~,Sh1]=min(Acc(:,1:5)');
% % [~,Sh2]=min(Acc(:,6:10)');
% % AccSh=[Sh1;Sh2+5]';

Acc=RippleAcc+ThetaAcc;
rTh=0.5;
tTh=0.5;

ValidI=ChRipple~=0&ChTheta~=0;
RTdist=zeros(size(ChRipple))+nan;
RTdist(ValidI)=ProbeXY(ChRipple(ValidI),2)-ProbeXY(ChTheta(ValidI),2);
RTdistVI=RTdist<0&RTdist>-150;   %%%Theta channel is 150 micro metter above pyrmidal channel
RTdistVI2=RTdist<0;              %%%Theta channel is above pyrmidal channel

RVI=RippleAcc<=rTh;   %%Here RippleAcc is RMSE of the gaussian fit of ripple power
TVI=ThetaAcc<=tTh;    %%Theta Acc is RMSE of Sigma fit of theta lag

ValidVI=(RVI+TVI+RTdistVI)==3;
% ValidVI2=(RVI+TVI+RTdistVI)==3;

sum(ValidVI')
% sum(ValidVI2')

% % ValidVI(sum(ValidVI')==0,:)=1;
% % ValidVI=ValidVI+RTdistVI==2;
% Invalid=find(sum(ValidVI')==0)
% Invalid=find(sum(ValidVI2')==0)

% ValidVI(Invalid,:)=RTdistVI(Invalid,:);

% % AccValid=Acc;
% % AccValid(~ValidVI)=nan;

[~,Sh1]=min(Acc(:,1:5)');
[~,Sh2]=min(Acc(:,6:10)');
% AccSh=[Sh1;Sh2+5]';

Sh1=zeros(length(ThetaList),1);
Sh2=zeros(length(ThetaList),1);
% 
% for i=1:length(ThetaList)
%     temp1=find(ValidVI(i,1:5)>0);
%     if ~isempty(temp1)
%     [~,temp2]=min(Acc(i,temp1));
%     Sh1(i)=temp1(temp2);
%     end
%     temp1=find(ValidVI(i,6:10)>0);
%     if ~isempty(temp1)
%     [~,temp2]=min(Acc(i,temp1));
%     Sh2(i)=temp1(temp2)+5;
%     end
%     
%     
%     
%     
% 
% end
% AccSh=[Sh1 Sh2];
%%%%%%%%%

for i=1:length(ThetaList)
    temp1=find(ValidVI(i,1:5)>0);
    if ~isempty(temp1)
    [~,temp2]=find(abs(RTdist(i,temp1))==min(abs(RTdist(i,temp1))));
    temp2=1:length(temp1);
    [~,temp3]=min(Acc(i,temp1(temp2)));
    Sh1(i)=temp1(temp2(temp3));
    end
    temp1=find(ValidVI(i,6:10)>0);
    if ~isempty(temp1)
    [~,temp2]=find(abs(RTdist(i,temp1+5))==min(abs(RTdist(i,temp1+5))));
        temp2=1:length(temp1);

%     [~,temp2]=min(abs(RTdist(i,temp1+5)));
    [~,temp3]=min(Acc(i,temp1(temp2)+5));
    Sh2(i)=temp1(temp2(temp3));
    end

end
AccSh=[Sh1 Sh2+5];
AccSh(AccSh(:,2)==5,2)=0;

%%%%For not satisfying the criteria, choose the best Acc channel
IV=find(sum(AccSh,2)==0);
for i=1:length(IV)
    [~,a]=min(Acc(IV(i),:));
    if a<=5
       AccSh(IV(i),1)=a; 
    else
       AccSh(IV(i),2)=a; 
       
    end
    
end


figure;
for i=1:length(ThetaList)
%         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
%         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','ticmes new roman','horizontalalignment','left','verticalalignment','top');
      ThetaList(i).MaxLag=[];
       if isempty(ThetaList(i).CA1CCAll)
          continue
       end
ThetaList(i).CA1ThChIndex=[];
ThetaList(i).CA1PyrChIndex=[];

for iCol=1:length(ChsortII)
   subplotLU(length(ChsortII),length(ThetaList),iCol,i,P);
% %    set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4],'ylim',[0 0.006]);

   if isempty(ThetaList(i).CA1CCAll{iCol})
      continue;
   end

    plot(ProbeXY(ChsortII{iCol},2),RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500,'r.');   
    hold on;
    if ismember(iCol,AccSh(i,:))
    plot(xf,Rippleyf(:,i,iCol),'r-','linewidth',1);
     text(xf(20),Rippleyf(20),['rmse' showNum(RippleAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
    plot(ProbeXY(ChRipple(i,iCol),2),RippleList(i).CA1rippleRMSAve(ChRipple(i,iCol))/500,'r*');
   
    ThetaList(i).CA1PyrChIndex=[ThetaList(i).CA1PyrChIndex ChRipple(i,iCol)];    
    
    else
     plot(xf,Rippleyf(:,i,iCol),'r-','linewidth',0.5);
%      text(xf(20),Rippleyf(20),['rmse' showNum(RippleAcc(i,iCol),4)],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','k')
     
    end
    
%     text(xf(20),yf(20),['rmse' num2str(RippleVecFit(1,iCol,i))],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
    
% %    figure;
    hold on;
    
     [~,ThetaList(i).MaxLag(ChsortII{iCol})]=max(ThetaList(i).CA1CCAll{iCol});

    if std(ThetaList(i).MaxLag(ChsortII{iCol}))<0.00000000001
       continue;
        
    end

    
    plot(ProbeXY(ChsortII{iCol},2),ThetaList(i).MaxLag(ChsortII{iCol})*bin_width,'b.');
    hold on;
    if ismember(iCol,AccSh(i,:))
    plot(xf,Thetayf(:,i,iCol),'b-','linewidth',2);
     ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
        text(xf(20),Thetayf(20),['rmse' showNum(ThetaAcc(i,iCol),4)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
   
    plot(ProbeXY(ChTheta(i,iCol),2),ThetaList(i).MaxLag(ChTheta(i,iCol))*bin_width,'b*');
    ThetaList(i).CA1ThChIndex=[ThetaList(i).CA1ThChIndex ChTheta(i,iCol)];    

   
    else
        plot(xf,Thetayf(:,i,iCol),'b-','linewidth',0.5);
%      ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
%         text(xf(20),Thetayf(20),['rmse' showNum(ThetaAcc(i,iCol),4)],...
%        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','k')
    
    end

   if iCol==1
      text(0,max(RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500),['Data' num2str(i)],...
       'fontsize',8,'fontname','times new roman','horizontalalignment','center');
   end

end

%     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);

end
papersizePX=[0 0 3*length(ThetaList) 3*10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[SaveData 'ThetaRippleShankColFitAccOpt'],'tiff'); 


save([SaveData 'ShankColAlignPyrOriCh.mat'],'ThetaList','AccSh','RippleAcc','RippleVecFit','ThetaAcc','ThetaVecFit','Acc','Rippleyf','Thetayf',...
    'ProbeXY','ChsortII','ChsortSh','ChTheta','ChRipple');



% % %%%%%%%%%%%
% % figure;
% % for i=1:length(ThetaList)
% % %         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
% % %         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','ticmes new roman','horizontalalignment','left','verticalalignment','top');
% %       ThetaList(i).MaxLag=[];
% %        if isempty(ThetaList(i).CA1CCAll)
% %           continue
% %        end
% % 
% % for iCol=1:length(ChsortII)
% %    subplotLU(length(ChsortII),length(ThetaList),iCol,i,P);
% % % %    set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4],'ylim',[0 0.006]);
% % 
% %    if isempty(ThetaList(i).CA1CCAll{iCol})
% %       continue;
% %    end
% % 
% %     plot(ProbeXY(ChsortII{iCol},2),RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500,'r.');   
% %     hold on;
% %     plot(xf,Rippleyf(:,i,iCol),'r-');
% %      text(xf(20),yf(20),['rmse' num2str(RippleAcc(i,iCol))],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
% % 
% % %     text(xf(20),yf(20),['rmse' num2str(RippleVecFit(1,iCol,i))],...
% % %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
% %     
% % % %    figure;
% %     hold on;
% %     
% %      [~,ThetaList(i).MaxLag(ChsortII{iCol})]=max(ThetaList(i).CA1CCAll{iCol});
% % 
% %     if std(ThetaList(i).MaxLag(ChsortII{iCol}))<0.00000000001
% %        continue;
% %         
% %     end
% % 
% %     
% %     plot(ProbeXY(ChsortII{iCol},2),ThetaList(i).MaxLag(ChsortII{iCol})*bin_width,'b.');
% %     hold on;
% %     plot(xf,Thetayf(:,i,iCol),'b-');
% %      ThetaAcc(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaList(i).MaxLag(ChsortII{iCol})*bin_width);
% %         text(xf(20),yf(20),['rmse' num2str(ThetaAcc(i,iCol))],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
% %    if iCol==1
% %       text(0,max(RippleList(i).CA1rippleRMSAve(ChsortII{iCol})/500),['Data' num2str(i)],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center');
% %    end
% % 
% % 
% % end
% % 
% % %     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);
% % 
% % end
% % papersizePX=[0 0 3*length(ThetaList) 3*10];
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % saveas(gcf,[SaveData 'ThetaRippleComparisonShankColFitAcc'],'tiff'); 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % 
% % load([SaveData 'ShankAlign.mat'],'ThetaAlignList');
% % ThetaListSh=ThetaAlignList;
% % clear ThetaAlignList
% % 
% %    RippleVecFitSh=zeros(5,length(ThetaList),length(ChsortSh))+nan;
% %    ThetaVecFitSh=zeros(5,length(ThetaList),length(ChsortSh))+nan;
% %    RippleAccSh=zeros(length(ThetaList),length(ChsortSh))+nan;
% %    ThetaAccSh=zeros(length(ThetaList),length(ChsortSh))+nan;
% % 
% % %%%%%%%%%%%
% % figure;
% % for i=1:length(ThetaListSh)
% % %         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
% % %         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','ticmes new roman','horizontalalignment','left','verticalalignment','top');
% % %     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);
% %     ThetaListSh(i).MaxLag=[];
% %     if isempty(ThetaListSh(i).CA1CCAll)
% %        continue
% %     end
% % 
% % for iCol=1:length(ChsortSh)
% %      subplotLU(length(ChsortSh),length(ThetaListSh),iCol,i,P);
% % % %    set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4],'ylim',[0 0.006]);
% % 
% %      if isempty(ThetaListSh(i).CA1CCAll{iCol})
% %         continue;
% %      end
% % 
% %     plot(ProbeXY(ChsortSh{iCol},2),RippleList(i).CA1rippleRMSAve(ChsortSh{iCol})/500,'r.');
% %     [RippleVecFitSh(:,iCol,i),ChRippleSh(i,iCol),xf,yf]=ripplePowerVsDepth(ChsortSh{iCol},RippleList(i).CA1rippleRMSAve(ChsortSh{iCol})/500,ProbeXY(ChsortSh{iCol},2));
% %     RippleAccSh(i,iCol)=RippleVecFit(1,iCol,i)/std(RippleList(i).CA1rippleRMSAve(ChsortSh{iCol})/500);
% %     
% %     hold on;
% %     plot(xf,yf,'r-');
% %      text(xf(20),yf(20),['rmse' num2str(RippleAccSh(i,iCol))],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
% % 
% % %     text(xf(20),yf(20),['rmse' num2str(RippleVecFit(1,iCol,i))],...
% % %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','r')
% %     
% % % %    figure;
% %     hold on;
% %     if isempty(ThetaListSh(i).CA1CCAll{iCol})
% %        continue; 
% %     end
% %     [~,ThetaListSh(i).MaxLag(ChsortSh{iCol})]=max(ThetaListSh(i).CA1CCAll{iCol});
% %     if std(ThetaListSh(i).MaxLag(ChsortSh{iCol}))<0.00000000001
% %        continue;
% %         
% %     end
% % 
% %     
% %     plot(ProbeXY(ChsortSh{iCol},2),ThetaListSh(i).MaxLag(ChsortSh{iCol})*bin_width,'b.');
% %     [ThetaVecFitSh(:,iCol,i),ChThetaSh(i,iCol),xf,yf]=thetaLagVsDepth(ChsortSh{iCol},ThetaListSh(i).MaxLag(ChsortSh{iCol})*bin_width,ProbeXY(ChsortSh{iCol},2));
% %     hold on;
% %     plot(xf,yf,'b-');
% %      ThetaAccSh(i,iCol)=ThetaVecFit(1,iCol,i)/std(ThetaListSh(i).MaxLag(ChsortSh{iCol})*bin_width);
% %         text(xf(20),yf(20),['rmse' num2str(ThetaAccSh(i,iCol))],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
% % 
% % %         text(xf(20),yf(20),['rmse' num2str(ThetaVecFit(1,iCol,i))],...
% % %        'fontsize',8,'fontname','times new roman','horizontalalignment','center','color','b')
% %     
% % 
% %     
% %    if iCol==1
% %       text(0,max(RippleList(i).CA1rippleRMSAve(ChsortSh{iCol})/500),['Data' num2str(i)],...
% %        'fontsize',8,'fontname','times new roman','horizontalalignment','center');
% %    end
% % 
% % 
% % end
% % 
% % %     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);
% % % %     set(gca,'xlim',[0 160],'visible','on','yticklabel',[]);
% % 
% % end
% % papersizePX=[0 0 3*length(ThetaListSh) 3*2];
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % 
% % saveas(gcf,[SaveData 'ThetaRippleComparisonShankFitAcc'],'tiff'); 
% % 
% % 
% % 
% % 
% % 
% % 
% % figure;
% % for i=1:length(ThetaList)
% %       ThetaList(i).MaxLag=[];
% %        if isempty(ThetaList(i).CA1CCAll)
% %           continue
% %        end
% % 
% % for iSh=1:length(ChsortSh)
% %             subplotLU(length(ChsortSh),length(ThetaList),iSh,i,P)
% %     if ~isempty(ThetaList(i).CA1CCAll{iSh})
% % 
% % %         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
% % %         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','times new roman','horizontalalignment','left','verticalalignment','top');
% % %        set(gca,'xlim',[-100 300],'visible','off','clim',[0 0.4]);
% % 
% %     plot(ProbeXY(ChsortSh{iSh},2),RippleList(i).CA1rippleRMSAve(ChsortSh{iSh})/500,'b.');
% %     hold on;
% %     [~,ThetaList(i).MaxLag(ChsortSh{iSh})]=max(ThetaList(i).CA1CCAll{iSh});
% %     plot(ProbeXY(ChsortSh{iSh},2),ThetaList(i).MaxLag(ChsortSh{iSh})*bin_width,'r.');
% % 
% %     set(gca,'xlim',[0 180],'visible','on','yticklabel',[]);
% % 
% %     end
% % end
% % 
% % %     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);
% % 
% % end
% % papersizePX=[0 0 3*length(ThetaList) 6];
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % 
% % saveas(gcf,[SaveData 'CA1ThetaRippleCompareShankRows'],'tiff'); 
% % 
% % 
% % 
% % figure;
% % InvalidFile=[];
% % for i=1:length(ThetaList)
% % 
% % % close all
% % % for i=11:20
% % % %     figure;
% % %  figure;
% %     for iSh=1:length(ChsortSh)
% % 
% %     subplotLU(length(ChsortSh),length(ThetaList),iSh,i,P)
% % %%     subplotLU(length(ChsortSh),1,iSh,1,P)
% % 
% % %         set(gca,'ylim',[-300 1],'xlim',[-0.16 0.16],'visible','off');
% % %         text(-0.15,1,[num2str(i)],'fontsize',5,'fontname','times new roman','horizontalalignment','left','verticalalignment','top');
% %     set(gca,'xlim',[-10 200],'visible','off','clim',[0 0.4]);
% %        if isempty(ThetaList(i).CA1CCAll)
% %           continue
% %        end
% % % figure;
% % 
% % %%%%%%%%%%%%phase offset fitting with sigmond function
% % xf=[-10:200];
% % ft = fittype( 'a+b/(1+exp(-c*(x-d)))', 'independent', 'x', 'dependent', 'y' );
% % opts=fitoptions('method','NonlinearLeastSquares','Robust','off','StartPoint',[min(ThetaList(i).MaxLag(ChsortSh{iSh}))-10 (max(ThetaList(i).MaxLag(ChsortSh{iSh}))-min(ThetaList(i).MaxLag(ChsortSh{iSh})))*1 1 80],...
% %     'lower',[min(ThetaList(i).MaxLag(ChsortSh{iSh}))-5 (max(ThetaList(i).MaxLag(ChsortSh{iSh}))-min(ThetaList(i).MaxLag(ChsortSh{iSh})))*0.2 -0.1 50],...
% %    'upper',[min(ThetaList(i).MaxLag(ChsortSh{iSh}))+25 (max(ThetaList(i).MaxLag(ChsortSh{iSh}))-min(ThetaList(i).MaxLag(ChsortSh{iSh})))*1.2 4 140]);
% % 
% % [fitresult, gof] = fit(ProbeXY(ChsortSh{iSh},2), ThetaList(i).MaxLag(ChsortSh{iSh})', ft, opts );
% % %%%%%%%%%%%%phase offset fitting with sigmond function
% % 
% % 
% % 
% % rmse(i)=gof.rmse;
% % a(i)=fitresult.a;
% % b(i)=fitresult.b;
% % c(i)=fitresult.c;
% % d(i)=fitresult.d;
% % yf=feval(fitresult,xf);
% % 
% % %     for iCol=1:length(ChsortII)
% % %         tempp=zscore(RippleList(i).RippleRMSBaseLine(ChsortII{iCol}))/200;
% % %         plot(ProbeXY(ChsortII{iCol},2),tempp+0.04,'-','color',ColorCol(iCol,:));
% % %         
% % % %         [~,maxI]=max(tempp);
% % % %         
% % % %         Ripplemax=ProbeXY(ChsortII{iCol}(maxI),2);
% % % %         plot([Ripplemax Ripplemax],[0 30]*bin_width,':','color',ColorCol(iCol,:));
% % % 
% % %     end
% % %     plot(ProbeXY(:,2),RippleList(i).RippleRMSBaseLine(:)/500,'r.');hold on
% %    
% % %%%%%%%%%%99 percentile point of the phase offset function
% %     ThetaChValue90(i,iSh)=a(i)+b(i)*0.8;
% % %%%%%%%%%%99 percentile point of the phase offset function
% % 
% % %%%%%%%%%%99 percentile point of the phase offset function
% %     ThetaChValue10(i,iSh)=a(i)+b(i)*0.10;
% % %%%%%%%%%%99 percentile point of the phase offset function
% % 
% %     
% %     %%%%%%%Channel closed to the 99 percentile point of phase offset function 
% % [tempValue,ChThetaI90(i,iSh)]=min(abs(ThetaList(i).MaxLag(ChsortSh{iSh})-ThetaChValue90(i,iSh)));
% % [tempValue,ChThetaI10(i,iSh)]=min(abs(ThetaList(i).MaxLag(ChsortSh{iSh})-ThetaChValue10(i,iSh)));
% % % % [tempValue,ChThetaI90(i,iSh)]=min(abs(ThetaList(i).MaxLag-ThetaChValue90(i,iSh)));
% % % % [tempValue,ChThetaI10(i,iSh)]=min(abs(ThetaList(i).MaxLag-ThetaChValue10(i,iSh)));
% %     
% % plot(xf,yf*bin_width,'b');hold on;
% % 
% % % %     for iCol=1:length(ChsortII)
% % % %         if ~isempty(intersect(I,ChsortII{iCol}))
% % % %            PyrCol(i)=iCol;
% % % %            
% % % %         tempp=zscore(RippleList(i).RippleRMSBaseLine(ChsortII{PyrCol(i)}))/200;
% % % % % %         plot(ProbeXY(ChsortII{PyrCol(i)},2),tempp+0.04,'-','color',ColorCol(PyrCol(i),:));
% % % %         
% % % %         [~,maxI]=max(tempp);
% % % %         
% % % %         Ripplemax=ProbeXY(ChsortII{PyrCol(i)}(maxI),2);
% % % %         plot([Ripplemax Ripplemax],[0 30]*bin_width,'-','color',ColorCol(PyrCol(i),:));
% % % %             break
% % % %         end
% % % %     end
% % 
% % %%%%%%%%%%%%%normlized ripple power;
% %         if max(ThetaList(i).CA1rippleRMSAve)==-1
% %        tempp=zscore(RippleList(i).CA1rippleRMSAve(ChsortSh{iSh}))/200;
% %         plot(ProbeXY(ChsortSh{iSh},2),tempp+0.04,'r.');
% % 
% %         else
% %         tempp=zscore(RippleList(i).CA1rippleRMSAve(ChsortSh{iSh}))/200;
% %         plot(ProbeXY(ChsortSh{iSh},2),tempp+0.04,'r.');
% %         end
% %         
% % %%%%%%%%%%%%%Channel with highest ripple power
% %      
% % %         TimeOffset(i)=(yfCh-ThetaChValue(i))*bin_width;
% %         Itheta=find((ProbeXY(ChsortSh{iSh},2)-ProbeXY(ChsortSh{iSh}(ChThetaI90(i,iSh)),2))>=10);
% %         if isempty(Itheta)
% %            ChThetaI= ChsortSh{iSh}(ChThetaI90(i,iSh));
% % 
% %         else
% %         [x,minI]=min(abs(ProbeXY(ChsortSh{iSh}(Itheta),2)-ProbeXY(ChsortSh{iSh}(ChThetaI90(i,iSh)),2)));
% %         ChThetaI=ChsortSh{iSh}(Itheta(minI));
% % 
% %         end
% %         
% %         
% % %         IChrippleNeed=find(ProbeXY(:,2)>=ProbeXY(ChThetaI10(i),2)&ProbeXY(:,2)<=ProbeXY(ChThetaI90(i),2));
% %         IChrippleNeed=find(ProbeXY(ChsortSh{iSh},2)>=ProbeXY(ChsortSh{iSh}(ChThetaI10(i,iSh)),2)&ProbeXY(ChsortSh{iSh},2)<=ProbeXY(ChThetaI,2));
% % %         IChrippleNeed=ChsortSh{iSh}(IChrippleNeed);
% % %         IChrippleNeed=find(ProbeXY(:,2)>=ProbeXY(ChThetaI10(i),2));
% % %         IChrippleNeed=find(ProbeXY(:,2)<ProbeXY(ChThetaI10(i),2));
% % 
% %         if isempty(IChrippleNeed)
% %            IChrippleNeed=1:length(ChsortSh{iSh});
% %         end
% %         %%%%%%%%%%%%%normlized ripple power;
% %        
% % %%%%%%%%%%%%%Channel with highest ripple power
% %         [rmax,ChRippleI]=max(tempp(IChrippleNeed));  %%%%%%%%%highest ripple power;
% %         ChRippleI=ChsortSh{iSh}(IChrippleNeed(ChRippleI));
% %         if length(find(tempp>rmax))>10
% %            InvalidFile=[InvalidFile;i];
% %         end
% % %%%%%%%%%%%%%Channel with highest ripple power
% % 
% %         RipplemaxLag=feval(fitresult,ProbeXY(ChRippleI,2));
% % 
% %         
% %         
% %        
% % % hold on;plot([xf(I) xf(I)],[0 30]*bin_width,'b:');
% % hold on;plot([ProbeXY(ChThetaI,2) ProbeXY(ChThetaI,2)],[0 30]*bin_width,'b-');
% % % % hold on;plot([ProbeXY(ChThetaI90(i,iSh),2) ProbeXY(ChThetaI90(i,iSh),2)],[0 30]*bin_width,'r:');
% % 
% % hold on;plot([ProbeXY(ChsortSh{iSh}(ChThetaI90(i,iSh)),2) ProbeXY(ChsortSh{iSh}(ChThetaI90(i,iSh)),2)],[0 30]*bin_width,'b:');
% % 
% % 
% % hold on;plot([ProbeXY(ChRippleI,2) ProbeXY(ChRippleI,2)],[0 30]*bin_width,'r-');
% % 
% % axis xy;
% % ThetaList(i).ThChIndex(iSh)=ChThetaI;
% % ThetaList(i).PyrChIndex(iSh)=ChRippleI;
% % 
% % % ThRippleOFF(i)=ProbeXY(I,2)-Ripplemax;
% % 
% % %       ThRippleOFF(i)=(ThetaList(i).MaxLag(ChTheta)-RipplemaxLag)*bin_width;
% %       ThetaRippleOFF(i,iSh)=(ThetaList(i).MaxLag(ChThetaI)-ThetaList(i).MaxLag(ChRippleI))*bin_width;
% % 
% % hold on;
% % plot((ProbeXY(ChsortSh{iSh},2)),ThetaList(i).MaxLag(ChsortSh{iSh})*bin_width,'b.');
% % text(0,26*bin_width,['rmse=' showNum(rmse(i),2)],...
% %     'fontsize',6,'fontname','times new roman','horizontalalignment','center');
% % 
% % text(0,24*bin_width,['b=' showNum(b(i),3)],...
% %     'fontsize',6,'fontname','times new roman','horizontalalignment','center');
% % 
% % text(0,28*bin_width,['Data' num2str(i)],...
% %     'fontsize',6,'fontname','times new roman','horizontalalignment','center');
% % 
% % % set(gca,'ylim',[min(ThetaList(i).MaxLag)-1 max(ThetaList(i).MaxLag)+1],'yticklabel',[])
% % set(gca,'ylim',[0 30*bin_width],'yticklabel',[],'xlim',[-20 200])
% % 
% % % barplot(-400,max(ThetaList(i).MaxLag(:)),rmse(i)*100,3,[1 0 0],1);
% % % % for iCol=1:length(ChsortII)
% % % %     plot(ProbeXY(ChsortII{iCol},2),RippleList(i).RippleRMSBaseLine(ChsortII{iCol}),'-','color',ColorCol(iCol,:));
% % % %     hold on;
% % % %     [~,ThetaList(i).MaxLag(ChsortII{iCol})]=max(ThetaList(i).CA1CCAll{iCol});
% % % %     plot(ProbeXY(ChsortII{iCol},2),ThetaList(i).MaxLag(ChsortII{iCol}),':','color',ColorCol(iCol,:));
% % % %     
% % % % end
% % % % 
% % %     set(gca,'xlim',[-300 1],'visible','off','clim',[0 0.4]);
% %     set(gca,'xlim',[0 180],'visible','on','yticklabel',[]);
% % end
% % end
% % papersizePX=[0 0 3*length(ThetaList) 6];
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % 
% % saveas(gcf,[SaveData 'CA1ThetaRippleFittingComparisonAll'],'tiff'); 
% % 
% % % figure;
% % % hist(ThRippleOFF,20)
% % % 
% % % figure;
% % % hist(TimeOffset)
% % 
% % 
% % save([SaveData 'PyrOriCh.mat']);
% % load([SaveData 'PyrOriCh.mat'])
% % 
% % 
% % 
% % % % % rippleDepth=[];
% % % % % thetaDepth=[];
% % % % % 
% % % % % for i=1:length(ThetaList)
% % % % %     rippleDepth=[rippleDepth;ProbeXY(ThetaList(i).PyrChIndex,2)];
% % % % %     thetaDepth=[thetaDepth;ProbeXY(ThetaList(i).ThChIndex,2)];
% % % % % 
% % % % % end
% % % % % 
% % % % % bth=10;
% % % % % find(b<=bth&rippleDepth<=-230)
% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%
% % % % % P.yBottom=0.1;
% % % % % P.xLeft=0.1;
% % % % % P.xInt=0.04;
% % % % % 
% % % % % figure;
% % % % % DiffDepth=thetaDepth-rippleDepth;
% % % % % Depth=[0 200];
% % % % % subplotLU(1,2,1,1,P)
% % % % % hist(DiffDepth,20);hold on;
% % % % % Invalid0=find(DiffDepth<=Depth(1)|DiffDepth>=Depth(2))  %%%%%%%Depth
% % % % % plot([Depth(1) Depth(1)],[0 15],'k:');
% % % % % plot([Depth(2) Depth(2)],[0 15],'k:');
% % % % % % barplot(DiffDepth(1),)
% % % % % barplotLu(Depth(1),0,diff(Depth),15,[0 1 0],0.2);
% % % % % 
% % % % % 
% % % % % xlabel('Depth difference (ThetaCh-RippleCh)')
% % % % % ylabel('Counts');
% % % % % set(gca,'ylim',[0 15],'ytick',0:5:15)
% % % % % subplotLU(1,2,1,2,P)
% % % % % hist(rmse,50);hold on;
% % % % % rmseth=4;    %%%%%%%%%fitting quality (root mean square error);
% % % % % Invalid1=find(rmse>=rmseth);
% % % % % plot([rmseth rmseth],[0 15],'k:');
% % % % % xlabel('rmse')
% % % % % set(gca,'ylim',[0 15],'ytick',0:5:15,'yticklabel',[])
% % % % % barplotLu(0,0,rmseth,15,[0 1 0],0.2);
% % % % % 
% % % % % papersizePX=[0 0 12 6];
% % % % % set(gcf, 'PaperUnits', 'centimeters');
% % % % % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % % % % saveas(gcf,[SaveData  'IncludedData.tiff'],'tiff'); 
% % % % % saveas(gcf,[SaveData  'IncludedData.eps'],'epsc'); 
% % % % % 
% % % % % InvalidFile=union(InvalidFile,Invalid0);
% % % % % InvalidFile=union(InvalidFile,Invalid1);
% % % % % InvalidFile=union(InvalidFile,9);
% % % % % save([SaveData,'BadRippleThetaChFile.mat'],'InvalidFile','ThetaList','a','b','c','d','rmse','Depth','thetaDepth','rippleDepth','rmseth','Depth');
% % % % % 

