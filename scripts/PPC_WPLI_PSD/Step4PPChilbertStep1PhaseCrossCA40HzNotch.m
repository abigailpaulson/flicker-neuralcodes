clear all
% load('\\ad.gatech.edu\bme\labs\singer\LuZhang\Project1Aging\Results\step1\rippleRMSInfo75aboveRippleCh.mat');

LoadData='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step1\';
load([LoadData 'StimInfoList.mat'],'StimFileList');
%%%%%%%Channels chosen by best Col among 5 Cols of one Shank
% load([LoadData 'DetermineCA1BoardPyrRunList.mat'],'RunList');
load('Y:\singer\LuZhang\Project2ChronicF\Results\step1\ThetaAlign\ShankColAlignCA1PyrOriCA3PyrCh.mat','ThetaList');
RunList=ThetaList; clear ThetaList;
%%%%%%%Channels chosen by best Col among 5 Cols of one Shank

% % DataStored='Y:\singer\ProcessedData\VR_AnnularTrack_Lu\';
% % ThetaDataStored='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project1Aging\Results\step1\';

LoadSpiking='\\ad.gatech.edu\bme\labs\singer\Abby\chronicflicker_annulartrack\decoding\bayesian\data\';
% % load([LoadSpiking 'NeuronData.mat'],'FileList','NeuronFile');       
load([LoadSpiking 'TrialInfo.mat']);       
load([LoadSpiking 'NeuronData_all.mat']);       

%%%Checking Abby extracting file list match to Lu Zhang file list
for i=1:length(FileList)
    if (FileList(i).Subj~=RunList(i).Subj)||(FileList(i).Session~=RunList(i).Session)
       disp([num2str(i) 'not match']);
    end
end
%%%Checking Abby extracting file list match to Lu Zhang file list


SavePath='\\ad.gatech.edu\bme\labs\singer\LuZhang\Project2ChronicF\Results\step4\';
mkdir(SavePath);
% % SavePath=[SavePath '\NoSpeedThFourierPSDBoardPyr\'];
% mkdir(SavePath);


% mkdir(SavePath)
ChList=0:63;   %%%%ChID for data storing
% 
% 
samprate=2000;

psdParameter.Fs=samprate;
psdParameter.window=2048;
psdParameter.noverlap=0;
psdParameter.nfft=4096;
nFre=psdParameter.nfft/2+1;
% 

F0 = 40;   % Interference is at 40 Hz
Fs = 2000; % Sampling frequency is 3000 Hz
BW = 1;    % Choose a bandwidth factor of 6Hz
[num1,den1] = iirnotch(F0/(Fs/2),BW/(Fs/2));
fvtool(num1,den1,'Fs',Fs);
set(gca,'xlim',[0 100])
close all
% y = filter(num1,den1,data);


   
RunTh=10;
RewardTh=5;

 TrialName={'Correct' 'Wrong'};
%  for itest=1:size(CompariGroup,2) 
paramM.lowSpeedThreshold=0; %%%%%5degree/s
paramM.highSpeedThreshold=360.00;%%%%%50cm/s
paramM.binWidth=10;
paramM.sampleTime=median(diff(TrialInfo(1).Trial(1).PosData.timeSmooth));
paramM.minBinTime=0;

SpaBinNum=round(180/paramM.binWidth);



ComparisonGroup=3;

clear passband
% Gamma=[20:40:140];
% 
% passband(1,:)=[2 6 10 Gamma(1:end-1)];   
% passband(2,:)=[6 10 20 Gamma(2:end)]; 
%%%Current Use for boarder gamma bands, simplified the results.
Gamma=[25 60 100;55 100 140];
passband(1,:)=[2 6 10 Gamma(1,:)];   
passband(2,:)=[6 10 20 Gamma(2,:)]; 
%%%Current Use for boarder gamma bands, simplified the results.


phaseMethod='hilbert';
powerThresh=[];
DownR=4;
SampleR=2000/DownR;
MorseParam=[3 20 4];  %%%%[gamma beta D]



%%%%%%Checking
%%%%%%%using morsespace(gamma,beta,D) to determin frequency points within
%%%%%%%each band
for j=1:size(passband,2)
% %         nvoice = 12;
% %         freqlist= 2.^(log2(passband(1,i)):1/nvoice:log2(passband(2,i)));
% %         length(freqlist)
        
        gamma=MorseParam(1);
        beta=MorseParam(2);
        D=MorseParam(3);
        dt=1/SampleR;
        HIGH=passband(2,j)*dt*2*pi;
        LOW=passband(1,j)*dt*2*pi;   %%%%%%frequencies in radius
        Frad=morsespace(gamma,beta,HIGH,LOW,D);
        passbandNum(j)=length(Frad);
        Frad
end
%%%%%%%using morsespace(gamma,beta,D) to determin frequency points within
% %%%%%%%each band
% %%%%%%Checking


NeuronFileS1=NeuronFile;
NeuronFileS2=NeuronFile;

    

%%%%%%%%%%%%%%consider all trial
for i=1:length(FileList)
        tic
        
        clear temp2CA1 temp2CA3
           temp1CA1=[FileList(i).PathFile 'CA1\'];
           ChPyrCA1=ChList(RunList(i).CA1PyrChIndex);
           for iCh=1:length(ChPyrCA1)
           temp2CA1{iCh}=[temp1CA1 num2str(ChPyrCA1(iCh)) '\'];    %%%%%%%%CA1 Shank iCh
           end

           temp1CA3=[FileList(i).PathFile 'CA3\'];
           ChPyrCA3=ChList(RunList(i).CA3PyrChIndex);
           for iCh=1:length(ChPyrCA3)
           temp2CA3{iCh}=[temp1CA3 num2str(ChPyrCA3(iCh)) '\'];    %%%%%%%%CA3 Shank iCh
           end

           
   temp3CA1=[temp1CA1 'sorted\'];
   temp3CA3=[temp1CA3 'sorted\'];

   if     isempty(dir(temp3CA1))&&isempty(dir(temp3CA3))
          continue
   end
           
%         RippleRMS(iCH);
%         tempP=[];
           tempPower=[];
           maxChCA1=[];
           maxShCA1=[];
           maxSh2CA1=[];

           if ~isempty(find(NeuronFile(i).brainReg==1))
           load([temp3CA1 'clustermetrics.mat']);
           for iCell=1:length(clustermetrics)
               maxChCA1(iCell)=clustermetrics(iCell).maxChan;
           end
           maxShCA1=[maxChCA1>=32]+1;   %%%%%%%Shank 1 for CH0-31; Shank 2 for Ch32-63;
           IndSh2CA1=find(maxShCA1==2);   %%%%%%%Cell Index for Shank2
           end
           
           maxChCA3=[];
           maxShCA3=[];
           maxSh2CA3=[];
           if ~isempty(find(NeuronFile(i).brainReg==2))
           load([temp3CA3 'clustermetrics.mat']);
           for iCell=1:length(clustermetrics)
               maxChCA3(iCell)=clustermetrics(iCell).maxChan;
           end
           maxShCA3=[maxChCA3>=32]+1;   %%%%%%%Shank 1 for CH0-31; Shank 2 for Ch32-63;
           IndSh2CA3=find(maxShCA3==2);   %%%%%%%Cell Index for Shank2
           end
           
           maxSh=[maxShCA1(:);maxShCA3(:)];
           IndSh2=find(maxSh==2);
           
           
           
           
               for iFile=1:length(FileList(i).File)
%                             CountF=0;    %%%%%%%whether count the file or not?0 not count yet;1 counted already
%                    SS1=SS(sampleFID==iFile);
%                    TS1=sampleTS(iFile).ThetaTS;
                 clear LFPTempCA1
                 for iSh=1:2
                     if iSh>length(temp2CA1)
                        if length(dir(temp2CA1{iSh-1}))>=1
                           load([temp2CA1{iSh-1} 'eeg' num2str(FileList(i).File(iFile))]);
                           y = FiltFiltM(num1,den1,eeg{end}{end}{end}.data(:));
                           LFPTempCA1{iSh}=decimate(y,DownR);

%                            LFPTempCA1{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);
                           continue
                        end

                     end

                     if length(dir(temp2CA1{iSh}))>=1
                         load([temp2CA1{iSh} 'eeg' num2str(FileList(i).File(iFile))]);
                         y = FiltFiltM(num1,den1,eeg{end}{end}{end}.data(:));
                         LFPTempCA1{iSh}=decimate(y,DownR);

%                          LFPTempCA1{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);
                     end
                 end
                 
                 clear LFPTempCA3
                 for iSh=1:2
                     if iSh>length(temp2CA3)
                        if length(dir(temp2CA3{iSh-1}))>=1
                         load([temp2CA3{iSh-1} 'eeg' num2str(FileList(i).File(iFile))]);
                         y = FiltFiltM(num1,den1,eeg{end}{end}{end}.data(:));
                         LFPTempCA3{iSh}=decimate(y,DownR);

%                          LFPTempCA3{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);
                         continue
                        end
                         
                     end

                     if length(dir(temp2CA3{iSh}))>=1
                         load([temp2CA3{iSh} 'eeg' num2str(FileList(i).File(iFile))]);
                         y = FiltFiltM(num1,den1,eeg{end}{end}{end}.data(:));
                         LFPTempCA3{iSh}=decimate(y,DownR);

%                          LFPTempCA3{iSh}=decimate(eeg{end}{end}{end}.data(:),DownR);
                     end
                 end
               
                 CA1I=find(NeuronFile(i).brainReg==1);
                 CA1UnitID=NeuronFile(i).TSIDtotal(CA1I);
                 if length(CA1I)~=length(maxShCA1)
                    disp(['Non matching CA1 Unit in ' num2str(i) ' SubjFile']); 
                 end
                 
                 
                 CA3I=find(NeuronFile(i).brainReg==2);
                 CA3UnitID=NeuronFile(i).TSIDtotal(CA3I);
                 if length(CA3I)~=length(maxShCA3)
                    disp(['Non matching CA3 Unit in ' num2str(i) ' SubjFile']); 
                 end

                 
                 CAUnitID{1}=CA1UnitID;
                 CAUnitID{2}=CA3UnitID;
                
                 LFPTempCA{1}=LFPTempCA1;
                 LFPTempCA{2}=LFPTempCA3;
                 
                           if length(CA3I)~=length(maxShCA3)
                             disp(['Non matching CA3 Unit in ' num2str(i) ' SubjFile']); 
                          end
                          if length(CA1I)~=length(maxShCA1)
                             disp(['Non matching CA1 Unit in ' num2str(i) ' SubjFile']); 
                          end

                 for iSh=1:2
                     TempPhase1{iSh}=zeros(length(NeuronFile(i).Data(iFile).TS),size(passband,2))+nan;
                     for iRegion=1:2
                                          
                         [localCAI,~]=ismember(NeuronFile(i).Data(iFile).TSID,CAUnitID{iRegion});
                         jRegion=setdiff([1 2],iRegion);        %%%%%%%%Cross region analysis, using LFP from CA1/CA3 for CA3/CA1 units
                     if ~isempty(LFPTempCA{jRegion}{iSh})
                        tempA= GT_Spike2Phase(NeuronFile(i).Data(iFile).TS(localCAI),LFPTempCA{jRegion}{iSh},SampleR,passband,phaseMethod,powerThresh,MorseParam);
                        TempPhase1{iSh}(localCAI,:)=tempA;
                     end
                     end
                 end
          
                 
                  for iSh=1:2
                     if isempty(TempPhase1{iSh})
                        TempPhase1{iSh}=TempPhase1{setdiff([1 2],iSh)};
                     end
                 end             
% %                  [~,I1,I2]=intersect(NeuronFile(i).Data(iFile).TSID,IndSh2);
                 NeuronFileS1(i).Data(iFile).Phase=TempPhase1{1};
                 NeuronFileS2(i).Data(iFile).Phase=TempPhase1{2};



% %                  if isempty(I1)
% %                  NeuronFile(i).Data(iFile).Phase(I1)= TempPhase1{2}(I1);
% %                  end
                     

              
              end
%                end
                % %                clear TempTsNumOcc TempTimeOcc TempTsNumTran TempTimeTran
          toc   
end

   NeuronDataSh{1}=NeuronFileS1;
   NeuronDataSh{2}=NeuronFileS2;

   save([SavePath num2str(size(passband,2)) 'bandPPCHilbertCrossCAShankColAlignphase.mat'],'NeuronDataSh','FileList','passband','MorseParam','-v7.3');       
%    
%    
%    load([SavePath num2str(size(passband,2)) 'bandPPCHilbertphase.mat']);       

% % NeuronDataSh{1}=NeuronDataS1;
% % NeuronDataSh{2}=NeuronDataS2;

clear NeuronData;

for jSh=1:2
  
    NeuronData=NeuronDataSh{jSh};
clear PhaseTrial


for i=1:length(StimFileList)
%  for i=37:37

    if isempty(NeuronData(i).TSIDtotal)
       continue; 
    end
   
    tic
% %    end
       if ismember(1,StimFileList(i).Stim) 
          FlickerFileID=find(StimFileList(i).RType==2);
          PreFlickerFileID=find(StimFileList(i).RType==1);
          PostFlickerFileID=find(StimFileList(i).RType==3);

%        continue
          
% %           PreFlickerFileID=[1:(min(FlickerFileID)-1)];
% %           PostFlickerFileID=[(max(FlickerFileID)+1):length(StimFileList(i).Stim)];
       elseif (~ismember(1,StimFileList(i).Stim))&&StimFileList(i).RType(1)==0
          PreFlickerFileID=find(StimFileList(i).RType==0);
          PostFlickerFileID=[];


       else
          PreFlickerFileID=[1:length(StimFileList(i).Trial)];
          PostFlickerFileID=[];
          FlickerFileID=[];
%           continue
       end
       
       PrePostFileID={PreFlickerFileID PostFlickerFileID};


       for iPrePost=1:length(PrePostFileID)
              if isempty(PrePostFileID{iPrePost})
                 continue;
              end
              
              
              

                        for iCell=1:length(NeuronData(i).TSIDtotal)
                        UnitPhase{iCell}(1).Phase=[];
                        UnitPhase{iCell}(1)=[];
                        end

  
              
              TrialTypeTotal=[];

           for jFile=1:length(PrePostFileID{iPrePost})
               iFile=PrePostFileID{iPrePost}(jFile);
               if iFile>length(StimFileList(i).Trial)
                  continue;
               end
               
               
               
               
               
                    CountF=0;    %%%%%%%whether count the file or not?0 not count yet;1 counted already
                
                    clear TrialPeriodtemp TrialRewardtemp TrialZonetemp
                    TrialPeriodtemp=[StimFileList(i).Trial(iFile).Zone1TrialTs StimFileList(i).Trial(iFile).Zone2TrialTs];
                    TrialRewardtemp=[StimFileList(i).Trial(iFile).Zone1Type StimFileList(i).Trial(iFile).Zone2Type];
                    TrialZonetemp=[zeros(size(StimFileList(i).Trial(iFile).Zone1Type))+1 zeros(size(StimFileList(i).Trial(iFile).Zone2Type))+2];
                    
                    
                    if isempty(TrialPeriodtemp)
                    continue
                    end

                    [~,s1I]=sort(TrialPeriodtemp(1,:));
                    TrialPeriodtemp=TrialPeriodtemp(:,s1I);
                    TrialRewardtemp=TrialRewardtemp(s1I);
                    TrialZonetemp=TrialZonetemp(s1I);

                    TrialTypeTotal=[TrialTypeTotal;TrialRewardtemp(:)];

                      
                    for iTrial=1:length(TrialRewardtemp)
                        [~,I1]=TsInTimerange(NeuronData(i).Data(iFile).TS,TrialPeriodtemp(:,iTrial));
                        
                        if ~isempty(I1)
                           I1ID=NeuronData(i).Data(iFile).TSID(I1);
                           I1Phase=NeuronData(i).Data(iFile).Phase(I1,:);

                           I1IDall=unique(I1ID);
                        for iCell=1:length(NeuronData(i).TSIDtotal)
                            I2=intersect(I1IDall,NeuronData(i).TSIDtotal(iCell));
                             UnitPhase{iCell}(end+1).Phase=I1Phase(I1ID==NeuronData(i).TSIDtotal(iCell),:);
                        end
                        else
                            for iCell=1:length(NeuronData(i).TSIDtotal)
                                 UnitPhase{iCell}(end+1).Phase=[];
                            end
                        end
                    end
                    


               
           end
       PhaseTrial(i,iPrePost).UnitPhase=UnitPhase;
       PhaseTrial(i,iPrePost).TrialReward=TrialTypeTotal;
                clear UnitPhase;    

       end


    toc
end

PhaseVRSh{jSh}=PhaseTrial;


clear PhaseTrial
FlickerTrialT=60;   %%%%%define the trial duration as 60 seconds for flicker period 

for i=1:length(StimFileList)
%  for i=37:37
       if isempty(NeuronData(i).TSIDtotal)
       continue; 
    end
   

    tic

               FileCount=0;
% %                PSDFile(i).PSD=zeros(nFre,numBin,3);
% %    end
       if ismember(1,StimFileList(i).Stim) 
          FlickerFileID=find(StimFileList(i).RType==2);
          FlickerSalineFileID=find(StimFileList(i).RType==4);

       else
          FlickerFileID=[];
          FlickerSalineFileID=[];
       end
       
       PrePostFileID={FlickerFileID FlickerSalineFileID};


       for iPrePost=1:length(PrePostFileID)
              
              if isempty(PrePostFileID{iPrePost})
                 continue;
              end
              
              
             for iCell=1:length(NeuronData(i).TSIDtotal)
                 UnitPhase{iCell}(1).Phase=[];
                 UnitPhase{iCell}(1)=[];
             end

              
           for jFile=1:length(PrePostFileID{iPrePost})
               iFile=PrePostFileID{iPrePost}(jFile);
                CountF=0;    %%%%%%%whether count the file or not?0 not count yet;1 counted already

                PeriodTask=StimFileList(i).FlickerInfo(iFile).Period;

                if isempty(PeriodTask)
                   continue
                end
                  
                NTrial=floor(diff(PeriodTask)/FlickerTrialT);                     
                
                clear PeriodTrial
                PeriodTrial(1,:)=([1:NTrial]-1)*FlickerTrialT;
                PeriodTrial(2,:)=[1:NTrial]*FlickerTrialT;
       
                PeriodTrial=PeriodTrial+PeriodTask(1);

               
                
                                      
                for iTrial=1:size(PeriodTrial,2)
                        [~,I1]=TsInTimerange(NeuronData(i).Data(iFile).TS,PeriodTrial(:,iTrial));
                        
                        if ~isempty(I1)
                           I1ID=NeuronData(i).Data(iFile).TSID(I1);
                           I1Phase=NeuronData(i).Data(iFile).Phase(I1,:);

                           I1IDall=unique(I1ID);
                        for iCell=1:length(NeuronData(i).TSIDtotal)
                            I2=intersect(I1IDall,NeuronData(i).TSIDtotal(iCell));
                             UnitPhase{iCell}(end+1).Phase=I1Phase(I1ID==NeuronData(i).TSIDtotal(iCell),:);
                        end
                        else
                            for iCell=1:length(NeuronData(i).TSIDtotal)
                                 UnitPhase{iCell}(end+1).Phase=[];
                            end
                        end
                    end

                
               
           end
           
             PhaseTrial(i,iPrePost).UnitPhase=UnitPhase;
             clear UnitPhase;    


       end
       

    toc
end

PhaseFlickerSh{jSh}=PhaseTrial;
end

save([SavePath num2str(size(passband,2)) 'bandPPCHilbertCrossCAShankColAlignphaseVRFlicker40Notch.mat'],'PhaseVRSh','PhaseFlickerSh','NeuronDataSh','FileList','passband','phaseMethod','MorseParam','-v7.3');       
    
clear all