function [syncedArray] = triggerSyncALP(index, vrFolder, intanFolder)
%% triggerSync: syncs Intan trigger times to VR trigger times - MKA
%   MUST RUN FUNCTION AS ADMINISTRATOR in MATLAB 2015
%   Outputs an array listing the synced time/samples that each file records 
%   the triggers turning OFF or ON
%   Both intanFile and saveDatCopy should be either loaded in the workspace or
%   in path.
%   intanFile may need to be extracted using read_Intan_RHD2000_file.m 
%   The name of the intan file may look like "board_dig_in_data" or tempdata.
%   The name of the matlab file may look like "saveDatCopy".
%   Updated for abby's track ALP 12/18/19

%% Directory 
% Double-check the Drive Letters e.g. Z:, Y: etc.
    
%     intanFolder = 'Z:\singer\SyncingToIntan\testData\extractedIntan\';
%     vrFolder = 'Z:\singer\SyncingToIntan\testData\VRtestdata\';
    
%% Import data arrays from intan and VR files
    load([intanFolder,'trigger',num2str(index(3)),'.mat'])
    
%     intanData = trigger{index(1)}{index(2)}{index(3)}.data;
%     
%     while isa(intanData{1},'cell')
%         intanData = intanData{1};
%     end
    
    isRecording = trigger{index(1)}{index(2)}{index(3)}.recordingtrigger; %marked by the start and end point triggers
    intanAutotrigs = trigger{index(1)}{index(2)}{index(3)}.autotrigger;
    
    if exist([vrFolder, '\dataWithLickometer.mat'])
        vrStruc = load([vrFolder, '\dataWithLickometer.mat']);
    else
       vrStruc = load([vrFolder, '\virmenDataRaw.mat']); 
    end
    
    vrCell = struct2cell(vrStruc);
    while isa(vrCell{1},'cell') 
        vrCell = vrCell{1};
    end
    vrData = vrCell{1};
    
    % get array from struct. Not all tracks save vrData as an array
    if isa(vrData,'struct');
        vrData = vrData.data;    
    end
    
    %Last 4 of vrData are the vr.IntanTrigger measurments
    %vrData(end-3) - isRecording
    %vrData(end-2) = ?
    %vrData(end-1) = autoTrigger
    %vrData(end) = Reward pump
    
    isVRrecording = vrData(:,end-3); %updated from 10 for Abby recordings 
    vrAutotrigs = vrData(:,end-1); %updated from 12 for abby recordings
    
    %% Get Triggers From INTAN File
    p = -1; % previous Intan index value
    sai = 1; %synced array index
    
    %traverses thru intan file's trigger data (first row)
    for i = 1:size(intanAutotrigs,2)
        if isRecording(i) == 1
            if intanAutotrigs(i) ~= p
                %checks to see if trigger value has changed since last
                %iteration
                syncedArray(sai,1) = intanAutotrigs(i); %saves trigger value
                syncedArray(sai,2) = i; %saves sample number to new array
                sai = sai + 1; % increments sync array index
                p = intanAutotrigs(i); %changes previous index value
            end
        end
    end

%% me testing trigger sync stuff    
% isRecordingIdx = find(isRecording);
% intanTrigsDuringRec = intanAutotrigs(isRecordingIdx);
% trigOnOff = find(diff(intanTrigsDuringRec) ~= 0);
% testsync = intanTrigsDuringRec(trigOnOff);
% lengthOnOff = diff(trigOnOff);
% 
% isVRRecordingIdx = find(isVRrecording);
% VRTrigsDuringRec = vrAutotrigs(isVRRecordingIdx);
% VRtrigOnOff = find(diff(VRTrigsDuringRec) ~= 0);
% testsyncVR = VRTrigsDuringRec(VRtrigOnOff);
% lengthOnOffVR = diff(VRtrigOnOff);
%% Get Triggers from VR File
    
    pp = -1; %previous VR index value
    saj = 1; %syncArray index
    
    for j = 1:size(vrData,1) % j = current VR index
        if isVRrecording(j) == 1  
            if vrAutotrigs(j) ~= pp
                syncedArray(saj,3) = vrData(j,1);
                triggerCheck(saj,1) = vrAutotrigs(j); %#ok<AGROW>
                saj = saj + 1;
                pp = vrAutotrigs(j);
            end
        end
    end    
    
%% 
    %checks both VR and intan data for matching trigger statuses
    % they should both be the same unless recording was not properly set-up
    
    if ~isequal(syncedArray(:,1), triggerCheck(:,1));
        disp('Error: Trigger Mismatch - Recording Problem Detected')
        if size(syncedArray(:,1),1)-size(triggerCheck(:,1),1) >= 1
            earlyEnd = size(triggerCheck(:,1),1);
            syncedArray = (syncedArray(1:earlyEnd-1,:));
        else
            triggerStopped = find(diff(syncedArray(:,1)) == 0);
            earlyEnd = triggerStopped(1);
            syncedArray = (syncedArray(1:earlyEnd-1,:));
        end
        disp('Trigger Rematched - check final results')
    end
        
    disp(['Intan & VR successfully synced for Test',num2str(index)]); 
    
end
    

    

