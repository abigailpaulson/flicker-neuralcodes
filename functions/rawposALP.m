function rawposALP(index, vrDataFolder, intanFolder)
    %% Raw Position Data Structure Construction  7/27/2017 - MKA
    % updated 11.13.18 to correct for misalignments of VR and intan data
    % Extract data (positon, # of licks, rewards) from VR files and 
    % creates a data structure from it.
    % Function to be used in preprocess_Intan script.
    % from rawposSMP - abby 12/18/19
    %updated ALP 2/17/22 to save manual reward information
  
    
    %% ------- Directory --------------------------------------------- 
%     vrDataFolder = 'Z:\singer\SyncingToIntan\testData\VRtestdata\';
%     savedatadir= 'Z:\singer\SyncingToIntan\rawPosStructs\';
     
    %% ---------- Import VR Data and Array of Synced Triggers ------------
    disp(['extracting raw position index: ', num2str(index)])
    %import SaveDatCopy from VR file
    if exist([vrDataFolder, '\dataWithLickometer.mat'])
        load([vrDataFolder, '\dataWithLickometer']);
    elseif exist([vrDataFolder, '\virmenDataRaw.mat']) %after 9/11/20
        load([vrDataFolder, '\virmenDataRaw.mat']) 
    else
        disp('dataWithLickometer missing')
        return
    end
    
    %gets data matrix from struct, if saveDatCopy was saved as struct
    if ~exist('saveDatCopy', 'var') %added ALP 9/25/20 for updated saving structure 
        saveDatCopy = virmenData.data; 
    end
    
    if isa(saveDatCopy,'struct'); %#ok<NODEF>
        saveDatCopy = saveDatCopy.data;
    end
        
    vrTimes(:,1) = saveDatCopy(:,1); %#ok<NODEF>   %VRvalues 
    vrTimes(:,2) = (vrTimes(:,1) - vrTimes(1,1)) * 24 * 60 * 60; 
        % ^Converts VR values to seconds elapsed
    syncedArray = triggerSyncALP(index, vrDataFolder, intanFolder); % returns array of synced triggers 
    trigList = (syncedArray(:,3) - syncedArray(1,3)) * 24 * 60 * 60;
    % A column vector of all VR trigger values
    % starts from 0 where 0 is when the isrecording trigger was clicked and
   
    
    %Do the same as above to get VR Indices
    vrIndices(:,1) = saveDatCopy(:,1); %#ok<NODEF>   %VRvalues 
    vrIndices(:,2) = (vrTimes(:,1) - vrTimes(1,1)); 
    % ^Converts VR values to indices elapsed
    trigListIndices = (syncedArray(:,3) - syncedArray(1,3));
    % A column vector of all VR trigger values
    
    
    %% ------------- Convert all VR times to Intan ------------------------
    % Finds range of vrTime when the isRecording trigger is on
    
    i = 1;
    
    %Skip thru the vrTimes to the first VR value after the first auto trigger
    %need to compare raw times instead of elapsed bc vrTimes and
    %syncedArray start at different time points
    while vrTimes(i,1) < syncedArray(2,3) 
        i = i+1;
    end
    
    iStart = i; %  VR index when intan starts (first isRecording and autoTrigger trigger point)
    %this iStart time is after the first AUTOtrigger, whileas the first
    %value in trigList would be when the firs
    for t = 2:length(trigList)-1 
        % Find best-fit line for approximation
        % trig list is the VR times from isRecording trigger (time elapsed)
        % synced array column 2 is the Intan index from isRecording trigger
        n = polyfit(trigList(t:t+1),syncedArray(t:t+1,2),1); 
        while vrTimes(i,1) <= syncedArray(t+1,3)
            % Convert VR values to Intan 
            % vrTimes here is in seconds elapsed since virmen was started
            % but it starts at the iStart value so it actually is starting
            % at the first autoTrigger
            % newIntan now starts from the first autoTrigger trigger and is converted
            % from virmen time to intan time via the synced array
            newIntan(i-iStart+1,1) = n(1) * (vrTimes(i,2)-vrTimes(iStart,2)) + n(2); %#ok<AGROW>  
            i = i+1;
        end
    end
    
    iEnd = i-1; % VR index when intan ends (final trigger point)
    
    %this section results in the newIntan vector, which now should match the
    %length of the rest of the raw pos data and should be the intan time
    %(ie the time corresponding to all of the ephys data files) when all of this data occurred 
    
     %% ------------ Convert all VR indices to Intan indices --------------
%     
%   
%      i = 1;
%     
%     %Skip thru the vrTimes to the first VR value after the first trigger
%     while vrIndices(i,2) < trigListIndices(2) 
%         i = i+1;
%     end
%     
%     iStart = i; %  VR index when intan starts (first trigger point)
%     
%     for t = 2:length(trigListIndices)-1 
%         % Find best-fit line for approximation
%         n = polyfit(trigListIndices(t:t+1),syncedArray(t:t+1,2),1); 
%         while vrIndices(i,2) <= trigListIndices(t+1)
%             % Convert VR values to Intan 
%             newIntanIndices(i-iStart+1,1) = n(1) * vrIndices(i,2) + n(2); %#ok<AGROW>    
%             i = i+1;
%         end
%     end
%     
%     iEndInd = i-1; % VR index when intan ends (final trigger point)
    
    
    

    %% ------------ Create & Save Raw Position Data Structure -------------
    %two different sizes because when I updated the way that the virmen
    %info is saved for in progress saving during the session 
    %ALP 2/17/2022 to extract the manual reward vector as well
 if size(saveDatCopy,2) == 11
    rawpos{index(1)}{index(2)}{index(3)} = struct(...
        'vrtimes',vrTimes(iStart:iEnd,1),...
        'vrind', vrIndices(iStart:iEnd,1),...
        'indices', newIntan,...%NOTE THIS OUTPUT IS IDX NOT TIME - newIntan is now the corresponding intan time for every VR time
        'index',index,...
        'theta',saveDatCopy(iStart:iEnd,2),...
        'licks',saveDatCopy(iStart:iEnd,6),...
        'reward',saveDatCopy(iStart:iEnd,5),...
        'manualreward', saveDatCopy(iStart:iEnd,7), ...
        'descrip', 'Raw position data relative to Intan time or Intan indices.');
 elseif size(saveDatCopy,2) == 13
     rawpos{index(1)}{index(2)}{index(3)} = struct(...
        'vrtimes',vrTimes(iStart:iEnd,1),...
        'vrind', vrIndices(iStart:iEnd,1),...
        'indices', newIntan,...%NOTE THIS OUTPUT IS IDX NOT TIME - newIntan is now the corresponding intan time for every VR time
        'index',index,...
        'theta',saveDatCopy(iStart:iEnd,2),...
        'licks',saveDatCopy(iStart:iEnd,6),...
        'reward',saveDatCopy(iStart:iEnd,5),...
        'manualreward', saveDatCopy(iStart:iEnd,8), ...
        'descrip', 'Raw position data relative to Intan time or Intan indices.');
 end
    
    %this is a sanity check for anyone who ever needs one SP 11.13.18 
%     figure; plot(trigList,syncedArray(:,2),'o'); hold on;
%     plot((rawpos{index(1)}{index(2)}{index(3)}.vrtimes-rawpos{index(1)}{index(2)}{index(3)}.vrtimes(1))*60*60*24,rawpos{index(1)}{index(2)}{index(3)}.indices)
%     
    save([intanFolder,'rawpos',num2str(index(3))],'rawpos');
    disp('Raw Position Data Structure Saved.');
    
    % keyboard; %debugging purposes (remove line in final version)
    
end