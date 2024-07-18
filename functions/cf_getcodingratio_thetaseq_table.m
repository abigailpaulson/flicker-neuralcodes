function out = cf_getcodingratio_thetaseq_table(decodedData, posvals, params)
%cf_getcodingratio_thetaseq()
%
%ALP 1/18/2023

qX = params.qX;
qY = params.qY; 
dec_edges = params.dec_edges; %bins for calculating the PCR over the track

%%% get per trial PCR
tmpNTrials = max(unique(cell2mat(decodedData.trialID))); 
allDecodedSeq = cat(3, decodedData.trialSeq{:});
allTroughPos = cell2mat(decodedData.troughPos);
allTrialNum = cell2mat(decodedData.trialID);

PCR_pos_trial = []; PCR_loc_trial = []; QR_pos_trial = []; mn_seq_reg_trial = [];
for t = 1:tmpNTrials
    isT = allTrialNum == t;
    tSeq = allDecodedSeq(:,:,isT);
    tPos = allTroughPos(isT);
    
    %get the correct indices for the RRZ, RZ, AZ, etc. based on the largest
    %position value
    if  min(tPos) > 150 %360 trials with position from 180 to 360
        posvals{1} = [216 244];
        posvals{2} = [244+9 244+9+54];
        posvals{3} = posvals{1}+90;
    elseif min(tPos) < 0 %distance 2 reward trials
        posvals{1} = [-18 9]; %RRZ
        posvals{2} = [18 72]; %middle of the track speed split area
        posvals{3} = posvals{1}+90; %control zone
    elseif min(tPos) < 150 && min(tPos) > 0 %360 trials with position from 0 to 180
        posvals{1} = [36 64];
        posvals{2} = [64+9 64+9+54];
        posvals{3} = posvals{1}+90; 
    else
        posvals{1} = [-18 9]; %RRZ
        posvals{2} = [18 72]; %middle of the track speed split area
        posvals{3} = posvals{1}+90; %control zone
    end
    
    %%% get PCR and average sequence in the position chunks defined above
    for p = 1:length(posvals)
        isRZTrough = tPos >= posvals{p}(1) & tPos < posvals{p}(2);
        tRZSeq = mean(tSeq(:,:,isRZTrough), 3, 'omitnan');
        
        PCR_loc_trial{p}(t,1) = cf_calcquadrantratio(qX, qY, tRZSeq, 'mod');
        mn_seq_pos_trial{p}{t} = tRZSeq;
        sig_other_pos{p}(t) = getZoneAboveThresh(tRZSeq, params.nDeg);
        
        clear isRZTrough tRZSeq
    end
    
    %%% get PCR over space
    [trough_counts_trial(t,:),~,binID] = histcounts(tPos, dec_edges);
    for b = 1:length(dec_edges)-1
        isB = binID == b;
        bSeq = tSeq(:,:,isB);
        PCR_bpos(b) = cf_calcquadrantratio(qX,qY,mean(bSeq,3,'omitnan'), 'mod');
        QR_pos(b) = cf_calcquadrantratio(qX,qY,mean(bSeq,3,'omitnan'), 'full');
    end
    PCR_pos_trial = [PCR_pos_trial; PCR_bpos];
    QR_pos_trial = [QR_pos_trial; QR_pos];

    clear PCR_bpos isT tSeq tPos PCR_bpos QR_pos tRZSeq
end
mn_seq_reg_trial = mn_seq_reg_trial';
PCR_RRZ_trial = PCR_loc_trial{1};
PCR_speed_trial = PCR_loc_trial{2};
PCR_ctrl_trial = PCR_loc_trial{3};

mn_seq_RRZ_trial = mn_seq_pos_trial{1}';
mn_seq_speed_trial = mn_seq_pos_trial{2}';
mn_seq_ctrl_trial = mn_seq_pos_trial{3}';

sig_otherpos_RRZ = sig_other_pos{1}';
sig_otherpos_speed = sig_other_pos{2}';
sig_otherpos_ctrl = sig_other_pos{3}';

%out = table(PCR_pos_trial, PCR_loc_trial, trough_counts_trial, QR_pos_trial, mn_seq_reg_trial); 
out = table(PCR_pos_trial, PCR_RRZ_trial, PCR_speed_trial, PCR_ctrl_trial, sig_otherpos_RRZ, ...
    sig_otherpos_speed, sig_otherpos_ctrl, trough_counts_trial, ...
    QR_pos_trial, mn_seq_RRZ_trial, mn_seq_speed_trial, mn_seq_ctrl_trial);

end


function out = getZoneAboveThresh(seq, nDeg)

%create a threshold to determine if there is "significant" decoding of
%the other reward zone
vectseq = reshape(seq, 1, []); %unwrap to get the mean
mProb = mean(vectseq); % overall mean
sdevProb = std(vectseq); %overall std dev

% does anything in the other reward region cross the threshold?
%not sure exactly where this should be... going to do 1 zone at the
%beginning and end
if nDeg == 360
    otherRZ = seq([1:9, 172:end],:);
elseif nDeg == 180
    otherRZ = seq([1:3, end-2:end],:);
end
vectOtherRZ = reshape(otherRZ,1,[]);
out = any(vectOtherRZ > (2*sdevProb+mProb));

end
