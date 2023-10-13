function out = cf_getcodingratio_thetaseq(decodedData, posvals, params)
%cf_getcodingratio_thetaseq()
%
%ALP 1/18/2023

qX = params.qX;
qY = params.qY; 
dec_edges = params.dec_edges; %bins for calculating the PCR over the track

%%% get per trial PCR
tmpNTrials = max(unique(decodedData.trialNum))); 

PCR_pos_trial = []; PCR_loc_trial = [];
for t = 1:tmpNTrials
    isT = decodedData.trialID == t;
    tSeq = decodedData.allDecodedSeq(:,:,isT);
    tPos = decodedData.troughPos(isT);
    
    %%% get PCR chunk around RZ
    isRZTrough = tPos >= posvals{1}(1) & tPos < posvals{1}(2);
    tRZSeq = mean(tSeq(:,:,isRZTrough), 3, 'omitnan');
    PCR_loc_trial(t) = cf_calcquadrantratio(qX, qY, tRZSeq, 'mod');
    
    %%% get PCR over space
    [trough_counts_trial(t,:),~,binID] = histcounts(tPos, dec_edges);
    for b = 1:length(dec_edges)-1
        isB = binID == b;
        bSeq = tSeq(:,:,isB);
        PCR_bpos(b) = cf_calcquadrantratio(qX,qY,mean(bSeq,3,'omitnan'), 'mod');
    end
    PCR_pos_trial = [PCR_pos_trial; PCR_bpos];
    clear PCR_bpos isT tSeq tPos PCR_bpos
end

%%% get PCR in large chunks in particular areas of the track
PCR_loc = NaN(2,1);
for p = 1:length(posvals)
    isTrough = decodedData.troughPos >= posvals{p}(1) & decodedData.troughPos < posvals{p}(2);
    tmpSeq = decodedData.allDecodedSeq(:,:,isTrough);
    mnTmpSeq = mean(tmpSeq,3, 'omitnan');
    PCR_loc(p) = cf_calcquadrantratio(qX, qY, mnTmpSeq, 'mod');
    
    clear isTrough tmpSeq mnTmpSeq
end

%%% get PCR in small bins over the whole track
PCR_pos = NaN(1,length(dec_edges)-1);
[trough_counts_day, ~, binTrough] = histcounts(decodedData.troughPos, dec_edges);
for b = 1:length(dec_edges)-1
    isTrough = binTrough == b;
    tmpSeq = decodedData.allDecodedSeq(:,:,isTrough);
    mnTmpSeq = mean(tmpSeq, 3, 'omitnan');
    
    PCR_pos(b) = cf_calcquadrantratio(qX,qY, mnTmpSeq, 'mod');
    clear isTrough tmpSeq mnTmpSeq
end

%%% collect stuff for output

out.PCR_overposition_trial = PCR_pos_trial; 
out.PCR_overposition_day = PCR_pos; 
out.PCR_prerew_trial = PCR_loc_trial; 
out.PCR_prerew_day = PCR_loc(1); %1 is the approaching reward and in the reward zone
out.trough_counts_trial = trough_counts_trial; 
out.trough_counts_day = trough_counts_day; 

end

