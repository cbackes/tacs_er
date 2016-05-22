expt = 'tacs_enc_xdiva';

% load behavioral data
dataPath    = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
load([dataPath 'Summary/BehavSummary.mat'])
addpath CircStats
%%
out                         = [];
out.nEncTrials              = 300;
out.SubjWithUniqueSeq       = sum(triu(corr(behav_out.Trials.StimTypeEnc')>0.99))==1;
switch expt
    case 'tacs_enc_xdiva'
        nSubjs =32;
        out.EventPhaseCond          = behav_out.Trials.EncTrialPhaseCond;
        out.EventEncCondAtRet       = behav_out.Trials.EncCondAtRet;
        out.EventPhase              = behav_out.Trials.EncTrialPhase/360*2*pi; % phase at encoding
    case 'tacs_enc'
        nSubjs = 11;
        out.EventPhase              = behav_out.Trials.EncTrialPhase; % phase at encoding
    otherwise
        error('expt not valid')
end
out.nSubjs                  = nSubjs;
out.EncStimAtRet            = behav_out.Trials.EncStimIDAtRet;
out.RetStimAtEnc            = behav_out.Trials.RetStimIDAtEnc;

out.datMatColumnNames      = ...
    {'StimType','PercepResp','FcScCorrect','RTs1','PhaseCond','PhaseDeg',...
    'Hit','Miss','Confidence','MemScore','RTs2','RetPos','ConfWeStimVects','RTWeStimVects'};
% matrix with all conditions and scores.
% 1st   column: face=1/scene=2
% 2nd   column: perceptual response (1/2)
% 3rd   column: correct face/scene categorization (bool)
% 4th   column: reaction times for perceptual decision
% 5rd   column: Phase Condition
% 6th   column: Phase (rads)
% 7th   column: Hits [(binary) -> subsequently remember
% 8th   column: Misses (binary -> subsequenty forgotten
% 9th   column: Confidence-> subsequent confidence in hits /misses
% 10th   column: MemScore (re-scale confidence misses*-1)
% 11th   column: reaction times on memory decision
% 12th  column: position of item at retrieval
% 13th  column: confidence weighted phase trial vectors conf*exp(1j*phase)
% 14th  column: -log10(RT) weighted phase trial vectors rt*exp(1j*phase)

out.datMat                  = nan(out.nSubjs,out.nEncTrials,numel(out.datMatColumnNames));
out.datMat(:,:,1)           = behav_out.Trials.StimTypeEnc;
out.datMat(:,:,5)           = behav_out.Trials.EncTrialPhaseCond;
out.datMat(:,:,6)           = out.EventPhase;
out.datMat(:,:,12)          = out.RetStimAtEnc;
out.HitMissEncTrialIDs      = nan(nSubjs,2,300);
out.ValidEncTrialIDs        = nan(nSubjs,300);

% first column is hits, second misses
out.HitMissPhases           = cell(nSubjs,2);   % phases for each trial cond
out.HitMissMePhase          = nan(nSubjs,2);    % mean phase per cond
out.HitMissAbPhase          = nan(nSubjs,2);    % mean resulting vector length per condition
out.HitMissCWMePhase        = nan(nSubjs,2);    % mean confidence weighted phase
out.HitMissCWAbPhase        = nan(nSubjs,2);    % mean confidence vector length
out.HitMissPhaseCond        = cell(nSubjs,2);   % phase condition (1 to 5)
out.HitMissRTSubj           = cell(nSubjs,2);   % RTs per condition
out.PhasesByConf            = cell(nSubjs,3);   % phases per confidence
out.HitMissPhaseConf        = cell(nSubjs,2,3);  % phases per confidence and condition
out.HitMissMePhaseConf      = nan(nSubjs,2,3);   % mean phase per conf and condition
out.HitMissAbPhaseConf      = nan(nSubjs,2,3);   % mean vec length per conf and cond
out.HitMissCMTestConf       = nan(nSubjs,3);     % circular test for difference in phases.

% first column is faces, second scenes
out.FaScnPhases             = cell(nSubjs,2);   % phases for each faces and scenes
out.FaScnMePhases           = nan(nSubjs,2);    % mean phase per cond
out.FaScnAbPhases           = nan(nSubjs,2);    % mean resulting vector length per condition

out.FaScnHitPhases          = cell(nSubjs,2);   % phases for hits for separated by face and scenes
out.FaScnMissPhases         = cell(nSubjs,2);   % phases for misses for separated by face and scenes
out.FaScnHitMePhases        = nan(nSubjs,2);    % mean phases for hits for separated by face and scenes
out.FaScnMissMePhases       = nan(nSubjs,2);    % mean phases for misses for separated by face and scenes
out.FaScnHitAbPhases        = nan(nSubjs,2);    % vec length phases for hits for separated by face and scenes
out.FaScnMissAbPhases       = nan(nSubjs,2);    % vec length phases for misses for separated by face and scenes

out.MemScore                = nan(nSubjs,1);
out.MemScoreM               = nan(nSubjs,1);
out.MemScoreStimType        = nan(nSubjs,2);
out.MemScoreByPhase         = nan(nSubjs,5);
out.nHitsByPhase            = nan(nSubjs,5);
out.propHitsByPhase         = nan(nSubjs,5);
out.nMissByPhase            = nan(nSubjs,5);
out.propMissByPhase         = nan(nSubjs,5);

% circular distribution stats
out.HitMiss_DistFromUniform     = nan(nSubjs,2);
out.Conf_DistFromUniform        = nan(nSubjs,3);
out.HitsConf_DistFromUniform    = nan(nSubjs,3);
out.MissConf_DistFromUniform    = nan(nSubjs,3);
out.DiffInDist                  = nan(nSubjs,2);

for ss = 1:out.nSubjs
    try
        % Get retrieval RTs
        out.datMat(ss,:,11)         = behav_out.retSubj{ss}.RTs(out.RetStimAtEnc(ss,:));
        
        % Trial IDs for hits and misses.
        out.datMat(ss,:,7)      = behav_out.retSubj{ss}.Hits(out.RetStimAtEnc(ss,:));   % hits IDs at encoding
        out.datMat(ss,:,8)      = behav_out.retSubj{ss}.Misses(out.RetStimAtEnc(ss,:)); % miss IDs at encoding
        
        % Get MemScore
        out.MemScore(ss)  = sum(out.datMat(ss,:,10),'omitnan');
        out.MemScoreM(ss) = mean(out.datMat(ss,:,10),'omitnan');
        
        % GetMemScore by stim type: face/scene
        for jj=1:2
            trials = squeeze(out.datMat(ss,:,1)==jj)==1;
            out.MemScoreStimType(ss,jj)= mean(out.datMat(ss,trials,10),'omitnan');
        end
        
        % Get Confidence and assign sign based on hit/miss
        Confidence                  = behav_out.retSubj{ss}.Confidence;
        Confidence(Confidence==0)   = nan;
        out.datMat(ss,:,9)          = Confidence(out.RetStimAtEnc(ss,:));
        out.datMat(ss,:,10)         = out.datMat(ss,:,9);
        out.datMat(ss,out.datMat(ss,:,8)==1,10) = -1*out.datMat(ss,out.datMat(ss,:,8)==1,10);
        
        out.datMat(ss,:,13) = out.datMat(ss,:,9).*exp(1i*out.datMat(ss,:,6));
        
        % Weight Phases by RT
        out.datMat(ss,:,14) = (1./out.datMat(ss,:,11)).*exp(1i*out.datMat(ss,:,6));
        
        % Get phases per confidence
        for co =1:3
            trials = out.datMat(ss,:,9)==co;
            x = exp(1j*out.EventPhase(ss,trials));
            x(isnan(x))=[];
            out.ConfByPhaseMeVec(ss,co) = angle(mean(x));
            out.ConfByPhaseAbVec(ss,co) = abs(mean(x));
            [~,out.ConfByPhaseRStatVec(ss,co)] = circ_rtest(angle(x));
        end
        
        % Get Mean Phase and Mean Vector length for hits and phases
        for jj = 1:2
            % Get true phase for Hits and Misses
            out.HitMissPhases{ss,jj} = out.EventPhase(ss,out.datMat(ss,:,6+jj)==1);
            x = exp(1i*out.HitMissPhases{ss,jj}');
            xm = mean(x);
            out.HitMissMePhase(ss,jj) = angle(xm);
            out.HitMissAbPhase(ss,jj) = abs(xm);
            % distance
            [~,out.HitMiss_DistFromUniform(ss,jj)]   = circ_rtest(out.HitMissPhases{ss,jj});
        end % hits/misses 
        
        % Face/Scenes Phases
        for jj = 1:2
            % all independent of memory
            trials = out.datMat(ss,:,1)==jj;
            out.FaScnPhases{ss,jj} = out.EventPhase(ss,trials);
            x = mean(exp(1j*out.FaScnPhases{ss,jj}));
            out.FaScnMePhases(ss,jj) = angle(x);
            out.FaScnAbPhases(ss,jj) = abs(x);
            
            % for hits
            trials2 = trials & out.datMat(ss,:,7)==1;
            x = mean(exp(1j*out.EventPhase(ss,trials2)));
            out.FaScnHitMePhases(ss,jj)     = angle(x);
            out.FaScnHitAbPhases(ss,jj)     = abs(x);
            out.FaScnHitN(ss,jj)            = sum(trials2);
            
            % for misses
            trials3 = trials & out.datMat(ss,:,8)==1;
            x = mean(exp(1j*out.EventPhase(ss,trials3)));
            out.FaScnMissMePhases(ss,jj)     = angle(x);
            out.FaScnMissAbPhases(ss,jj)     = abs(x);
            out.FaScnMissN(ss,jj)            = sum(trials3);
        end % Face/Scenes
        
        % Get Mean Phases weighted by confidence.
        for jj=1:2
            trials = squeeze(out.datMat(ss,:,6+jj))==1;
            x = out.datMat(ss,trials,13);
            x(isnan(x))=[];
            out.HitMissCWMePhase(ss,jj) = angle(nanmean(x));
            out.HitMissCWAbPhase(ss,jj) = abs(nanmean(x));
            [~,out.HitMiss_CWDist(ss,jj)]   = circ_rtest(angle(x));
        end
        
        % Obtain mean phase and mean vector length per confidence level by
        % subject.
        for co=1:3
            try
                for jj=1:2
                    trials = squeeze(out.datMat(ss,:,9)==co & out.datMat(ss,:,6+jj)==1);
                    out.HitMissPhaseConf{ss,jj,co} = out.EventPhase(ss, trials);
                    x = exp(1i*out.HitMissPhaseConf{ss,jj,co}');
                    out.HitMissMePhaseConf(ss,jj,co)= angle(nanmean(x));
                    out.HitMissAbPhaseConf(ss,jj,co)= abs(nanmean(x));
                    [~,out.HitMiss_DistFromUniformConf(ss,jj,co)]   = circ_rtest(angle(x));
                end % for hits/misses
            catch
            end
        end % for confidence level.
        
        % low confidence hit miss
        for jj=1:2
            out.HitMissMePhaseLoConf(ss,jj)= out.HitMissMePhaseConf(ss,jj,1);
            out.HitMissAbPhaseLoConf(ss,jj)= out.HitMissAbPhaseConf(ss,jj,1);
        end
        
        % Obtain mean phase and mean vector length for confidence >=2 mid/high
        x=[];
        for jj=1:2
            trials = squeeze(out.datMat(ss,:,9)>1 & out.datMat(ss,:,6+jj)==1);
            phases = out.EventPhase(ss, trials);
            x{jj} = exp(1i*phases');
            out.HitMissMePhaseHMedConf(ss,jj)= angle(mean(x{jj}));
            out.HitMissAbPhaseHMedConf(ss,jj)= abs(mean(x{jj}));
            
            [~,out.HitMiss_DistFromUniformHMedConf(ss,jj)]   = circ_rtest(angle(x{jj}));
        end % for hits/misses
        
        %----%----%----%----%----%----%----%----%----%----%----%----%----%----
        %----%----%----%----%----%----%----%----%----%----%----%----%----%----
        if strcmp(expt,'tacs_enc_xdiva')
            % Get proportion of hits and misses by phase
            for pp = 1:5
                trials = squeeze(out.datMat(ss,:,5)==pp)&out.datMat(ss,:,7)==1;
                out.nHitsByPhase(ss,pp)= sum(trials);
                
                trials = squeeze(out.datMat(ss,:,5)==pp)&out.datMat(ss,:,8)==1;
                out.nMissByPhase(ss,pp)= sum(trials);
                out.HR_ByPhase(ss,pp) = out.nHitsByPhase(ss,pp)/(out.nMissByPhase(ss,pp) + out.nHitsByPhase(ss,pp));
                
                % for high and med confidence only.
                trials = squeeze(out.datMat(ss,:,5)==pp & out.datMat(ss,:,9)>1 & out.datMat(ss,:,7)==1);
                out.nHitsHMedConfByPhase(ss,pp)= sum(trials);
                
                trials = squeeze(out.datMat(ss,:,5)==pp & out.datMat(ss,:,9)>1 & out.datMat(ss,:,8)==1);
                out.nMissHMedConfByPhase(ss,pp)= sum(trials);
                out.HR_HMedConf_ByPhase(ss,pp) = out.nHitsHMedConfByPhase(ss,pp)/(out.nMissHMedConfByPhase(ss,pp) + out.nHitsHMedConfByPhase(ss,pp));
                
                % for low confidence trials.
                trials1 = squeeze(out.datMat(ss,:,5)==pp & out.datMat(ss,:,9)==1 & out.datMat(ss,:,7)==1);
                trials2 = squeeze(out.datMat(ss,:,5)==pp & out.datMat(ss,:,9)==1 & out.datMat(ss,:,8)==1);
                
                out.HR_LoConf_ByPhase(ss,pp) = sum(trials1)/(sum(trials2)+sum(trials1));
            end
            out.propHitsByPhase(ss,:) = out.nHitsByPhase(ss,:)/sum(out.nHitsByPhase(ss,:));
            out.propMissByPhase(ss,:) = out.nMissByPhase(ss,:)/sum(out.nMissByPhase(ss,:));
            % GetMemScore by phase
            for pp = 1:5
                trials = squeeze(out.datMat(ss,:,5)==pp)==1;
                out.MemScoreByPhase(ss,pp)= mean(out.datMat(ss,trials,10),'omitnan');
            end
            % RT analyses.
            % Get RTs by phase and condition
            for pp=1:5
                trials = squeeze(out.datMat(ss,:,5)==pp)&out.datMat(ss,:,7)==1;
                out.HitMissRTsByPhase(ss,pp,1)= mean(-log10(out.datMat(ss,trials,11)));
                trials = squeeze(out.datMat(ss,:,5)==pp)&out.datMat(ss,:,8)==1;
                out.HitMissRTsByPhase(ss,pp,2)= mean(-log10(out.datMat(ss,trials,11)));
            end
        end
        
        % Weight Phases by RT
        for jj=1:2
            trials = squeeze(out.datMat(ss,:,6+jj))==1;
            x = mean(out.datMat(ss,trials,14),'omitnan');
            out.HitMissRTWMePhase(ss,jj) = angle(x);
            out.HitMissRTWAbPhase(ss,jj) = abs(x);
        end
    catch msg
        warning(sprintf('on subject %i',ss))
        disp(msg)
    end
end

save([dataPath 'Summary/PhaseDependentAnalyses.mat'],'out')
