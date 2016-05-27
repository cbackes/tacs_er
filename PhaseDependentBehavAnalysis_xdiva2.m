expt = 'tacs_enc_xdiva';

% load behavioral data
dataPath    = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/'];
load([dataPath 'Summary/BehavSummary.mat'])
addpath CircStats
%%
out                         = [];
out.nEncTrials              = 300;
out.StimPhases              = 0:2*pi/5:2*pi-0.01;
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
    'Hit','Miss','Confidence','MemScore','RTs2','RetPos','ConfWeStimVects',...
    'RTWeStimVects','MemScore2','ConfWeStimVects2'};
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
% 15th   column: MemScore2 (demeaned memscore such that it is zero for each subject)
% 16th   column: MemScore2*phase

out.datMat                  = nan(out.nSubjs,out.nEncTrials,numel(out.datMatColumnNames));
out.datMat(:,:,1)           = behav_out.Trials.StimTypeEnc;
out.datMat(:,:,5)           = behav_out.Trials.EncTrialPhaseCond;
out.datMat(:,:,6)           = out.EventPhase;
out.datMat(:,:,12)          = out.RetStimAtEnc;

% Pre-allocation
out.HitMissEncTrialIDs      = nan(nSubjs,2,300);
out.ValidEncTrialIDs        = nan(nSubjs,300);
out.HitMissRTSubj           = cell(nSubjs,2);   % RTs per condition

% Splitting of Phases: first column is hits, second misses
out.HM_Phases               = cell(nSubjs,2);   % phases for each trial cond
out.HM_PhaseCond            = cell(nSubjs,2);   % phase condition (1 to 5)
out.ConfPhases              = cell(nSubjs,3);   % phases per confidence
out.HM_PhaseConf            = cell(nSubjs,2,3);  % phases per confidenc
out.FaScnPhases             = cell(nSubjs,2);   % phases for each faces and scenes
out.HM_FaScn_Phases         = cell(nSubjs,2,2);   % phases for hits for separated by face and scenes

% Mean Vectors
out.HM_MeVects              = nan(nSubjs,2,2);    % mean phase per cond
out.HM_CW_MeVects           = nan(nSubjs,2,2);    % mean confidence weighted phase
out.HM_Conf_MeVects         = nan(nSubjs,2,3,2);  % mean phase per conf and condition
out.FaScn_MeVects           = nan(nSubjs,2,2);      % mean phase per cond
out.HM_FaScn_MeVects        = nan(nSubjs,2,2,2);  % mean phases for hits for separated by face and scenes
out.HM_Conf_FaScn_MeVects   = nan(nSubjs,3,2,2,2);  % mean phases for hits for separated by face and scenes
out.Conf_MeVects            = nan(nSubjs,3,2);
out.CW_MeVects              = nan(nSubjs,2);

% Memory Measures
out.HR                      = nan(nSubjs,1);
out.HRStimType              = nan(nSubjs,2);
out.MemScore                = nan(nSubjs,1);
out.MemScoreM               = nan(nSubjs,1);
out.MemScoreStimType        = nan(nSubjs,2);

% Memory Measures by Phase
out.MemScoreByPhase         = nan(nSubjs,5);
out.MemScore2ByPhase        = nan(nSubjs,5);
out.HRByPhase               = nan(nSubjs,5);
out.nHitsByPhase            = nan(nSubjs,5);
out.propHitsByPhase         = nan(nSubjs,5);
out.nMissByPhase            = nan(nSubjs,5);
out.propMissByPhase         = nan(nSubjs,5);

out.QuantileMarks           = [0.33 0.66];
%%
for ss = 1:out.nSubjs
    try
        %% Additional entries for data matrix
        %Get retrieval RTs
        out.datMat(ss,:,11)         = behav_out.retSubj{ss}.RTs(out.RetStimAtEnc(ss,:));
        
        % Trial IDs for hits and misses.
        HM      = false(300,2);
        HM(:,1) = behav_out.retSubj{ss}.Hits(out.RetStimAtEnc(ss,:));   % hits IDs at encoding
        HM(:,2) = behav_out.retSubj{ss}.Misses(out.RetStimAtEnc(ss,:)); % miss IDs at encoding
        out.datMat(ss,:,7:8)      = HM;
        
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
        out.datMat(ss,HM(:,2),10)   = -1*out.datMat(ss,HM(:,2),10);
        out.datMat(ss,:,15)         = out.datMat(ss,:,10) - nanmean(out.datMat(ss,:,10));
        
        % Confidence Weighted Phases
        out.datMat(ss,:,13) = out.datMat(ss,:,9).*exp(1i*out.datMat(ss,:,6));
        out.datMat(ss,:,16) = out.datMat(ss,:,15).*exp(1i*out.datMat(ss,:,6));
        
        % Weight Phases by RT
        out.datMat(ss,:,14) = (-log10(out.datMat(ss,:,11))).*exp(1i*out.datMat(ss,:,6));
        
        % Face/Scenes
        FS      = false(300,2);
        FS(:,1) =out.datMat(ss,:,1)==1;
        FS(:,2) =out.datMat(ss,:,1)==2;
        
        % Confidence
        CO      = false(300,3);
        CO(:,1) = out.datMat(ss,:,9)==1; % Low
        CO(:,2) = out.datMat(ss,:,9)==2; % Med
        CO(:,3) = out.datMat(ss,:,9)==3; % Hgh
        
        % Phases
        Phases = out.EventPhase(ss,:);
        % Retrieval RTs Quantiles
        RetRTs      = out.datMat(ss,:,11);
        RetRTsQIDs  = false(2,3,300);
        for jj = 1:2
            RetRTsQ     = quantile(RetRTs(HM(:,jj)),out.QuantileMarks);
            [c,i]=histc(RetRTs(HM(:,jj)),[0 RetRTsQ inf]);
            for qq = 1:3
                RetRTsQIDs(jj,qq,HM(:,jj)) = i==qq;
            end
            out.RetRTsQ(ss,jj,:) = RetRTsQ;
        end
        
        %RetRTsQIDs  = [RetRTs<RetRTsQ(1) RetRTs>RetRTsQ(1)]
        
        %% Get phases by confidence
        for co =1:3
            trials = CO(:,co);
            x = exp(1j*Phases(trials));
            mx = nanmean(x);
            out.ConfPhases{ss,co} = x;
            out.ConfPhases_N(ss,co) = numel(x);
            out.Conf_MeVects(ss,co,:) = [angle(mx), abs(mx)];
        end
        
        %% Get Mean Phase and Mean Vector length for hits and misses
        xm  = zeros(2,1);
        xmR = zeros(2,1);
        for jj = 1:2
            out.HM_Phases{ss,jj} = Phases(HM(:,jj));
            x = exp(1i*out.HM_Phases{ss,jj}');
            xm(jj) = nanmean(x);
            out.HM_Phases_N(ss,jj)  = numel(x);
            out.HM_MeVects(ss,jj,:) = [angle(xm(jj)),abs(xm(jj))];
            out.HM_R(ss,jj)         = numel(x)*abs(xm(jj));
            xmR(jj)                 = out.HM_R(ss,jj)*exp(1j*angle(xm(jj)));
        end
        out.HM_Z(ss,:)    = [angle(xm(1)-xm(2)),abs(xm(1)-xm(2))];
        out.HM_ZR(ss,:) = [angle(xmR(1)-xmR(2)),abs(xmR(1)-xmR(2))];

        
        xm  = zeros(2,3);
        xmR = zeros(2,3);
        for jj=1:2
            for co = 1:3
                x = exp(1j*Phases( HM(:,jj) & CO(:,co)));
                xm(jj,co) = nanmean(x);
                out.HM_Conf_N(ss,jj,co) = numel(x);
                out.HM_Conf_MeVects(ss,jj,co,:) = [angle(xm(jj,co)), abs(xm(jj,co))];
                out.HM_Conf_R(ss,jj,co)         = numel(x)*abs(xm(jj,co));
                xmR(jj,co)                      = out.HM_Conf_R(ss,jj,co)*exp(1j*angle(xm(jj,co)));
            end % By Confidence
        end % Hit/Miss
        for co = 1:3
            out.HM_Conf_Z(ss,co,1)  = angle(xm(1,co)-xm(2,co));
            out.HM_Conf_ZR(ss,co,1) = angle(xmR(1,co)-xmR(2,co));
            out.HM_Conf_Z(ss,co,2)  = abs(xm(1,co)-xm(2,co));
            out.HM_Conf_ZR(ss,co,2) = abs(xmR(1,co)-xmR(2,co));
        end
        
        %% Face/Scenes Phases
        xm = zeros(2,2);
        ym = zeros(2,2,3);
        for kk = 1:2
            % all independent of memory
            
            out.FaScnPhases{ss,kk} = Phases(FS(:,kk));
            x = nanmean(exp(1j*out.FaScnPhases{ss,kk}));
            out.FaScn_MeVects(ss,kk,:) = [angle(x),abs(x)];
            for jj=1:2
                HM_trials                       = FS(:,kk) & HM(:,jj);
                out.HM_FaScn_Phases{ss,jj,kk}   = Phases(HM_trials);
                
                x                                   = exp(1j*Phases(HM_trials));
                xm(jj,kk)                           = nanmean(x);
                out.HM_FaScn_N(ss,jj,kk)            = numel(x);
                out.HM_FaScn_MeVects(ss,jj,kk,:)    = [angle(xm(jj,kk)) abs(xm(jj,kk))];
                out.HM_FaScn_R(ss,jj,kk)            = numel(x)*abs(xm(jj,kk));
                for co = 1:3
                    Conf_trials         = HM_trials & CO(:,co);
                    x                   = exp(1j*Phases(Conf_trials));
                    xm2                 = nanmean(x);
                    ym(jj,co,kk)        = xm2;
                    out.HM_Conf_FaScn_MeVects(ss,jj,co,kk,:)    = [angle(xm2) abs(xm2)];
                    out.HM_Conf_FaScn_N(ss,jj,co,kk)            = numel(x);
                    out.HM_Conf_FaScn_R(ss,jj,co,kk)            = numel(x)*abs(xm2);
                end
            end % for hits/misses
            
        end % Face/Scenes
        for kk = 1:2
            out.HM_FaScn_Z(ss,kk) = abs(xm(1,kk)-xm(2,kk));
            for co = 1:3
                out.HM_Conf_FaScn_Z(ss,co,kk) = abs(ym(1,co,kk)-ym(2,co,kk));
            end
        end
        %% Confidence Weighted Vectors
        % Get Mean Phases weighted by confidence.
        xm = zeros(2,1);
        for jj=1:2
            x = out.datMat(ss,HM(:,jj),13);
            xm(jj) = nanmean(x);
            out.HM_CW_MeVects(ss,jj,:) = [angle(xm(jj)),abs(xm(jj))];
        end
        out.HM_CW_Z(ss) = abs(xm(1)-xm(2));
        
        x = nanmean(out.datMat(ss,:,16));
        out.CWMeVec(ss,:) = [angle(x) abs(x)];
        
        %% Get number of trials per phase condition
        PhaseConds=out.EventPhaseCond(ss,:)';
        r = histc(PhaseConds(HM(:,1)),1:5);
        out.nHitsByPhase(ss,:) = r;
        f = histc(PhaseConds(HM(:,2)),1:5);
        out.nMissByPhase(ss,:) = f;
        out.HRByPhase(ss,:)  = r./(r+f);
        
        % Hit Rate by Phase for stim categories
        for kk = 1:2
            r = histc(PhaseConds(HM(:,1)&FS(:,kk)),1:5); r=r(:);
            f = histc(PhaseConds(HM(:,2)&FS(:,kk)),1:5); f=f(:);
            out.HR_FS_Phase(ss,kk,:)  = r./(r+f);
        end
        
        % MemScore by Phase
        for pp = 1:5
            trials  = PhaseConds==pp;
            H_trials = trials & HM(:,1);
            M_trials = trials & HM(:,2);
            out.MemScoreByPhase(ss,pp)= mean(out.datMat(ss,trials,10),'omitnan');
            out.MemScore2ByPhase(ss,pp)= mean(out.datMat(ss,trials,15),'omitnan');
        end

        % HR by Confidence and Phase
        for co = 1:3
            H_trials = HM(:,1) & CO(:,co);
            r = histc(PhaseConds(H_trials),1:5); r=r(:);
            out.nHitsConfByPhase(ss,co,:) = r;
            
            M_trials = HM(:,2) & CO(:,co);
            f = histc(PhaseConds(M_trials),1:5); f=f(:);
            out.nMissConfByPhase(ss,co,:) = f;
            out.HR_Conf_Phase(ss,co,:)  = r./(r+f);
        end

        
        %% Get mean vectors per retrieval Rts quantiles
        
        xm = zeros(2,3);
        for jj=1:2
            for qq = 1:3
                trials = RetRTsQIDs(jj,qq,:);
                x           = exp(1j*Phases(trials));
                xm(jj,qq)   = nanmean(x);
                out.HM_RetRTsQ_MeVects(ss,jj,qq,:)  = [angle(xm(jj,qq)),abs(xm(jj,qq))];
                out.HM_RetRTsQ_N(ss,jj,qq)          = numel(x);
                out.HM_RetRTsQ_R(ss,jj,qq)          = numel(x)*abs(xm(jj,qq));                 
            end
            out.HM_RetRTsQ_Z(ss,:,:) = xm;
        end
        
        %%
        
        %         % Get R-F differences
        %         x=out.EventPhaseCond(ss,:);
        %         r=histc(x(out.datMat(ss,:,7)==1),1:5);
        %         f=histc(x(out.datMat(ss,:,8)==1),1:5);
        %         out.SubjRemNums(ss,:)       = r;
        %         rr=r./sum(r);
        %         out.SubjRemRate(ss,:)       = rr;
        %          z = mean(rr.*exp(1j*(out.StimPhases)));
        %         out.SubjRemRateVec(ss,:)     = [angle(z), abs(z)];
        %
        %         out.SubjForgNums(ss,:)      = f;
        %         ff=f/sum(f);
        %         out.SubjForgRate(ss,:)       = ff;
        %         z = mean(ff.*exp(1j*(out.StimPhases)));
        %         out.SubjForgRateVec(ss,:)     = [angle(z), abs(z)];
        %
        %         z = mean((r./(f+r)).*exp(1j*(out.StimPhases)));
        %         out.SubjRFMeanVec(ss,:)       = [angle(z) abs(z)];
        %
        %         % Get RF Differences High Confidence
        %         trials  =  out.datMat(ss,:,9)>1;
        %         r=histc(x(out.datMat(ss,:,7)==1 & trials),1:5);
        %         f=histc(x(out.datMat(ss,:,8)==1 & trials),1:5);
        %         out.SubjHMCRemNums(ss,:)       = r;
        %         rr = r./sum(r);
        %         out.SubjHMCRemRate(ss,:)       = rr;
        %         z = mean(rr.*exp(1j*(out.StimPhases)));
        %         out.SubjHMCRemRateVec(ss,:)     = [angle(z), abs(z)];
        %
        %         out.SubjHMCForgNums(ss,:)      = f;
        %         ff = f./sum(f);
        %         out.SubjHMCForgRate(ss,:)       = ff;
        %         z = mean(ff.*exp(1j*(out.StimPhases)));
        %         out.SubjHMCForgRateVec(ss,:)     = [angle(z), abs(z)];
        %
        %         z = mean((r./(f+r)).*exp(1j*(out.StimPhases)));
        %         out.SubjHMC_HR_MeanVec(ss,:)       = [angle(z) abs(z)];
        
    catch msg
        warning(sprintf('on subject %i',ss))
        disp(msg)
    end
end

save([dataPath 'Summary/PhaseDependentAnalyses.mat'],'out')
