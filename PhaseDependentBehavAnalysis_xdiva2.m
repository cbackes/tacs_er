
% load behavioral data
dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
load([dataPath 'Summary/BehavSummary.mat']) 
addpath CircStats
%%
nSubjs =30;

out                         = [];
out.nEncTrials              = 300;
out.SubjWithUniqueSeq       = sum(triu(corr(behav_out.Trials.EncCondAtRet')>0.99))==1;
out.nSubjs                  = nSubjs;
out.EventPhase              = behav_out.Trials.EncTrialPhase/360*2*pi; % phase at encoding
out.EventPhaseCond          = behav_out.Trials.EncTrialPhaseCond;
out.EventEncCondAtRet       = behav_out.Trials.EncCondAtRet; 
out.EncStimAtRet            = behav_out.Trials.EncStimIDAtRet;
out.RetStimAtEnc            = nan(nSubjs,out.nEncTrials);

out.datMatColumnNames      = ...
    {'StimType','PercepResp','FcScCorrect','RTs1','PhaseCond','PhaseDeg',...
    'Hit','Miss','Confidence','MemScore','RTs2','RetPos','ConfWeStimVects','RTWeStimVects'};
% matrix with all conditions and scores.
% 1st   column: face=1/scene=2
% 2nd   column: perceptual response (1/2)
% 3rd   column: correct face/scene categorization (bool)
% 4th   column: reaction times for perceptual decision
% 5rd   column: Phase Condition 
% 6th   column: Phase (degrees)
% 7th   column: Hits (binary) -> subsequently remember
% 8th   column: Misses (binary -> subsequenty forgotten
% 9th   column: Confidence-> subsequent confidence in hits /misses
% 10th   column: MemScore (re-scale confidence misses*-1)
% 11th   column: reaction times on memory decision
% 12th  column: position of item at retrieval
% 13th  column: confidence weighted phase trial vectors conf*exp(1j*phase)
% 14th  column: -log10(RT) weighted phase trial vectors rt*exp(1j*phase)

out.datMat                 =  nan(out.nSubjs,out.nEncTrials,numel(out.datMatColumnNames));
out.datMat(:,:,1)          =  behav_out.Trials.StimTypeEnc;
out.datMat(:,:,5)          =  behav_out.Trials.EncTrialPhaseCond;
out.datMat(:,:,6)          =  out.EventPhase;

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
        
    % Sort trials by their position at encoding.   
    EncIDsAtRet             = out.EncStimAtRet(ss,:);
    [s,i]                   = sort(EncIDsAtRet);
    out.RetIDsAtEnc(ss,:)   = i(s>0);
    out.datMat(ss,:,12)     = i(s>0);
    % Get retrieval RTs
    out.datMat(ss,:,11)         = behav_out.retSubj{ss}.RTs(out.RetIDsAtEnc(ss,:));
    
    % Trial IDs for hits and misses.
    out.datMat(ss,:,7)      = behav_out.retSubj{ss}.Hits(out.RetIDsAtEnc(ss,:));   % hits IDs at encoding
    out.datMat(ss,:,8)      = behav_out.retSubj{ss}.Misses(out.RetIDsAtEnc(ss,:)); % miss IDs at encoding
    
    % Get Confidence and assign sign based on hit/miss
    Confidence                  = behav_out.retSubj{ss}.Confidence;   
    Confidence(Confidence==0)   = nan;
    out.datMat(ss,:,9)          = Confidence(out.RetIDsAtEnc(ss,:));      
    out.datMat(ss,:,10)         = out.datMat(ss,:,9);
    out.datMat(ss,out.datMat(ss,:,8)==1,10) = -1*out.datMat(ss,out.datMat(ss,:,8)==1,10);
    
    out.datMat(ss,:,13) = out.datMat(ss,:,9).*exp(1i*out.datMat(ss,:,6));
    
    % Weight Phases by RT
    out.datMat(ss,:,14) = -log10(out.datMat(ss,:,11)).*exp(1i*out.datMat(ss,:,6));
    
    % Get phases per confidence
    for co =1:3
        trials = out.datMat(ss,:,9)==co;
        x = exp(1j*out.EventPhase(ss,trials));
        x(isnan(x))=[];
        out.ConfByPhaseMeVec(ss,co) = angle(mean(x));
        out.ConfByPhaseAbVec(ss,co) = abs(mean(x));
        [~,out.ConfByPhaseRStatVec(ss,co)] = circ_rtest(angle(x));
    end
    
    % Get true phase for Hits and Misses
    out.HitMissPhases{ss,1} = out.EventPhase(ss,out.datMat(ss,:,7)==1);    
    out.HitMissPhases{ss,2} = out.EventPhase(ss,out.datMat(ss,:,8)==1);    
    
    % Get phase condition for Hits and Misses
    out.HitMissPhaseCond{ss,1} = out.EventPhaseCond(ss,out.datMat(ss,:,7)==1);    
    out.HitMissPhaseCond{ss,2} = out.EventPhaseCond(ss,out.datMat(ss,:,8)==1);    
    % Get Mean Phase and Mean Vector length
    for jj = 1:2        
        x = exp(1i*out.HitMissPhases{ss,jj}');
        xm = mean(x);
        out.HitMissMePhase(ss,jj) = angle(xm);
        out.HitMissAbPhase(ss,jj) = abs(xm);
    end
    % Hit/Miss distance from uniform   
    [~,out.HitMiss_DistFromUniform(ss,1)]   = circ_rtest(out.HitMissPhases{ss,1});
    [~,out.HitMiss_DistFromUniform(ss,2)]   = circ_rtest(out.HitMissPhases{ss,2});
    
    % hit/miss difference in distributions
    [~,k]=circ_kuipertest(out.HitMissPhases{ss,1}',out.HitMissPhases{ss,2}');
    out.HitMiss_DistDiff(ss) = k;
    
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
            [~,out.HitMissConf_DistDiff(ss,co)]  = circ_kuipertest(out.HitMissPhaseConf{ss,1,co}',out.HitMissPhaseConf{ss,2,co}');
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
    [~,out.HitMissHMedConf_DistDiff(ss)]  = circ_kuipertest(angle(x{1}),angle(x{2}));

    %----%----%----%----%----%----%----%----%----%----%----%----%----%----
    %----%----%----%----%----%----%----%----%----%----%----%----%----%----
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
    
    % Get MemScore
    out.MemScore(ss)  = sum(out.datMat(ss,:,10),'omitnan');
    out.MemScoreM(ss) = mean(out.datMat(ss,:,10),'omitnan');
    
    % GetMemScore by stim type: face/scene
    for jj=1:2
        trials = squeeze(out.datMat(ss,:,1)==jj)==1;
        out.MemScoreStimType(ss,jj)= mean(out.datMat(ss,trials,10),'omitnan');
    end
    
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
    
    % Weight Phases by RT    
    for jj=1:2
        trials = squeeze(out.datMat(ss,:,6+jj))==1;
        x = mean(out.datMat(ss,trials,14),'omitnan');
        out.HitMissRTWMePhase(ss,jj) = angle(x);
        out.HitMissRTWAbPhase(ss,jj) = abs(x);
    end    
end

%%

% subject by subject phased shift
% th = mod(out.HitMissMePhase,2*pi);
% rho = out.HitMissAbPhase;
% z = rho.*exp(1i*th);
% 
% out2 =[];
% out2.HMissTheta = mod(out.HitMissMePhase,2*pi);
% out2.HMissRho   =  out.HitMissAbPhase;
% out2.HMiss
% 


save([dataPath 'Summary/PhaseDependentAnalyses.mat'],'out') 


%%

% %%
% % circularStatsByCond       = cell(nSubjs,2);
% % HitsFixedEffects        = []; 
% % CRsFixedEffects         = [];
% 
% % startEdge               = tacs_er.EncPhases;
% % binWidth                = 2*pi/10;
% % %phaseEdges              = (-pi:binWidth:pi);
% % nPhaseBins              = numel(startEdge);
% % conds                   = 1:5;
% % %phaseEdges
% % ConfBySubj              = cell(nSubjs,1);
% % PhasesBySubj            = cell(nSubjs,1);
% % PhasesBySubjCounts      = zeros(nSubjs,nPhaseBins);
% % PhaseBySubjCond1        = zeros(nSubjs,nPhaseBins);
% % PhaseBySubjCond2        = zeros(nSubjs,nPhaseBins);
% 
% % propHits        = nan(nSubjs,nPhaseBins);
% % HitRateByPhase  = nan(nSubjs,nPhaseBins);
% % P               = nan(nSubjs,nPhaseBins);
% % ConfByPhase     = cell(nSubjs,1);
% % ConfByPhaseNorm = cell(nSubjs,1);
% 
% for ss = 1:nSubjs
%     
%     EncStimAtRet    = behav_out.retSubj{ss}.EncStimIDAtRet;
%     hits            = behav_out.retSubj{ss}.Hits;
%     misses          = behav_out.retSubj{ss}.Misses;
%     
%     Confidence                  = behav_out.retSubj{ss}.Confidence;
%     Confidence(misses)          = Confidence(misses)*-1;
%     Confidence(~(hits|misses))  = [];            
%     
%     EncStimAtRet(~(hits|misses))=[];
%     ConfBySubj{ss} = Confidence;
%     
%     % Hits
%     PhasesPerSubj{ss,1} = EventPhase(ss,behav_out.EncRet.EncHitTrialIDs{ss});    
%     CondPerSubj{ss,1} = EventPhaseCond(ss,behav_out.EncRet.EncHitTrialIDs{ss});    
%     %temp = histc(PhasesPerCondSubj{ss,1},-pi:binWidth:pi);
%     %PhaseBySubjCond1(ss,:) = temp(1:end-1);    
%     
%     %Misses
%     %PhasesPerCondSubj{ss,2} = subjEventPhases{ss}.out.TrueAngleStims(behav_out.EncRet.EncMissTrialIDs{ss});
%     PhasesPerSubj{ss,2} = EventPhase(ss,behav_out.EncRet.EncMissTrialIDs{ss});
%     CondPerSubj{ss,2} = EventPhaseCond(ss,behav_out.EncRet.EncMissTrialIDs{ss});    
%     %temp = histc(PhasesPerCondSubj{ss,2},-pi:binWidth:pi);
%     %PhaseBySubjCond2(ss,:) = temp(1:end-1);
% 
%     %
%     PhasesBySubj{ss} = EventPhase(ss,EncStimAtRet);   
%     CondBySubj{ss}   = EventPhaseCond(ss,EncStimAtRet);    
%     %temp = histc(PhasesBySubj{ss},-pi:binWidth:pi);
%     %PhasesBySubjCounts(ss,:)   = temp(1:end-1);
%     
%     [~,HCRDistFromUniform(ss,1)]   = circ_rtest(PhasesPerSubj{ss,1});
%     [~,HCRDistFromUniform(ss,2)]   = circ_rtest(PhasesPerSubj{ss,2});
%     
%     circularStatsByCond{ss,1}      = circ_stats(PhasesPerSubj{ss,1});
%     circularStatsByCond{ss,2}      = circ_stats(PhasesPerSubj{ss,2});
%     
%     me(ss,1)= circularStatsByCond{ss,1}.mean;
%     me(ss,2)= circularStatsByCond{ss,2}.mean;    
%        
%     CondRTSubj{ss,1}    = behav_out.EncRet.EncHitRTs{ss};
%     CondRTSubj{ss,2}    = behav_out.EncRet.EncMissRTs{ss};
%     
%     nHitsByPhase(ss,:)        = histc(CondPerSubj{ss,1},conds);
%     nMissesByPhase(ss,:)      = histc(CondPerSubj{ss,2},conds);
%     temp = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
%     HitRateByPhase(ss,:) = temp;
%     
%     temp = nHitsByPhase/numel(CondPerSubj{ss,1});
%     propHits(ss,:) = temp;
%     
% end
% 
% %% trial counts by phase
% 
% figure(1); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% 
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhasesBySubjCounts','-','color',[0.8 0.8 0.8])
% plot([30:60:360],mean(PhasesBySubjCounts), 'k','linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts' )
% 
% print(gcf, '-dpdf', '../plots/trialCountsByPhase_xdiva');
% 
% %% hits and misses raw #s
% 
% figure(2); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2);
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhaseBySubjCond1','-','color',[255 180 150]/255)
% plot([30:60:360],mean(PhaseBySubjCond1), 'color',[240 80 40]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts (hits) ' )
% 
% print(gcf, '-dpdf', '../plots/HitsCountsByPhase');
% 
% %
% 
% figure(3); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond1./repmat(sum(PhaseBySubjCond1,2),[1,6]);
% plot([30:60:360],X','-','color',[255 180 150]/255)
% plot([30:60:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0 0.3])
% ylabel(' proportion (hits) ' )
% 
% print(gcf, '-dpdf', '../plots/HitsPropByPhase');
% 
% %% hits and misses, proportions
% figure(4); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% t = linspace(0,1,1000);
% x = cos(2*pi*t-pi);
% xa = angle(hilbert(x));
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% 
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% plot([30:60:360],PhaseBySubjCond2','-','color',[180 250 250]/255)
% plot([30:60:360],mean(PhaseBySubjCond2), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylabel(' trial counts (misses) ' )
% 
% print(gcf, '-dpdf', '../plots/MissesCountsByPhase');
% 
% figure(5); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond2./repmat(sum(PhaseBySubjCond2,2),[1,6]);
% plot([30:60:360],X','-','color',[180 250 250]/255)
% plot([30:60:360],mean(X), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0 0.3])
% ylabel(' proportion (misses) ' )
% 
% print(gcf, '-dpdf', '../plots/MissesPropByPhase');
% 
% figure(6); clf;
% 
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% 
% a2 = axes('position',[0.12 0.1 0.8 0.15]);
% axes(a2)
% plot((xa+pi)./pi*180,x,'k','linewidth',4)
% axis tight; 
% set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:60:360])
% xlabel(' Encoding Phase (deg)')
% grid on
% 
% a1 = axes('position', [0.12 0.3 0.8 0.6]);
% axes(a1); hold on;
% X = PhaseBySubjCond1./repmat(sum(PhaseBySubjCond1,2),[1,6]);
% Y = PhaseBySubjCond2./repmat(sum(PhaseBySubjCond2,2),[1,6]);
% plot([30:60:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
% plot([30:60:360],mean(Y), 'color',[120 200 200]/255,'linewidth',5)
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% set(gca,'xtick',[0:30:360],'xTickLabel','')
% xlim([0 360])
% ylim([0.1 0.25])
% 
% ylabel(' proportion of trials ' )
% 
% print(gcf, '-dpdf', '../plots/HMPropByPhase');
% 
% %% distance from uniform
% figure(8); clf;
% set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
% 'paperposition',[0.1 0.2 0.6 0.4],'paperunits','normalized');
% a2 = axes('position',[0.15 0.1 0.3 0.8]); hold on;
% for ss = 1:nSubjs
%     h1 = plot([1 2],[me(ss,1) me(ss,2)]/pi*180,'o-','linewidth',2,'color',[0.5 0.5 0.5],...
%         'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
% end
% set(a2,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'hits','misses'})
% set(a2,'ytick',[-pi:pi/3:pi]./pi*180)
% ylabel(' mean Phase (º) ')
% plot([0.8 1.2], [1 1]*mean(me(:,1))/pi*180,'linewidth',5,'color','k')
% plot([1.8 2.2], [1 1]*mean(me(:,2))/pi*180,'linewidth',5,'color','k')
% grid on
% xlim([0.6 2.4])
% 
% a3 = axes('position',[0.6 0.1 0.35 0.8]); hold on;
% for ss = 1:nSubjs
%     h1 = plot([1 2],[HCRDistFromUniform(ss,1) HCRDistFromUniform(ss,2)],'o-','linewidth',2,'color',[0.5 0.5 0.5],...
%         'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
% end
% set(a3,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'D(H)','D(M)'})
% set(a3,'ytick',[0:0.5:4])
% ylim([0 4])
% ylabel(' Z-Stat ')
% plot([0.8 1.2], [1 1]*mean(HCRDistFromUniform(:,1)),'linewidth',5,'color','k')
% plot([1.8 2.2], [1 1]*mean(HCRDistFromUniform(:,2)),'linewidth',5,'color','k')
% 
% xlim([0.6 2.4])
% grid on
% print(gcf, '-dpdf', '../plots/DistFromUniformByCond');
% 
% %% 
% 
% %% modulation hit rate and memory score
% 
% 
% %%
% stepSize   = pi/6;
% phaseEdges = 0:%(-pi:stepSize:pi);
% nPhaseBins = 5;%numel(phaseEdges);
% 
% propHits        = nan(nSubjs,nPhaseBins);
% propMisses      = nan(nSubjs,nPhaseBins);
% HitRateByPhase  = nan(nSubjs,nPhaseBins);
% H_Hits          = nan(nSubjs,1);
% P               = nan(nSubjs,nPhaseBins);
% ConfByPhase     = cell(nSubjs,1);
% ConfByPhaseNorm = cell(nSubjs,1);
% 
% confScores = [-3 -2 -1 1 2 3]; nConfScores = numel(confScores);
% phaseMemScore = nan(nSubjs-1,nPhaseBins-1);
% for ss=1:nSubjs
%     nHitsByPhase        = histc(PhasesPerCondSubj{ss,1},conds);
%     nMissesByPhase      = histc(PhasesPerCondSubj{ss,2},conds);
%     HitRateByPhase(ss-1,:) = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
%     
%     propHits(ss-1,:) = nHitsByPhase/numel(PhasesPerSubj{ss,1});
%     P(ss-1,:) = HitRateByPhase(ss-1,:)-nanmean(HitRateByPhase(ss-1,:));    
%     
%     propMisses(ss-1,:) = nMissesByPhase/numel(PhasesPerSubj{ss,2});
%     
%     counter = 1;
%     for cc=confScores
%         ConfByPhase{ss-1}(counter,:) = histc(PhasesBySubj{ss}(ConfBySubj{ss}==cc),conds);                
%         counter = counter+1;
%     end    
%     [~, Idx ]=histc(PhasesBySubj{ss},phaseEdges);
%     for ph = 1:(nPhaseBins-1)
%         phaseMemScore(ss-1,ph)=sum(ConfBySubj{ss}(Idx==ph));
%     end
%     
%      ConfByPhase{ss-1}(:,nPhaseBins)=[];
%     ConfByPhaseNorm{ss-1}=ConfByPhase{ss-1}./repmat(sum(ConfByPhase{ss-1}),[nConfScores 1]);
%     ConfByPhaseNorm{ss-1} =  ConfByPhaseNorm{ss-1}./repmat(sum(ConfByPhaseNorm{ss-1},2),[1 nPhaseBins-1]);
% %     
% end
% 
% phaseMemScoreNorm=phaseMemScore-repmat(mean(phaseMemScore,2),[1 nPhaseBins-1]);
% 
% P(:,nPhaseBins)=[];
% HitRateByPhase(:,nPhaseBins)=[];
% 
% propHits(:,nPhaseBins)=[];
% propMisses(:,nPhaseBins)=[];
% 
% pu = 1/(nPhaseBins-1)*ones(nPhaseBins-1,1);
% sum(-pu.*log2(pu))
% %X = propHits-propMisses;
% 
% %plot(diff(cumsum(phaseEdges)),HitRateByPhase-repmat(mean(HitRateByPhase,2),[1,nPhaseBins-1]),'-*')
