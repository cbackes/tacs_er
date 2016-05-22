dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc/';
subjs       = 2:11;
nSubjs      = numel(subjs);

% load behavioral data
load([dataPath 'Summary/BehavSummary.mat']) 

subjEventPhases =cell(nSubjs,1);
for ss = 1:nSubjs
    load([dataPath 's' num2str(subjs(ss)) '/EventsPhase.mat']);
    subjEventPhases{ss} = out;
end

%%

PhasesPerCondSubj   = cell(nSubjs,2); % first column is hits, second misses
PhasesPerConfSubj   = cell(nSubjs,3); % high to low
HCRDistFromUniform  = nan(nSubjs,2);
ConfDistFromUniform  = nan(nSubjs,2);
DiffInDist          = nan(nSubjs,2);
CondRTSubj          = cell(nSubjs,2); % first column is hits, second misses

circularStatsByCond       = cell(nSubjs,2);
HitsFixedEffects        = []; 
CRsFixedEffects         = [];

binWidth                = pi/3;
phaseEdges = (-pi:binWidth:pi);
nPhaseBins = numel(phaseEdges)-1;

ConfBySubj              = cell(nSubjs,1);
PhasesBySubj            = cell(nSubjs,1);
PhasesBySubjCounts      = zeros(nSubjs,nPhaseBins);
PhaseBySubjCond1        = zeros(nSubjs,nPhaseBins);
PhaseBySubjCond2        = zeros(nSubjs,nPhaseBins);

propHits        = nan(nSubjs,nPhaseBins);
HitRateByPhase  = nan(nSubjs,nPhaseBins);
P               = nan(nSubjs,nPhaseBins);
ConfByPhase     = cell(nSubjs,1);
ConfByPhaseNorm = cell(nSubjs,1);

%%
for ss = 1:nSubjs
    subj = subjs(ss);
    EncStimAtRet    = behav_out.retSubj{subj}.EncStimIDAtRet;
    hits            = behav_out.retSubj{subj}.Hits;
    misses          = behav_out.retSubj{subj}.Misses;
    
    Confidence                  = behav_out.retSubj{subj}.Confidence;
    Confidence(misses)          = Confidence(misses)*-1;
    Confidence(~(hits|misses))  = [];        
    
    
    EncStimAtRet(~(hits|misses))=[];
    ConfBySubj{subj} = Confidence;
    
    PhasesPerCondSubj{ss,1} = subjEventPhases{ss}.out.TrueAngleStims(behav_out.EncRet.EncHitTrialIDs{subj});
    temp = histc(PhasesPerCondSubj{ss,1},-pi:binWidth:pi);
    PhaseBySubjCond1(ss,:) = temp(1:end-1);    
    
    PhasesPerCondSubj{ss,2} = subjEventPhases{ss}.out.TrueAngleStims(behav_out.EncRet.EncMissTrialIDs{subj});
    temp = histc(PhasesPerCondSubj{ss,2},-pi:binWidth:pi);
    PhaseBySubjCond2(ss,:) = temp(1:end-1);

    PhasesBySubj{ss}        = subjEventPhases{ss}.out.TrueAngleStims(EncStimAtRet);   
    temp = histc(PhasesBySubj{ss},-pi:binWidth:pi);
    PhasesBySubjCounts(ss,:)   = temp(1:end-1);
    
    [~,HCRDistFromUniform(ss,1)]   = circ_rtest(PhasesPerCondSubj{ss,1});
    [~,HCRDistFromUniform(ss,2)]   = circ_rtest(PhasesPerCondSubj{ss,2});
    
    circularStatsByCond{ss,1}      = circ_stats(PhasesPerCondSubj{ss,1});
    circularStatsByCond{ss,2}      = circ_stats(PhasesPerCondSubj{ss,2});
    
    me(ss,1)= circularStatsByCond{ss,1}.mean;
    me(ss,2)= circularStatsByCond{ss,2}.mean;    
       
    CondRTSubj{ss,1}    = behav_out.EncRet.EncHitRTs{subj};
    CondRTSubj{ss,2}    = behav_out.EncRet.EncMissRTs{subj};
    
    nHitsByPhase        = histc(PhasesPerCondSubj{ss,1},phaseEdges);
    nMissesByPhase      = histc(PhasesPerCondSubj{ss,2},phaseEdges);
    temp = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
    HitRateByPhase(ss,:) = temp(1:end-1);
    
    temp = nHitsByPhase/numel(PhasesPerCondSubj{ss,1});
    propHits(ss,:) = temp(1:end-1);
    
end
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
% print(gcf, '-dpdf', '../plots/trialCountsByPhase');
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
% phaseEdges = (-pi:stepSize:pi);
% nPhaseBins = numel(phaseEdges);
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
% for ss=2:nSubjs
%     nHitsByPhase        = histc(PhasesPerCondSubj{ss,1},phaseEdges);
%     nMissesByPhase      = histc(PhasesPerCondSubj{ss,2},phaseEdges);
%     HitRateByPhase(ss-1,:) = nHitsByPhase./(nHitsByPhase+nMissesByPhase);
%     
%     propHits(ss-1,:) = nHitsByPhase/numel(PhasesPerCondSubj{ss,1});
%     P(ss-1,:) = HitRateByPhase(ss-1,:)-nanmean(HitRateByPhase(ss-1,:));    
%     
%     propMisses(ss-1,:) = nMissesByPhase/numel(PhasesPerCondSubj{ss,2});
%     
%     counter = 1;
%     for cc=confScores
%         ConfByPhase{ss-1}(counter,:) = histc(PhasesBySubj{ss}(ConfBySubj{ss}==cc),phaseEdges);                
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
