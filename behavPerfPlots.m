


%% categorization performance 
mACC = behav_out.encSummary.meanAcc;
RTs  = behav_out.encSummary.meanRTs;
nSubjs = numel(mACC);

figure(1);clf;
set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
'paperposition',[0.2 0.2 0.6 0.4],'paperunits','normalized');

a1 = axes('position',[0.1 0.1 0.25 0.8]); hold on;
h1 = scatter(ones(nSubjs,1)+0.05*randn(nSubjs,1),mACC);
set(h1,'markerfacecolor',[0.5 0.5 0.5],'markeredgeColor','k','sizeData',70)
set(a1,'fontSize',16,'lineWidth',3,'xtick',[])
set(a1,'ytick',[0.5:0.1:1])
ylabel(' mean Accuracy ')

M_mACC = mean(mACC);

plot([0.9 1.1],[M_mACC M_mACC],'linewidth',4,'color','k')
xlim([0.8 1.2])

grid on

a2 = axes('position',[0.6 0.1 0.25 0.8]); hold on;
h1 = scatter(ones(nSubjs,1)+0.05*randn(nSubjs,1),RTs);
set(h1,'markerfacecolor',[0.5 0.5 0.5],'markeredgeColor','k','sizeData',70)
set(a2,'fontSize',16,'lineWidth',3,'xtick',[])
set(a2,'ytick',[0.5:0.1:1])
ylabel(' mean RTs (s) ')

M_mRTs = mean(RTs);

plot([0.9 1.1],[M_mRTs M_mRTs],'linewidth',4,'color','k')
xlim([0.8 1.2])
grid on

print(gcf, '-dpdf', '../plots/Encoding_Results');

%% performance split by category 

mACC = [behav_out.encSummary.FaceHR  behav_out.encSummary.SceneHR];
RTs  = [behav_out.encSummary.meanFaceRTs behav_out.encSummary.meanSceneRTs];

figure(2);clf;
set(gcf,'position',[100 100 600 400],'paperpositionmode','auto','color','w',...
'paperposition',[0.2 0.2 0.6 0.4],'paperunits','normalized');

a1 = axes('position',[0.1 0.1 0.3 0.8]); hold on;
for ss = 1:nSubjs
    h1 = plot([1 2],mACC(ss,:),'o-','linewidth',2,'color',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k');
end
set(a1,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'faces','scenes'})
set(a1,'ytick',[0.5:0.1:1])
ylabel(' mean Accuracy ')
plot([0.8 1.2], [1 1]*mean(mACC(:,1)),'linewidth',5,'color','k')
plot([1.8 2.2], [1 1]*mean(mACC(:,2)),'linewidth',5,'color','k')

xlim([0.6 2.4])
grid on

a1 = axes('position',[0.6 0.1 0.3 0.8]); hold on;
for ss = 1:nSubjs
    h1 = plot([1 2],RTs(ss,:),'o-','linewidth',2,'color',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k');
end
set(a1,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'faces','scenes'})
set(a1,'ytick',[0.5:0.1:1])
ylabel(' mean RTs (s) ')
plot([0.8 1.2], [1 1]*mean(RTs(:,1)),'linewidth',5,'color','k')
plot([1.8 2.2], [1 1]*mean(RTs(:,2)),'linewidth',5,'color','k')

xlim([0.6 2.4])
grid on
print(gcf, '-dpdf', '../plots/Encoding_ResultsFaceVsScene');

%% retrieval

dPrime = behav_out.retSummary.dPrime;
nSubjs = numel(dPrime);

figure(1);clf;
set(gcf,'position',[100 100 800 400],'paperpositionmode','auto','color','w',...
'paperposition',[0.1 0.2 0.8 0.4],'paperunits','normalized');

a1 = axes('position',[0.1 0.1 0.1 0.8]); hold on;
h1 = plot(ones(nSubjs,1)+0.05*randn(nSubjs,1),dPrime,'o');
set(h1,'markerfacecolor',[0.5 0.5 0.5],'markeredgeColor','k','markersize',10)
set(a1,'fontSize',16,'lineWidth',3,'xtick',[])
set(a1,'ytick',[0:0.2:1])
ylim([0 1])
ylabel(' d prime ')

M_dPrime = mean(dPrime);

plot([0.9 1.1],[M_dPrime M_dPrime],'linewidth',4,'color','k')
xlim([0.8 1.2])

grid on

a2 = axes('position',[0.35 0.1 0.2 0.8]); hold on;
accByConf = behav_out.retSummary.meanAccuracyByConf;

for ss = 1:nSubjs
    h1 = plot([1 2 3],[accByConf(ss,:)],'o-','linewidth',2,'color',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
end
set(a2,'fontSize',16,'lineWidth',3,'xtick',[1 2 3],'xticklabel',{'low','mid','high' })
set(a2,'ytick',[0:0.2:1])
ylim([0 1])
xlim([0.6 3.4])
plot([0.8 1.2], [1 1]*mean(accByConf(:,1)),'linewidth',5,'color','k')
plot([1.8 2.2], [1 1]*mean(accByConf(:,2)),'linewidth',5,'color','k')
plot([2.8 3.2], [1 1]*mean(accByConf(:,3)),'linewidth',5,'color','k')
ylabel(' acc by confidence ')
grid on

dPrimeFace = behav_out.retSummary.Face_dPrime;
dPrimeScene = behav_out.retSummary.Scene_dPrime;

a3 = axes('position',[0.7 0.1 0.2 0.8]); hold on;
for ss = 1:nSubjs
    h1 = plot([1 2],[dPrimeFace(ss) dPrimeScene(ss)],'o-','linewidth',2,'color',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
end
set(a3,'fontSize',16,'lineWidth',3,'xtick',[1 2],'xticklabel',{'faces','scenes'})
set(a3,'ytick',[0:0.4:1.6])
ylim([0 1.6])
ylabel(' d prime ')
plot([0.8 1.2], [1 1]*mean(dPrimeFace),'linewidth',5,'color','k')
plot([1.8 2.2], [1 1]*mean(dPrimeScene),'linewidth',5,'color','k')

xlim([0.6 2.4])
grid on

print(gcf, '-dpdf', '../plots/RetBehavPerf_Results');
%%
RTs = [ behav_out.retSummary.meanHit_RTs' behav_out.retSummary.meanCRs_RTs' ...
    behav_out.retSummary.meanMiss_RTs' behav_out.retSummary.meanFA_RTs'];


figure(2);clf;
set(gcf,'position',[100 100 800 400],'paperpositionmode','auto','color','w',...
'paperposition',[0.1 0.2 0.8 0.4],'paperunits','normalized');

a1 = axes('position',[0.15 0.1 0.4 0.8]); hold on;

for ss = 1:nSubjs
    h1 = plot([1 2 3 4]+0.05*randn(1,4),[RTs(ss,:)],'o','linewidth',2,'color',[0.5 0.5 0.5],...
        'MarkerFaceColor',[0.5 0.5 0.5],'markeredgecolor','k','MarkerSize',10);
end
set(a1,'fontSize',16,'lineWidth',3,'xtick',[1 2 3 4],'xticklabel',{'hits','CRs','Misses','FAs'})
set(a1,'ytick',[0:0.5:2.5])
%ylim([0 1])
xlim([0.6 4.4])
plot([0.8 1.2], [1 1]*mean(RTs(:,1)),'linewidth',5,'color','k')
plot([1.8 2.2], [1 1]*mean(RTs(:,2)),'linewidth',5,'color','k')
plot([2.8 3.2], [1 1]*mean(RTs(:,3)),'linewidth',5,'color','k')
plot([3.8 4.2], [1 1]*mean(RTs(:,4)),'linewidth',5,'color','k')
ylabel(' Retrieval RTs (s) ')
grid on

print(gcf, '-dpdf', '../plots/RetBehavRTs_Results');


