close all
clearvars

dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
load([dataPath 'Summary/BehavSummary.mat'])
load([dataPath 'Summary/PhaseDependentAnalyses.mat'])

uniqSeqSelection = 2;
if uniqSeqSelection==1
    subjs  = find(out.SubjWithUniqueSeq);
    SubjSelectStr = 'SSuniqueSeq';
elseif uniqSeqSelection==2
    subjs1  = find(out.SubjWithUniqueSeq);
    subjs2  = find(behav_out.encSummary.goodSubj);
    subjs   = intersect(subjs1,subjs2);
    SubjSelectStr = 'SSuniqueSeq_ValidEnc';
else
    subjs = 1:numel(behav_out.retSubj);
    SubjSelectStr = 'all';
end
nSubjs  = numel(subjs);
PhasesDeg = 36:72:359;
PhasesRad = PhasesDeg./180*pi;

%% Encoding Results: Categorization Task
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
HR = [ behav_out.encSummary.meanAcc(subjs) behav_out.encSummary.FaceHR(subjs) behav_out.encSummary.SceneHR(subjs)]*100;
Dstrs = {'ACC','Face','Scn'};
disp(array2table(mean(HR),'variablenames',Dstrs))
[~,p,~,t]=ttest(HR(:,2),HR(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn ACC'},'variablenames',{'T','P'}))

% Figure : (a) ACC
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = HR(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y); 
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim([90 100])
xlim([0 1])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)

% Figure 1: (b) acc by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = HR(:,2:3);
% 
% for ii =1:nSubjs
%     plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])   
% end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Faces
s=scatter(x,y(:,1)); 
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [100 200 100]/255;
s.MarkerFaceColor = [100 200 100]/255;

% Scenes
s=scatter(x+1,y(:,2)); 
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([90 100])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[90 100])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerf_' SubjSelectStr])

%%  Categorization by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
HR = behav_out.encSummary.HRByPhase(subjs,:,:)*100;

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 200 400 150]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = squeeze(HR(:,1,:));
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*mean(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [119,136,153]/255;
    s.MarkerFaceColor = [119,136,153]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([85 100])
xlim([0 5])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)


a2 = axes('units','points','position',[100 130 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhase_' SubjSelectStr])

%% Categorization By Phase / Stim Type
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(2); % for location reproducibility
HR = behav_out.encSummary.HRByPhase(subjs,:,:)*100;

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 200 400 150]); hold on;
x1 = randn(nSubjs,1)*0.02+0.3;
x2 = randn(nSubjs,1)*0.02+0.7;
y1 = squeeze(HR(:,2,:));
y2 = squeeze(HR(:,3,:));
nBars = 5;
for ii =1:nBars
    xx1 = x1+(ii-1);
    xx2 = x2+(ii-1);
    yy1 = y1(:,ii);
    yy2 = y2(:,ii);
    
    % Faces
    plot([0.15 0.45]+(ii-1), ones(1,2)*mean(yy1),'linewidth',4,'color','k')    
    s = scatter(xx1,yy1); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [100 200 100]/255;
    s.MarkerFaceColor = [100 200 100]/255;
    
    % Scenes
    plot([0.55 0.85]+(ii-1), ones(1,2)*mean(yy2),'linewidth',4,'color','k')
    s = scatter(xx2,yy2); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [200 100 200]/255;
    s.MarkerFaceColor = [200 100 200]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([85 100])
xlim([0 5])
set(gca,'ytick',[90 100])
ylabel(' ACC (%) ')
set(gca,'LineWidth',2)


a2 = axes('units','points','position',[100 130 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationPerfByPhaseStimType_' SubjSelectStr])

%% Categorization RTs:
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
RTs = [behav_out.encSummary.medianRTs(subjs)  ... 
    behav_out.encSummary.medianFaceRTs(subjs) behav_out.encSummary.medianSceneRTs(subjs)];
Dstrs = {'RTs','Face','Scn'};
disp(array2table(mean(RTs),'variablenames',Dstrs))
[~,p,~,t]=ttest(RTs(:,2),RTs(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn RTs'},'variablenames',{'T','P'}))

% Figure : (a) ACC
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = RTs(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y); 
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim([0.5 1])
xlim([0 1])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)

%
% Figure 1: (b) acc by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = RTs(:,2:3);

for ii =1:nSubjs
    if y(ii,1)>y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.8 0.2 0.2])
    elseif y(ii,1)<=y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])
    end
end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Faces
s=scatter(x,y(:,1)); 
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [100 200 100]/255;
s.MarkerFaceColor = [100 200 100]/255;

% Scenes
s=scatter(x+1,y(:,2)); 
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([0.5 1])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[0.5:0.2:1])
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTs_' SubjSelectStr])

%% Categorization RTs by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(1); % for location reproducibility
RTs = squeeze(behav_out.encSummary.medianRTsByPhase(subjs,1,:));

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 200 400 150]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = RTs;
nBars = 5;
for ii =1:nBars
    xx = x+(ii-1);
    yy = y(:,ii);
    plot([0.2 0.8]+(ii-1), ones(1,2)*mean(yy),'linewidth',4,'color','k')
    s = scatter(xx,yy); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [119,136,153]/255;
    s.MarkerFaceColor = [119,136,153]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([0.5 1])
xlim([0 5])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)


a2 = axes('units','points','position',[100 130 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTsByPhase_' SubjSelectStr])

%% Categorization RTs by Phase and category
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr PhasesDeg
rng(2); % for location reproducibility
RTs = behav_out.encSummary.medianRTsByPhase(subjs,:,:);

% Figure : ACC by Phase
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])
a = axes('units','points','position',[100 200 400 150]); hold on;
x1 = randn(nSubjs,1)*0.02+0.3;
x2 = randn(nSubjs,1)*0.02+0.7;
y1 = squeeze(RTs(:,2,:));
y2 = squeeze(RTs(:,3,:));
nBars = 5;
for ii =1:nBars
    xx1 = x1+(ii-1);
    xx2 = x2+(ii-1);
    yy1 = y1(:,ii);
    yy2 = y2(:,ii);
    
    for jj =1:nSubjs
        if yy1(jj)>yy2(jj)
            plot([xx1(ii) xx2(ii)], [yy1(jj) yy2(jj)],'-','color',[0.8 0.2 0.2])
        else
            plot([xx1(ii) xx2(ii)], [yy1(jj) yy2(jj)],'-','color',[0.6 0.6 0.6])
        end
    end

    % Faces
    plot([0.15 0.45]+(ii-1), ones(1,2)*mean(yy1),'linewidth',4,'color','k')    
    s = scatter(xx1,yy1); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [100 200 100]/255;
    s.MarkerFaceColor = [100 200 100]/255;
    
    % Scenes
    plot([0.55 0.85]+(ii-1), ones(1,2)*mean(yy2),'linewidth',4,'color','k')
    s = scatter(xx2,yy2); 
    s.MarkerFaceAlpha=0.5;
    s.MarkerEdgeAlpha=0.4;
    s.SizeData          = 120;
    s.MarkerEdgeColor = [200 100 200]/255;
    s.MarkerFaceColor = [200 100 200]/255;
end
set(gca,'fontsize',20,'xtick',[0.5:5],'xticklabel',[])%,'xtick',[0.5:5],'xticklabel',PhasesDeg)
ylim([0.5 1])
xlim([0 5])
set(gca,'ytick',[0.5:0.2:1])
ylabel(' RTs (s) ')
set(gca,'LineWidth',2)


a2 = axes('units','points','position',[100 130 400 50]); hold on;
xa = linspace(0,2*pi,1000); x = cos(xa);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight
set(gca,'ytick',[],'ycolor','w','fontsize',20,'box','off','lineWidth',2)
set(gca,'xtick',PhasesDeg)
xlabel(' Encoding Phase (deg)')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'CategorizationRTsByPhaseStimType_' SubjSelectStr])

%% dPrimes
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
rng(1); % for location reproducibility

Dstrs   = {'dPrime','Face_dPrime','Scene_dPrime',...
    'dPrime_C','Face_dPrime_C','Scene_dPrime_C'};

Dstrs2  = {'dP', 'Face dP', 'Scn dP', 'C', 'Face C', 'Scn C'};
Dstrs3  = {'dP', 'FaceDP', 'ScnDP', 'C', 'FaceC', 'ScnC'};
nD = numel(Dstrs);
D       = zeros(nSubjs,nD);
for ii = 1:nD
    D(:,ii) = behav_out.retSummary.(Dstrs{ii})(subjs);
end
disp(array2table(mean(D),'variablenames',Dstrs3))
[~,p,~,t]=ttest(D(:,2),D(:,3));
disp(table(t.tstat,p,'rownames',{'Face vs Scn DP'}))

yLims = [-0.1 1.5; 0 2; 0 2; -1.2 1.2; -1.2 1.2; -1.2 1.2];
yTicks1 = [0 1 2];

% Figure 1: (a) d-prime
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[200,200,600,400])

a = axes('units','points','position',[100 100 80 200]); hold on;
x = randn(nSubjs,1)*0.05+0.5;
y = D(:,1);
plot([0.2 0.8], ones(1,2)*mean(y),'linewidth',4,'color','k')
s=scatter(x,y); 

s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [119,136,153]/255;
s.MarkerFaceColor = [119,136,153]/255;

set(gca,'fontsize',20,'xtick',[])
ylim(yLims(1,:))
xlim([0 1])
set(gca,'ytick',yTicks1)
ylabel('dP')
set(gca,'LineWidth',2)

% Figure 1: (b) d-prime by category
a = axes('units','points','position',[300 100 180 200]); hold on;
y = D(:,2:3);

for ii =1:nSubjs
    if y(ii,1)>y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.8 0.2 0.2])
    elseif y(ii,1)<=y(ii,2)
        plot([x(ii) x(ii)+1], y(ii,:),'-','color',[0.6 0.6 0.6])
    end
end

plot([0.2 0.8], ones(1,2)*mean(y(:,1)),'linewidth',4,'color','k')
plot([0.2 0.8]+1, ones(1,2)*mean(y(:,2)),'linewidth',4,'color','k')

% Faces
s=scatter(x,y(:,1)); 
s.MarkerFaceAlpha=0.5;
s.MarkerEdgeAlpha=0.4;
s.SizeData          = 120;
s.MarkerEdgeColor = [100 200 100]/255;
s.MarkerFaceColor = [100 200 100]/255;

% Scenes
s=scatter(x+1,y(:,2)); 
s.MarkerFaceAlpha   = 0.5;
s.MarkerEdgeAlpha   = 0.4;
s.SizeData          = 120;
s.MarkerEdgeColor   = [200 100 200]/255;
s.MarkerFaceColor   = [200 100 200]/255;

xlim([-0.1 2.1])
ylim([-0.1 2.1])
set(gca,'fontsize',20,'xtick',[0.5 1.5],'xticklabel',{'F','S'})
set(gca,'ytick',[0 1 2])
ylabel(' dP ')
set(gca,'LineWidth',2)
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'dPrimes' SubjSelectStr])

%% Confidence
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% dPrime by confidence:
DPC = behav_out.retSummary.dPrimeConf(subjs,:);
figure(2); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,600,400])

hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPC(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
end
plot(1:3,mean(DPC),'linewidth',3,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',2,'ytick',[0:3])
ylabel(' dP ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/'  'Confidence_dPrime' SubjSelectStr])


disp(array2table(mean(DPC),'variablenames',{'Low','Med','High'}))
[~,p1,~,t1]=ttest(DPC(:,2),DPC(:,1));
[~,p2,~,t2]=ttest(DPC(:,3),DPC(:,2));
[~,p3,~,t3]=ttest(DPC(:,3),DPC(:,1));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'Mid>Low','Hi>Mid','Hi>Lo'}))

anova1(DPC(:),[ones(nSubjs,1);2*ones(nSubjs,1);3*ones(nSubjs,1)])

%% Confidence Face/Scenes
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
close all
% Faces
DPCF = behav_out.retSummary.Face_dPrimeConf(subjs,:);
figure(3); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,600,400])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPCF(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
end
plot(1:3,nanmean(DPCF),'linewidth',3,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',2,'ytick',[0:3])
ylabel(' dP ')

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'Confidence_dPrimeFaces' SubjSelectStr])
disp('Face dPrime by Confidence ')
disp(array2table(nanmean(DPCF),'variablenames',{'Low','Med','High'}))

% Scenes
DPCS = behav_out.retSummary.Scene_dPrimeConf(subjs,:);
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,600,400])
hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPCS(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);
end
plot(1:3,nanmean(DPCS),'linewidth',3,'color',[0.1 0.1 0.1]);
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'})
xlim([0.8 3.2])
ylim([-0.5 3.5])
set(gca,'LineWidth',2,'ytick',[0:3])
ylabel(' dP ')
print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'Confidence_dPrimeScn' SubjSelectStr])

disp('Scene dPrime by Confidence ')
disp(array2table(nanmean(DPCS),'variablenames',{'Low','Med','High'}))

%% RTs
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

strs    = {'medianHit_RTs','medianMiss_RTs','medianFA_RTs','medianCRs_RTs'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTs       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTs(:,ii) = behav_out.retSummary.(strs{ii})(subjs);
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 500 300])
set(gcf,'position',[100,100,500,300])

dx = 80;
xPosCore = [100:dx:500];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 50 60 200]);
end
yLims = [0.5 2];
yTicks1 = [0:0.5:2.5];
for ii=1:nRTConds
    axes(a(ii))
    %set(gca,'Position',[xPos(ii) 100 80 250])
    y = RTs(:,ii);
    
    for ss=1:nSubjs
        s = scatter(ss,y(ss)); hold on;
        s.MarkerFaceAlpha=0.8;
        s.MarkerEdgeAlpha=0.6;
        s.SizeData          = 50;
        s.MarkerEdgeColor = [0.45 [0.75,0.9]*ss/nSubjs];
        s.MarkerFaceColor = [0.45 [0.75,0.9]*ss/nSubjs];
    end
    
    set(gca,'fontsize',20,'xTick','')
    plot([1 nSubjs], ones(1,2)*mean(y),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0 nSubjs+1])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' RTs (s) ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'RTs' SubjSelectStr])

disp(array2table(mean(RTs),'variablenames',strs2))
[~,p1,~,t1]=ttest(RTs(:,1),RTs(:,2));
[~,p2,~,t2]=ttest(RTs(:,1),RTs(:,3));
[~,p3,~,t3]=ttest(RTs(:,1),RTs(:,4));
p = [p1;p2;p3];
t = [t1.tstat;t2.tstat;t3.tstat];
disp(table(t,p,'rownames',{'HvsM','HvFA','HvCRs'}))

%% RTs by Confidence
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
strs    = {'medianHit_RTsConf','medianMiss_RTsConf','medianFA_RTsConf','medianCRs_RTsConf'};
strs2   = {'Hits','Misses','FA','CRs'};
nRTConds = numel(strs);

RTsConf = zeros(nSubjs,3,nRTConds);
for ii = 1:nRTConds
    RTsConf(:,:,ii) = behav_out.retSummary.(strs{ii})(subjs,:);
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,150,500,250])

dx = 90;
xPosCore = [80:dx:600];
xPos = [xPosCore ];
a     = zeros(nRTConds,1);
for ii =1:nRTConds
    a(ii)=axes('units','points','position',[xPos(ii) 50 80 180]);
end
yLims = [0.5 2.5];
yTicks1 = [0:0.5:2.5];
%
for ii=1:nRTConds
    axes(a(ii));
    hold on;
    %set(gca,'Position',[xPos(ii) 100 80 250])
    Y = RTsConf(:,:,ii);
    
    for ss=1:nSubjs
        p = plot(1:3,Y(ss,:));
        p.Color = [0.7 0.7 0.7];
        p.LineWidth = 1;
    end
    set(gca,'fontsize',20,'xTick',1:3, 'xticklabel',{'L','M','H'} )
    plot(1:3,nanmean(Y),'linewidth',4,'color','k')
    xlabel(strs2{ii})
    ylim(yLims)
    xlim([0.8 3.2])
    if ii==1
        set(gca,'ytick',yTicks1)
        ylabel(' RTs (s) ')
    end
    if ii>1
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end

print(gcf,'-dpdf',['../plots/tacs_enc_xdiva/Behavior/' 'RTsByConf' SubjSelectStr])

%% PropHits PropMisses and Difference by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase(subjs,:);
plot([36:72:360],X','-','color',[255 180 150]/255)
plot([36:72:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (hits) ' )

print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsPropByPhase'  SubjSelectStr]);

% Misses
figure(5); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propMissByPhase(subjs,:);
plot([36:72:360],X','-','color',[150 220 220]/255)
plot([36:72:360],mean(X), 'color',[120 200 200]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (misses) ' )

print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MissPropByPhase'  SubjSelectStr]);

% Difference
figure(6); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase-out.propMissByPhase;
X = X(subjs,:);
plot([36:72:360],X','-','color',[180 180 180]/255)
plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([-0.15 0.15])
ylabel(' p(h)-p(m) ' )

print(gcf, '-dpdf',['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'H-MPropByPhase' SubjSelectStr]);

%% Hit Rate by phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

X = out.HR_ByPhase(subjs,:);
xa = linspace(0,2*pi,1000);
x = cos(xa);

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa./pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
plot([36:72:360],X','-','color',[180 180 180]/255)
plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([-1 361])
ylim([0 1])
ylabel(' Hit Rate ' )
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsRatepByPhase' SubjSelectStr]);

% mean Vectors
Z=X.*exp(1j*repmat([36:72:360]./180*pi,[nSubjs,1]));
mZ = mean(Z,2);
th = mod(angle(mZ),2*pi); rho = abs(mZ);
opts = [];
opts.colors = [180 180 180]/255;
opts.markerSize=behav_out.retSummary.dPrime(subjs)*150;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitsRateMeanVec' SubjSelectStr]);
disp(table(circ_rtest(th),'variablenames',{'rhao_test_pval'}))

% scatter of radii to dPrime
x = rho;
y = behav_out.retSummary.dPrime(subjs);
opts =[];
opts.colors = [150 150 150]/255;
opts.ylabel = ' dPrime ';
opts.xlabel = ' \rho ';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,y,'type','spearman')*100)/100)];
opts.xytext     = [0.04 1.1];
han = xyScatter(x,y,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'HitRateScatterRho-Dprime' SubjSelectStr]);

%% MemScore by Phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

xa = linspace(0,2*pi,1000);
x = cos(xa);

X = out.MemScoreByPhase(subjs,:);
figure(7); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot(xa/pi*180,x,'k','linewidth',4)
axis tight;
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360])
xlabel(' Encoding Phase (deg)')
xlim([-1 361])
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
Xm=(X-repmat(mean(X,2),[1,5]));
plot([36:72:360],Xm','-','color',[180 180 180]/255)
plot([36:72:360],mean(Xm), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([-1 361])
ylim([-1 1])
ylabel(' MemScore ' )
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MemScoreByPhase' SubjSelectStr] );

% polar plot
Z=X.*exp(1j*repmat([36:72:360]./180*pi,[nSubjs,1]));
mZ = mean(Z,2);
th = mod(angle(mZ),2*pi); rho = abs(mZ);
opts = [];
opts.colors = [180 180 180]/255;
opts.markerSize=behav_out.retSummary.dPrime(subjs)*150;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/'  'MemScoreMeanVec' SubjSelectStr]);
disp(table(circ_rtest(th),'variablenames',{'rhao_test_pval'}))

% scatter of radii to dPrime
x = rho;
y = behav_out.retSummary.dPrime(subjs);
opts =[];
opts.colors = [150 150 150]/255;
opts.ylabel = ' dPrime ';
opts.xlabel = ' \rho ';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,y,'type','spearman')*100)/100)];
opts.xytext     = [0.04 1.1];
han = xyScatter(x,y,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MemScoreRho-Dprime' SubjSelectStr]);

%% Hit Miss Phase distributions
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

opts = [];
opts.colors = [255 180 150; 150 220 220]/255;
opts.markerSize=squeeze(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,:,1:2),2));

% raw vectors
th  = mod(out.HitMissMePhase(subjs,:),2*pi);
rho = out.HitMissAbPhase(subjs,:);
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond' SubjSelectStr]);

%test and disply circular uniformity
disp('Rayleigh Test:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(th(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% using r-statistic as magnitude.
rhoR     = out.HitMiss_DistFromUniform(subjs,:);
han = PolarPlot(th,rhoR,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCondRStat' SubjSelectStr]);

% raw only theta (ignores individuals SS strength)
opts.maxR = 4/3;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/''MeanPhaseVectsByCond-Theta' SubjSelectStr]);

%% line between points
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.HitMissMePhase(subjs,:),2*pi);
rho = out.HitMissAbPhase(subjs,:);
z   = rho.*exp(1i*th);

opts = [];
opts.markerSize =squeeze(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,:,1:2),2));
opts.colors = [255 180 150;150 220 220]/255;
opts.connect = 1; opts.polarGrid =0; opts.meanVecs = 0;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-Connected'  SubjSelectStr]);

% line between hit-miss angles without magnitude
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid=0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ThetaConnected'  SubjSelectStr]);

% figure with angles re-centered
thc  = th-repmat(th(:,1),[1 2]);
opts.polarGrid = 1;
han = PolarPlot(thc,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ConnectedReCenter'  SubjSelectStr]);

% mean Vecs re-centered
opts.meanVecs = 1;
opts.connect  = 0;
han = PolarPlot(thc,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectsByCond-ReCenter'  SubjSelectStr]);

%% difference in hit phase and miss phase
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr
% mean differences
th = mod(out.HitMissMePhase(subjs,:),2*pi);
rho = out.HitMissAbPhase(subjs,:);
dTh= th(:,1)-th(:,2);

z = rho.*exp(1j*th);
R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiff' SubjSelectStr]);

% difference in hit phase and miss phase Rstat
rhoR = out.HitMiss_DistFromUniform(subjs,:);
z = rhoR.*exp(1j*th); RR = abs(z(:,1)-z(:,2));
han = PolarPlot(dTh,RR,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiff-Rstat' SubjSelectStr]);

% only Delta in Theta
opts.maxR = 4/3;
han = PolarPlot(dTh,ones(nSubjs,1),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'MeanPhaseVectDiffTheta' SubjSelectStr]);

% variance decreases as mean vector length increases
%d = 0.5*(1-cos(dTh));
x = R;
d = abs(pi-abs(dTh));
opts =[];
opts.colors = [150 150 150]/255;
opts.xlabel = ' \rho ';
opts.ylabel = ' |\pi-|\Delta \theta| |';
opts.polyfitN = 1;
opts.text       =['R = ' num2str(round(corr(x,d,'type','spearman')*100)/100)];
opts.xytext     = [0.1 0.5];
han = xyScatter(x,d,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'DiffAbAllHitMissesScatter' SubjSelectStr]);

% variance decreases as mean vector length increases R-stat
%d = 0.5*(1-cos(th(:,1)-th(:,2)));
x = RR;
opts.text       =['R = ' num2str(round(corr(x,d,'type','spearman')*100)/100)];
opts.xytext     = [1 0.5];
han = xyScatter(x,d,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_HitMiss/' 'Diff-RstatAllHitMissesScatter' SubjSelectStr]);

[p,u]=circ_rtest(th(:,1)-th(:,2));
disp( 'Rayleigh test for Hit-Miss Angle: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% Face / Scene differences independent of memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th      = mod(out.FaScnMePhases(subjs,:),2*pi);
rho     = out.FaScnAbPhases(subjs,:);
opts            = [];
opts.colors     = [100 200 100; 200 100 200]/255;
opts.markerSize = out.FaScnHitN(subjs,:)+out.FaScnMissN(subjs,:);
opts.maxR       = 1/100;
han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaceScn' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(th(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Faces','Scenes'}))

% difference
dTh  = th(:,1)-th(:,2);
z = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaceScnDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Face-Scene Angle: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Face-Scene'}))

%% Faces difference by memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% Hit/Miss for Faces
thF = mod([out.FaScnHitMePhases(subjs,1) , out.FaScnMissMePhases(subjs,1)],2*pi);
rhoF= [out.FaScnHitAbPhases(subjs,1), out.FaScnMissAbPhases(subjs,1)];
dThF = mod(thF(:,1)-thF(:,2),2*pi);
zF = rhoF.*exp(1j*thF); RF = abs(zF(:,1)-zF(:,2));

opts        = [];
opts.colors = [50 100 50;100 200 100]/255;
opts.markerSize = [out.FaScnHitN(subjs,1) out.FaScnMissN(subjs,1)];
han = PolarPlot(thF,rhoF,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMiss' SubjSelectStr]);

% hit/miss for faces projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(thF,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissTheta' SubjSelectStr]);

% hit/miss for faces lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(thF,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissThetaConnected' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test Faces:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(thF(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% Hit Miss Difference for Faces
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dThF,RF,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(thF(:,1)-thF(:,2));
disp( 'Rayleigh test for Hit-Miss Angle Faces: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% Scenes difference by memory
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

% Hit/Miss for Scenes
thS = mod([out.FaScnHitMePhases(subjs,2) , out.FaScnMissMePhases(subjs,2)],2*pi);
rhoS= [out.FaScnHitAbPhases(subjs,2), out.FaScnMissAbPhases(subjs,2)];
dThS = mod(thS(:,1)-thS(:,2),2*pi);
zS = rhoS.*exp(1j*thS); RS = abs(zS(:,1)-zS(:,2));

opts        = [];
opts.colors = [100 50 100;200 100 200]/255;
opts.markerSize = [out.FaScnHitN(subjs,2) out.FaScnMissN(subjs,2)];
han = PolarPlot(thS,rhoS,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMiss' SubjSelectStr]);

% hit/miss for scenes projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(thS,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMissTheta' SubjSelectStr]);

% hit/miss for scenes lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(thS,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFacesHitMissThetaConnected' SubjSelectStr]);

%tests of circular uniformity
disp('Rayleigh Test Scenes:')
p = zeros(2,1); u = p;
for ii = 1:2
    [p(ii),u(ii)]=circ_rtest(thS(:,ii));
end
disp(table(p,u,'VariableNames',{'P_Val','Z'},'rownames',{'Hits','Misses'}))

% Hit Miss Difference for Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dThS,RS,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsScnHitMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(dThS);
disp( 'Rayleigh test for Hit-Miss Scenes: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Hit-Miss'}))

%% Face / Scene Hits
close all
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.FaScnHitMePhases(subjs,:),2*pi);
rho = out.FaScnHitAbPhases(subjs,:);
dTh = mod(th(:,1)-th(:,2),2*pi);
z   = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [50 100 50; 100 50 100]/255;
opts.markerSize = out.FaScnHitN(subjs,:);

han = PolarPlot(th,rho,opts);
print(han, '-dpdf',['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHits' SubjSelectStr]);

% hits projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsTheta' SubjSelectStr]);

% hits for connected lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsThetaConnected' SubjSelectStr]);

% Hits Difference for Faces and Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnHitsDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Hits Angle Faces-Scn: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Fa-Scn'}))

%% Face / Scene Misses
close all
clearvars -except out behav_out subjs nSubjs dataPath SubjSelectStr

th  = mod(out.FaScnMissMePhases(subjs,:),2*pi);
rho = out.FaScnMissAbPhases(subjs,:);
dTh = mod(th(:,1)-th(:,2),2*pi);
z   = rho.*exp(1j*th); R = abs(z(:,1)-z(:,2));

opts = [];
opts.colors = [100 200 100; 200 100 200]/255;
opts.markerSize = out.FaScnMissN(subjs,:);

han = PolarPlot(th,rho,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMiss' SubjSelectStr]);

% hits projected to the unit circle
opts.maxR = 4/3;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissTheta' SubjSelectStr]);

% hits for connected lines
opts.maxR = 4/3;
opts.connect = 1; opts.meanVecs = 0; opts.polarGrid =0;
han = PolarPlot(th,ones(nSubjs,2),opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissThetaConnected' SubjSelectStr]);

% Hits Difference for Faces and Scenes
opts        = []; opts.colors = [150 150 150]/255;
han = PolarPlot(dTh,R,opts);
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Phase_FaceScene/' 'MeanPhaseVectsFaScnMissDiff' SubjSelectStr]);

[p,u]=circ_rtest(dTh);
disp( 'Rayleigh test for Miss Angle Faces-Scn: ')
disp(table(p,u,'variablenames',{'P_Val','Z'},'rownames',{'Fa-Scn'}))

%%
%
% %% hit and miss phases, high conf
% th = mod(out.HitMissMePhaseConf(subjs,:,3),2*pi);
% rho = out.HitMissAbPhaseConf(subjs,:,3);
% nh = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,3,1);
% nm = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,3,2);
% markerSize = [nh nm];
% colors = [255 180 150; 150 220 220]/255;
% % don't compute if nh/nm lower than 10 trials.
%
% ii= nh<10 | nm<10;
% th(ii,:) = []; rho(ii,:) = []; markerSize(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecHitMissHighConf' SubjSelectStr]);
%
% th = mod(out.HitMissMePhaseConf(subjs,:,3),2*pi);
% rho = out.HitMiss_DistFromUniformConf(subjs,:,3);
% th(ii,:) = []; rho(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecHitMissHighConfRStat' SubjSelectStr]);
%
% th = mod(out.HitMissMePhaseConf(subjs,:,3),2*pi);
% rho= out.HitMiss_DistFromUniformConf(subjs,:,3);
% th(ii,:) = []; rho(ii,:) = [];
% theta = th(:,1)-th(:,2);
% z = rho.*exp(1j*th);
% rho = abs(z(:,1)-z(:,2));
% colors = [150 150 150]/255;
%
% han = PolarPlot(theta,rho,colors,geomean(markerSize,2));
% print(han, '-dpdf', ['../plots/xdiva/MeanVecBySubj_HConfDist' SubjSelectStr]);
%
% p=circ_rtest(th(:,1)-th(:,2));
% disp( 'Difference in rtest: ')
% disp(table(p,'variablenames',{'rhao_test_pval'}))
%
% %% hit and miss phases, low conf
% th = mod(out.HitMissMePhaseConf(subjs,:,1),2*pi);
% rho = out.HitMissAbPhaseConf(subjs,:,1);
% nh = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1,1);
% nm = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1,2);
% markerSize=[nh nm];
% colors = [255 180 150; 150 220 220]/255;
% ii= nh<10 | nm<10;
% th(ii,:) = []; rho(ii,:) = []; markerSize(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf',[ '../plots/xdiva/MeanVecHitMissLoConf' SubjSelectStr]);
%
% %rstat
% th = mod(out.HitMissMePhaseConf(subjs,:,1),2*pi);
% rho = out.HitMiss_DistFromUniformConf(subjs,:,1);
% nh = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1,1);
% nm = behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1,2);
% markerSize=[nh nm];
% ii= nh<10 | nm<10;
% th(ii,:) = []; rho(ii,:) = []; markerSize(ii,:) = [];
% colors = [255 180 150; 150 220 220]/255;
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecHitMissLoConf' SubjSelectStr]);
%
% % low condifence difference
% th = mod(out.HitMissMePhaseConf(subjs,:,1),2*pi);
% rho = out.HitMiss_DistFromUniformConf(subjs,:,1);
% th(ii,:) = []; rho(ii,:) = [];
% thetaLC = th(:,1)-th(:,2);
% z = rho.*exp(1j*th);
% rhoLC = abs(z(:,1)-z(:,2));
% colors = [150 150 150]/255;
%
% han = PolarPlot(thetaLC,rhoLC,colors,geomean(markerSize,2));
% print(han, '-dpdf',[ '../plots/xdiva/MeanVecBySubj_LoConfDist'  SubjSelectStr]);
%
% [p,k]=circ_rtest(th(:,1)-th(:,2));
% [p1,k1]=circ_rtest(th(:,1));
% [p2,k2]=circ_rtest(th(:,2));
%
% disp( 'Rtests for Low Confidence: ')
% disp(table([p;p1;p2],[k;k1;k2],'variablenames',{'rhao_test_pval','rhao_test_k'},'rownames',{'diff','hits','misses'}))
%
% %% high-med confidence
%
% th = mod(out.HitMissMePhaseHMedConf(subjs,:),2*pi);
% rho = out.HitMiss_DistFromUniformHMedConf(subjs,:);
% nh = sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,2:3,1),2);
% nm = sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,2:3,2),2);
% colors = [255 180 150; 150 220 220]/255;
% markerSize=[nh nm];
% jj= nh<10 | nm<10;
% th(jj,:) = []; rho(jj,:) = []; markerSize(jj,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecBySubj_HMConfRstat'  SubjSelectStr]);
%
% % difference
% th = mod(out.HitMissMePhaseHMedConf(subjs,:),2*pi);
% rho = out.HitMiss_DistFromUniformHMedConf(subjs,:);
% th(jj,:) = []; rho(jj,:) = [];
% thetaHMC = th(:,1)-th(:,2);
% z = rho.*exp(1j*th);
% rhoHMC = abs(z(:,1)-z(:,2));
% colors = [150 150 150]/255;
%
% han = PolarPlot(thetaHMC,rhoHMC,colors,geomean(markerSize,2));
%
% print(han, '-dpdf', ['../plots/xdiva/MeanVecBySubj_HMConfDif'  SubjSelectStr]);
% [p,k]=circ_rtest(th(:,1)-th(:,2));
% [p1,k1]=circ_rtest(th(:,1));
% [p2,k2]=circ_rtest(th(:,2));
% disp( 'Rtests for High-Med Confidence: ')
% disp(table([p;p1;p2],[k;k1;k2],'variablenames',{'rhao_test_pval','rhao_test_k'},'rownames',{'diff','hits','misses'}))
%
% %% confidence weighted
%
% % hit miss polar
% th = out.HitMissCWMePhase(subjs,:);
% rho = out.HitMiss_CWDist(subjs,:);
% nh = sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1:3,1),2);
% nm = sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1:3,2),2);
% colors = [255 180 150; 150 220 220]/255;
% markerSize=[nh nm];
% ii= nh<10 | nm<10;
% th(ii,:) = []; rho(ii,:) = []; markerSize(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecBySubj_CWe'  SubjSelectStr]);
%
% % difference polar
% theta = th(:,1)-th(:,2);
% z = rho.*exp(1j*th);
% rhoZ = abs(z(:,1)-z(:,2));
% colors = [150 150 150]/255;
% han = PolarPlot(theta,rhoZ,colors,geomean(markerSize,2));
% print(han, '-dpdf', ['../plots/xdiva/MeanVecBySubj_CWeDif'  SubjSelectStr]);
% p=circ_rtest(th(:,1)-th(:,2));
% disp( 'Difference in rtest: ')
% disp(table(p,'variablenames',{'rhao_test_pval'}))
%
% % scatter of angle vs rho
% d = abs(pi-abs(th(:,1)-th(:,2)));
% x = rhoZ;
% opts =[];
% opts.colors = [150 150 150]/255;
% opts.xlabel = ' \rho ';
% opts.ylabel = ' |\pi-|\Delta \theta| |';
% opts.polyfitN = 1;
% opts.text       =['R = ' num2str(round(corr(x,d,'type','spearman')*100)/100)];
% opts.xytext     = [2 1.5];
% opts.markerSize=D(:,1)*150;
% han = xyScatter(x,d,opts);
% print(han, '-dpdf', ['../plots/xdiva/CWeScatterAng-RhoR' SubjSelectStr]);
%
% p=circ_rtest(th(:,1)-th(:,2));
% disp( 'Difference in rtest: ')
% disp(table(p,'variablenames',{'rhao_test_pval'}))
%
% %% plots of confidence phase independent of hit/miss status
% th  = out.ConfByPhaseMeVec(subjs,:);
% rho = out.ConfByPhaseAbVec(subjs,:);
% n1 = sum(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,1,1:2),3),2);
% n2 = sum(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,2,1:2),3),2);
% n3 = sum(sum(behav_out.retSummary.nH_nMiss_nFA_nCRs(subjs,3,1:2),3),2);
%
% colors = [46 204 204; 250 220 150; 255 51 51]/255;
% markerSize = [n1 n2 n3];
% ii=n1<10 | n2<10 | n3<10;
% th(ii,:) = []; rho(ii,:) = []; markerSize(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecAbBySubj_Conf' SubjSelectStr]);
% p =zeros(3,1); k = zeros(3,1);
% for kk = 1:3
%     [p(kk),k(kk)]=circ_rtest(th(:,kk));
% end
% disp( 'Difference in rtest for Confidence: ')
% disp(table(p,k,'variablenames',{'rhao_test_pval','rhao_test_k'},'rownames',{'Lo','Med','High'}))
%
% for kk = 1:3
%     [p(kk),k(kk)]=circ_rtest(th(:,kk)-th(:,mod(kk+1,3)+1));
% end
% disp( 'Difference in rtest Confidence Differences: ')
% disp(table(p,k,'variablenames',{'rhao_test_pval','rhao_test_k'},'rownames',{'H-L','M-L','H-M'}))
%
% th  = out.ConfByPhaseMeVec(subjs,:);
% rho = out.ConfByPhaseRStatVec(subjs,:);
% colors = [46 204 204; 250 220 150; 255 51 51]/255;
% th(ii,:) = []; rho(ii,:) = [];
% han = PolarPlot(th,rho,colors,markerSize);
% print(han, '-dpdf', ['../plots/xdiva/MeanVecRstatBySubj_Conf' SubjSelectStr]);
%%
% %% relationship of bias to magnitude of mean vector:
% th = mod(out.HitMissMePhase,2*pi);
% rho = out.HitMissAbPhase;
%
% z = rho.*exp(1i*th);
% %zD = z(:,1)./z(:,2);
% zD = abs((rho(:,1)+rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));
%
% figure(11); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% x=rho(:,1)-rho(:,2);
% y=D(:,4);
%
% s=scatter(x,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [180 180 180]/255;
% s.MarkerFaceColor = [180 180 180]/255;
%
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% xlabel(' \rho_{Hits}-\rho_{Misses} ' )
% ylabel(' Bias ')
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecRelBias');
% %%  relationship of d' to magnitude of mean vector:
% th = mod(out.HitMissMePhase,2*pi);
% rho = out.HitMissAbPhase;
%
% z = rho.*exp(1i*th);
% %zD = z(:,1)./z(:,2);
% zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));
%
% figure(11); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% x=abs(rho(:,1)+rho(:,2));
% y=D(:,1);
%
% s=scatter(x,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [180 180 180]/255;
% s.MarkerFaceColor = [180 180 180]/255;
%
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% xlabel(' |\rho_{Hits}+\rho_{Misses}| ' )
% ylabel(' dPrime ')
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecReldP')
%
% %%  relationship of d' to magnitude of mean vector for misses and hits
% th = mod(out.HitMissMePhase,2*pi);
% rho = out.HitMissAbPhase;
%
% z = rho.*exp(1i*th);
% %zD = z(:,1)./z(:,2);
% zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));
%
% figure(11); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% x1=rho(:,1);
% x2=rho(:,2);
% y=D(:,1);
%
% s=scatter(x1,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [255 180 150]/255;
% s.MarkerFaceColor = [255 180 150]/255;
%
% l=lsline;
% l(1).LineWidth= 5;
% l(1).Color=[240 120 100]/255;
%
% hold on
%
% s=scatter(x2,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [150 220 220]/255;
% s.MarkerFaceColor = [150 220 220]/255;
%
%
% l=lsline;
% l(1).LineWidth= 5;
% l(1).Color=[120 200 200]/255;
%
%
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% xlabel(' \rho ' )
% ylabel(' dPrime ')
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecHMReldP')
%
% %%
% %% MEeanVecDiff WEighted by confidence.
% th = mod(out.HitMissCWMePhase,2*pi);
% rho = out.HitMissCWAbPhase;
%
% z = rho.*exp(1i*th);
% %zD = z(:,1)./z(:,2);
% zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));
%
% figure(12); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% xLims = [-0.15 0.15];
% xlim(xLims);
% ylim(xLims); hold on;
% lCol = [0.8 0.8 0.8];
% plot(xlim,[0 0],'color',[0.8 0.8 0.8],'linewidth',2);
% plot([0 0],ylim,'color',lCol,'linewidth',2);
% plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
% plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);
%
% reZ = real(zD);
% imZ = imag(zD);
%
% s=scatter(reZ,imZ,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [180 180 180]/255;
% s.MarkerFaceColor = [180 180 180]/255;
%
% plot([0 mean(reZ)],[0 mean(imZ)],'linewidth',5,'color',[20 20 20]/255);
%
% axis off
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecCweBySubj3');
%
% %% Mean vector weighted by confidence relation to bias
% th = mod(out.HitMissCWMePhase,2*pi);
% rho = out.HitMissCWAbPhase;
%
% z = rho.*exp(1i*th);
% %zD = z(:,1)./z(:,2);
%
% figure(13); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% x=rho(:,1)-rho(:,2);
% y=D(:,4);
%
% s=scatter(x,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [180 180 180]/255;
% s.MarkerFaceColor = [180 180 180]/255;
%
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% xlabel(' log_{10}(|hit/miss|) ' )
% ylabel(' Bias ')
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecCWeRelBias');
%
% %% mean vectors (h/miss) weighted by conf relation to dP
% th = mod(out.HitMissCWMePhase,2*pi);
% rho = out.HitMissCWAbPhase;
%
% figure(11); clf;
% set(gcf,'paperpositionmode','auto','color','white')
% set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
% set(gcf,'position',[100,100,500,500])
%
% x1=rho(:,1);
% x2=rho(:,2);
% y=D(:,1);
%
% s=scatter(x1,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [255 180 150]/255;
% s.MarkerFaceColor = [255 180 150]/255;
%
% l=lsline;
% l(1).LineWidth= 5;
% l(1).Color=[240 120 100]/255;
%
% hold on
%
% s=scatter(x2,y,'o');
% s.MarkerFaceAlpha   = 0.9;
% s.MarkerEdgeAlpha   = 0.9;
% s.SizeData          = 100;
% s.MarkerEdgeColor = [150 220 220]/255;
% s.MarkerFaceColor = [150 220 220]/255;
%
%
% l=lsline;
% l(1).LineWidth= 5;
% l(1).Color=[120 200 200]/255;
%
%
% set(gca,'fontsize',16,'box','off','lineWidth',2)
% xlabel(' \rho ' )
% ylabel(' dPrime ')
% print(gcf, '-dpdf', '../plots/xdiva/MeanVecCWEHMReldP')
% %%
% % for ss = 1:nSubjs
% %     c=compass(z(ss,1));
% %     c.MarkerSize=10;
% %     c.LineWidth=1;
% %     c.Color=[255 180 150]/255;
%     hold on;
%
%     c=compass(z(ss,2));
%     c.MarkerSize=10;
%     c.LineWidth=1;
%     c.Color=[150 220 220]/255;
% end

%axis tight

%plot(1:3,mean(conf),'linewidth',3,'color',[0.1 0.1 0.1]);
%set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'})
%xlim([0.5 3.5])
%set(gca,'LineWidth',2)
%ylabel('Acc')

%print(gcf,'-dpdf',['../plots/xdiva/ConfidenceAcc'])





