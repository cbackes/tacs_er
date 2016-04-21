dataPath    = '~/Google Drive/Research/tACS/tACS_ER_task/data/tacs_enc_xdiva/';
load([dataPath 'Summary/BehavSummary.mat']) 
load([dataPath 'Summary/PhaseDependentAnalyses.mat'])

%% dPrimes
nSubjs  = numel(behav_out.retSubj);
Dstrs   = {'dPrime','Face_dPrime','Scene_dPrime',...
    'dPrime_C','Face_dPrime_C','Scene_dPrime_C'};

Dstrs2  = {'dP', 'Face dP', 'Scn dP', 'C', 'Face C', 'Scn C'};
nD = numel(Dstrs);
D       = zeros(nSubjs,nD);
for ii = 1:nD
    D(:,ii) = behav_out.retSummary.(Dstrs{ii});
end

figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,800,400])

xPosCore = [50 160 270];
xPos = [xPosCore xPosCore+400];
a     = zeros(nD,1);
for ii =1:nD
    a(ii)=axes('units','points','position',[xPos(ii) 100 80 250]);
end
yLims = [0 2; 0 2; 0 2; -1.2 1.2; -1.2 1.2; -1.2 1.2];
yTicks1 = [0:0.5:2];
for ii=1:nD
    axes(a(ii))
    set(gca,'Position',[xPos(ii) 100 80 250])
    y = D(:,ii);
    
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
    xlabel(Dstrs2{ii})
    ylim(yLims(ii,:))
    xlim([0 nSubjs+1])
    if ii==1
        set(gca,'ytick',yTicks1)
    end    
    if ismember(ii,[2 3 5 6])
        set(gca,'ycolor','none')
    end
    set(gca,'LineWidth',2)
end
print(gcf,'-dpdf',['../plots/xdiva/dPrimes'])

%% Confidence
conf = behav_out.retSummary.meanAccuracyByConf;

figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,500,300])

hold on;
for ss = 1:nSubjs
    p=plot(1:3,conf(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);    
end
plot(1:3,mean(conf),'linewidth',3,'color',[0.1 0.1 0.1]);    
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'}) 
xlim([0.5 3.5])
set(gca,'LineWidth',2)
ylabel('Acc')

print(gcf,'-dpdf',['../plots/xdiva/ConfidenceAcc'])

%% dPrime by confidence:

DPC = behav_out.retSummary.dPrimeConf;
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,500,300])

hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPC(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);    
end
plot(1:3,mean(DPC),'linewidth',3,'color',[0.1 0.1 0.1]);    
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'}) 
xlim([0.5 3.5])
set(gca,'LineWidth',2)
ylabel(' dPrime ')

print(gcf,'-dpdf',['../plots/xdiva/Confidence_dPrime'])

%% dPrime by confidence:

DPC = behav_out.retSummary.dPrimeConf;
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[50,500,500,300])

hold on;
for ss = 1:nSubjs
    p=plot(1:3,DPC(ss,:),'linewidth',1,'color',[0.8 0.8 0.8]);    
end
plot(1:3,mean(DPC),'linewidth',3,'color',[0.1 0.1 0.1]);    
set(gca,'fontsize',20,'xTick',1:3,'xticklabel',{'Low','Med','High'}) 
xlim([0.5 3.5])
set(gca,'LineWidth',2)
ylabel(' dPrime ')

print(gcf,'-dpdf',['../plots/xdiva/Confidence_dPrime'])

%% RTs
strs    = {'median_Hit_RTs','median_CRs_RTs','median_FA_RTs', 'median_Miss_RTs'};
strs2   = {'Hits','CRs','FA','Misses'};
nRTConds = numel(strs);

RTs       = zeros(nSubjs,nRTConds);
for ii = 1:nRTConds
    RTs(:,ii) = behav_out.retSummary.(strs{ii});
end
figure(1); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 500 300])
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

print(gcf,'-dpdf',['../plots/xdiva/RTs'])

%% PropHits PropMisses and Difference by Phase

t = linspace(0,1,1000);
x = cos(2*pi*t-pi);
xa = angle(hilbert(x));

% Hits
figure(4); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot((xa+pi)./pi*180,x,'k','linewidth',4)
axis tight; 
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase;
plot([36:72:360],X','-','color',[255 180 150]/255)
plot([36:72:360],mean(X), 'color',[240 80 40]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (hits) ' )

print(gcf, '-dpdf', '../plots/xdiva/HitsPropByPhase');

%% Misses
figure(5); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot((xa+pi)./pi*180,x,'k','linewidth',4)
axis tight; 
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propMissByPhase;
plot([36:72:360],X','-','color',[150 220 220]/255)
plot([36:72:360],mean(X), 'color',[120 200 200]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([0.1 0.3])
ylabel(' proportion (misses) ' )

print(gcf, '-dpdf', '../plots/xdiva/MissPropByPhase');

%% Difference
figure(6); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot((xa+pi)./pi*180,x,'k','linewidth',4)
axis tight; 
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.propHitsByPhase-out.propMissByPhase;
plot([36:72:360],X','-','color',[180 180 180]/255)
plot([36:72:360],mean(X), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([-0.15 0.15])
ylabel(' p(h)-p(m) ' )

print(gcf, '-dpdf', '../plots/xdiva/H-MPropByPhase');

%% MemScore by Phase
figure(7); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[1000 600],'paperposition',[0 0 1000 600])
set(gcf,'position',[100,100,600,400])

a2 = axes('units','points','position',[0.12*600 0.1*400 0.8*600 0.15*400]);
axes(a2)
plot((xa+pi)./pi*180,x,'k','linewidth',4)
axis tight; 
set(gca,'ytick',[],'ycolor','w','fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[0:72:360])
xlabel(' Encoding Phase (deg)')
grid on

a1 = axes('position', [0.12 0.3 0.8 0.6]);
axes(a1); hold on;
X = out.MemScoreByPhase;
Xm=(X-repmat(mean(X,2),[1,5]));
plot([36:72:360],Xm','-','color',[180 180 180]/255)
plot([36:72:360],mean(Xm), 'color',[100 100 100]/255,'linewidth',5)
set(gca,'fontsize',16,'box','off','lineWidth',2)
set(gca,'xtick',[36:72:360],'xTickLabel','')
xlim([0 360])
ylim([-1 1])
ylabel(' MemScore ' )

print(gcf, '-dpdf', '../plots/xdiva/MemScoreByPhase');

%% Hit Miss Phase distributions by subject

th = mod(out.HitMissMePhase,2*pi);
rho = out.HitMissAbPhase;

z = rho.*exp(1i*th);
zdot = rho(:,1).*rho(:,2).*cos(th(:,1)-th(:,2));
f1=@(x)(angle(mean(exp(1j*x))));
f2=@(x)angle(mean(x(:,1)))-angle(mean(x(:,2)))

CIdmth = bootci(1000,f2,z);
CId=bootci(1000,f1,th(:,1)-th(:,2));
CIh=bootci(1000,f1,th(:,1));
CIm=bootci(1000,f1,th(:,2));


figure(8); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

xLims = [-0.12 0.12];
xlim(xLims);
ylim(xLims); hold on;
lCol = [0.8 0.8 0.8];
plot(xlim,[0 0],'color',[0.8 0.8 0.8],'linewidth',2);
plot([0 0],ylim,'color',lCol,'linewidth',2);
plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);

reZ = real(z);
imZ = imag(z);

s=scatter(reZ(:,1),imZ(:,1),'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [255 180 150]/255;
s.MarkerFaceColor = [255 180 150]/255;

s=scatter(reZ(:,2),imZ(:,2),'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [150 220 220]/255;
s.MarkerFaceColor = [150 220 220]/255;

plot([0 mean(reZ(:,1))],[0 mean(imZ(:,1))],'linewidth',5,'color',[240 120 100]/255)
plot([0 mean(reZ(:,2))],[0 mean(imZ(:,2))],'linewidth',5,'color',[120 200 200]/255)

axis off

print(gcf, '-dpdf', '../plots/xdiva/MeanVecBySubj1');

%% line between points
figure(8); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

xLims = [-0.12 0.12];
xlim(xLims);
ylim(xLims); hold on;
lCol = [0.8 0.8 0.8];
plot(xlim,[0 0],'color',[0.8 0.8 0.8],'linewidth',2);
plot([0 0],ylim,'color',lCol,'linewidth',2);
plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);

reZ = real(z);
imZ = imag(z);

for ss=1:nSubjs
    plot(reZ(ss,:),imZ(ss,:),'color',[0.2 0.2 0.2])
end

s=scatter(reZ(:,1),imZ(:,1),'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [255 180 150]/255;
s.MarkerFaceColor = [255 180 150]/255;

s=scatter(reZ(:,2),imZ(:,2),'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [150 220 220]/255;
s.MarkerFaceColor = [150 220 220]/255;

axis off
print(gcf, '-dpdf', '../plots/xdiva/MeanVecBySubj2');
%% difference in hit phase and miss phase
th = mod(out.HitMissMePhase,2*pi);
rho = out.HitMissAbPhase;

z = rho.*exp(1i*th);
zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));

figure(10); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

xLims = [-0.1 0.1];
xlim(xLims);
ylim(xLims); hold on;
lCol = [0.8 0.8 0.8];
plot(xlim,[0 0],'color',[0.8 0.8 0.8],'linewidth',2);
plot([0 0],ylim,'color',lCol,'linewidth',2);
plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);

reZ = real(zD);
imZ = imag(zD);

s=scatter(reZ,imZ,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [180 180 180]/255;
s.MarkerFaceColor = [180 180 180]/255;

plot([0 mean(reZ)],[0 mean(imZ)],'linewidth',5,'color',[20 20 20]/255);

axis off
print(gcf, '-dpdf', '../plots/xdiva/MeanVecBySubj3');

%% relationship of bias to magnitude of mean vector:
th = mod(out.HitMissMePhase,2*pi);
rho = out.HitMissAbPhase;

z = rho.*exp(1i*th);
%zD = z(:,1)./z(:,2);
zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));

figure(11); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

x=rho(:,1)-rho(:,2);
y=D(:,4);

s=scatter(x,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [180 180 180]/255;
s.MarkerFaceColor = [180 180 180]/255;

set(gca,'fontsize',16,'box','off','lineWidth',2)
xlabel(' \rho_{Hits}-\rho_{Misses} ' )
ylabel(' Bias ')
print(gcf, '-dpdf', '../plots/xdiva/MeanVecRelBias');
%%  relationship of d' to magnitude of mean vector:
th = mod(out.HitMissMePhase,2*pi);
rho = out.HitMissAbPhase;

z = rho.*exp(1i*th);
%zD = z(:,1)./z(:,2);
zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));

figure(11); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

x=abs(rho(:,1)-rho(:,2));
y=D(:,1);

s=scatter(x,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [180 180 180]/255;
s.MarkerFaceColor = [180 180 180]/255;

set(gca,'fontsize',16,'box','off','lineWidth',2)
xlabel(' |\rho_{Hits}-\rho_{Misses}| ' )
ylabel(' dPrime ')
print(gcf, '-dpdf', '../plots/xdiva/MeanVecReldP')

%%  relationship of d' to magnitude of mean vector for misses and hits
th = mod(out.HitMissMePhase,2*pi);
rho = out.HitMissAbPhase;

z = rho.*exp(1i*th);
%zD = z(:,1)./z(:,2);
zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));

figure(11); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

x1=rho(:,1);
x2=rho(:,2);
y=D(:,1);

s=scatter(x1,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [255 180 150]/255;
s.MarkerFaceColor = [255 180 150]/255;

l=lsline;
l(1).LineWidth= 5;
l(1).Color=[240 120 100]/255;

hold on

s=scatter(x2,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [150 220 220]/255;
s.MarkerFaceColor = [150 220 220]/255;


l=lsline;
l(1).LineWidth= 5;
l(1).Color=[120 200 200]/255;


set(gca,'fontsize',16,'box','off','lineWidth',2)
xlabel(' \rho ' )
ylabel(' dPrime ')
print(gcf, '-dpdf', '../plots/xdiva/MeanVecHMReldP')

%% MEeanVecDiff WEighted by confidence.
th = mod(out.HitMissCWMePhase,2*pi);
rho = out.HitMissCWAbPhase;

z = rho.*exp(1i*th);
%zD = z(:,1)./z(:,2);
zD = abs((rho(:,1)-rho(:,2))).*exp(1j*(th(:,1)-th(:,2)));

figure(12); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

xLims = [-0.15 0.15];
xlim(xLims);
ylim(xLims); hold on;
lCol = [0.8 0.8 0.8];
plot(xlim,[0 0],'color',[0.8 0.8 0.8],'linewidth',2);
plot([0 0],ylim,'color',lCol,'linewidth',2);
plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);

reZ = real(zD);
imZ = imag(zD);

s=scatter(reZ,imZ,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [180 180 180]/255;
s.MarkerFaceColor = [180 180 180]/255;

plot([0 mean(reZ)],[0 mean(imZ)],'linewidth',5,'color',[20 20 20]/255);

axis off
print(gcf, '-dpdf', '../plots/xdiva/MeanVecCweBySubj3');

%% Mean vector weighted by confidence relation to bias
th = mod(out.HitMissCWMePhase,2*pi);
rho = out.HitMissCWAbPhase;

z = rho.*exp(1i*th);
%zD = z(:,1)./z(:,2);

figure(13); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

x=rho(:,1)-rho(:,2);
y=D(:,4);

s=scatter(x,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [180 180 180]/255;
s.MarkerFaceColor = [180 180 180]/255;

set(gca,'fontsize',16,'box','off','lineWidth',2)
xlabel(' log_{10}(|hit/miss|) ' )
ylabel(' Bias ')
print(gcf, '-dpdf', '../plots/xdiva/MeanVecCWeRelBias');

%% mean vectors (h/miss) weighted by conf relation to dP
th = mod(out.HitMissCWMePhase,2*pi);
rho = out.HitMissCWAbPhase;

figure(11); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[100,100,500,500])

x1=rho(:,1);
x2=rho(:,2);
y=D(:,1);

s=scatter(x1,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [255 180 150]/255;
s.MarkerFaceColor = [255 180 150]/255;

l=lsline;
l(1).LineWidth= 5;
l(1).Color=[240 120 100]/255;

hold on

s=scatter(x2,y,'o');
s.MarkerFaceAlpha   = 0.9;
s.MarkerEdgeAlpha   = 0.9;
s.SizeData          = 100;
s.MarkerEdgeColor = [150 220 220]/255;
s.MarkerFaceColor = [150 220 220]/255;


l=lsline;
l(1).LineWidth= 5;
l(1).Color=[120 200 200]/255;


set(gca,'fontsize',16,'box','off','lineWidth',2)
xlabel(' \rho ' )
ylabel(' dPrime ')
print(gcf, '-dpdf', '../plots/xdiva/MeanVecCWEHMReldP')
%%
% for ss = 1:nSubjs
%     c=compass(z(ss,1)); 
%     c.MarkerSize=10;
%     c.LineWidth=1;
%     c.Color=[255 180 150]/255;
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





