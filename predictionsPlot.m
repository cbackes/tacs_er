function predictionPlot()

%% PDFs
rng(4);
th = 0:0.01:2*pi-0.01;
N = numel(th);
dx = pi/2;

th1 = th-dx;
th2 = th-pi-dx;
kp1  = 1;
kp2  = 1;
rho1 = exp(kp1*cos(th1))/(2*pi*besseli(0,kp1));
rho2 = exp(kp2*cos(th2))/(2*pi*besseli(0,kp2));

opts = [];
opts.markerSize = [50 50];
opts.meanVecs   = 0;
opts.colors     = [255 180 150; 150 220 220]/255;
% raw vectors
han = PolarPlot([th' th'],[rho1' rho2'],opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])
set(gca,'units','points','position',[50 50 300 300]); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ProbabilityPlot_dx' ....
    num2str(round(dx/pi*180)) 'kp1_' num2str(kp1) 'kp2_' num2str(kp2) ]);

%% samples from PDFs
th = 0:0.05:2*pi-0.01;
N = numel(th);
dx = -pi/4;

th1 = th-dx;
th2 = th-pi-dx;
kp1  = 2;
kp2  = 0.5;
rho1 = exp(kp1*cos(th1))/(2*pi*besseli(0,kp1));
rho2 = exp(kp2*cos(th2))/(2*pi*besseli(0,kp2));

nn = 1000;
th1_r = randsample(th,nn,true,rho1);
th1_r= th1_r(:);
th2_r = randsample(th,nn,true,rho2);
th2_r = th2_r(:);

opts = [];
opts.markerSize = [50 50];
opts.meanVecs   = 1;
opts.colors     = [255 180 150; 150 220 220]/255;
opts.maxR       = 4/3;
opts.alpha      = 0.5;
% raw vectors
nn = numel(th1_r);
han = PolarPlot([th1_r th2_r],ones(nn,2)+randn(nn,2)*0.05,opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])
set(gca,'units','points','position',[50 50 300 300]); 
print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ProbabilityPlot2_dx' ....
    num2str(round(dx/pi*180)) 'kp1_' num2str(kp1) 'kp2_' num2str(kp2) ]);
 
%% Samples from fixed points
th = 0:2*pi/5:2*pi-0.1;
N = numel(th);
dx = -pi/4;

th1 = th-dx;
th2 = th-pi-dx;
kp1  = 0.1;
kp2  = 0.1;
rho1 = exp(kp1*cos(th1))/(2*pi*besseli(0,kp1));
rho2 = exp(kp2*cos(th2))/(2*pi*besseli(0,kp2));

nn = 150;
th1_r = randsample(th,nn,true,rho1);
th1_r= th1_r(:);
th2_r = randsample(th,nn,true,rho2);
th2_r = th2_r(:);

%
%
opts = [];
opts.markerSize = [50 50];
opts.meanVecs   = 1;
opts.colors     = [255 180 150; 150 220 220]/255;
opts.maxR       = 4/3;
opts.alpha      = 0.8;
% raw vectors
nn = numel(th1_r);
han = PolarPlot([th1_r th2_r]+unifrnd(-18/180*pi,18/180*pi,nn,2),ones(nn,2)+randn(nn,2)*0.05,opts);
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])
set(gca,'units','points','position',[50 50 300 300]); 
% print(han, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/ProbabilityPlot2_dx' ....
%     num2str(round(dx/pi*180)) 'kp1_' num2str(kp1) 'kp2_' num2str(kp2) ]);
 

%%
figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 400],'paperposition',[0 0 600 400])
set(gcf,'position',[100,100,400,200])

set(gca,'units','points','position',[50 80 300 80]); hold on;

N  = 500;
th = linspace(0,4*pi,N);
x = cos(th);
c = [255 180 150; 150 220 220]/255;
cols1 = [linspace(c(1,1),c(2,1),N/4)',linspace(c(1,2),c(2,2),N/4)',linspace(c(1,3),c(2,3),N/4)'];
cols  = [cols1;flipud(cols1);cols1;flipud(cols1)];
dth =diff(th);
cols  = circshift(cols,round( dx/dth(1)));
for jj = 1:N
    s = scatter(th(jj),x(jj),'o');   
    s.MarkerFaceAlpha   = 0.5;
    s.MarkerEdgeAlpha   = 0.5;
    s.SizeData          = 100;
    s.MarkerEdgeColor   = cols(jj,:);
    s.MarkerFaceColor   = cols(jj,:);
end
set(gca,'xtick',0:pi:4*pi,'xticklabel',(0:pi:4*pi)/pi*180)
set(gca,'YColor','w','fontsize',20,'linewidth',2)
axis tight
ylim([-1.1 1.1])
xlabel( ' \phi ')
print(gcf, '-dpdf', ['../plots/tacs_enc_xdiva/Predictions/OscillationPlot_dx' num2str(round(dx/pi*180))]);
%axis off



%% Demonstration of HR for one subject



