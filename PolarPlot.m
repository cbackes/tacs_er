function han = PolarPlot(th,rho,opts)
% make polar scatter plot for multiple categories

if nargin>=2
    if isempty(rho)
        rho = ones(size(th));
    end
    if any(size(th)==1)
        th = th(:);
        rho = rho(:);
        nCategories=1;
    else
        nCategories     = size(th,2);
    end
end

if nargin<2
    if any(size(th)==1)
        th = th(:);        
        nCategories=1;
    else
        nCategories     = size(th,2);
    end
    opts = [];
    opts.maxR = 4/3;
    rho = ones(size(th));    
end

if nargin==2
    opts=[];
end

N       = size(th,1);
z       = rho.*exp(1i*th);
reZ     = real(z);
imZ     = imag(z);
zM      = nanmean(z);

% colors by category
if ~isfield(opts,'colors')
    colors = unifrnd(0.5,1,[nCategories 3]);
else
    colors = opts.colors;
end

if ~isfield(opts,'alpha')
    alpha = 0.7;
else     
    alpha = opts.alpha;
end

% marker sizes
if ~isfield(opts,'markerSize')
    markerSize = ones(N,nCategories)*200;
elseif size(opts.markerSize,1)>1
    markerSize = opts.markerSize;
else
    markerSize = opts.markerSize(1)*ones(N,nCategories);
end

% max Radius
if ~isfield(opts,'maxR')
    xLims = [-1 1]*max(rho(:));
else
    xLims =  [-1 1]*opts.maxR;
end

% polar grid
if ~isfield(opts,'polarGrid')
    polarGrid = 1;
else
    polarGrid = opts.polarGrid;
end

if ~isfield(opts,'magText')
    magText = 1;
else
    magText = opts.magText;
end

% mean Vector
if ~isfield(opts,'meanVecs')
    meanVecs = 1;
else
    meanVecs = opts.meanVecs;
end

%  connect Categories
if ~isfield(opts,'connect')
    connect = 0;
else
    connect = opts.connect;
end
% set figure
han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[600 600],'paperposition',[0 0 600 600])
set(gcf,'position',[100,100,400,400])
xlim(xLims);
ylim(xLims); hold on;

% polar grid
if polarGrid
    lCol = [0.9 0.9 0.9];
    plot(xlim,[0 0],'color',lCol,'linewidth',2);
    plot([0 0],ylim,'color',lCol,'linewidth',2);
    plot(xLims*0.67,xLims*0.67,'color',lCol,'linewidth',2);
    plot(xLims*0.67,[xLims(2) xLims(1)]*0.67,'color',lCol,'linewidth',2);
    
    % reference circle for magnitude
    t = 0:0.01:2*pi;
    xx = xLims(1)*cos(t)*0.75;
    yy = xLims(1)*sin(t)*0.75;
    plot(xx,yy,'-','color',lCol,'linewidth',1);
    xpos = xLims(2)*0.75;
end

% scatter of dots by category
K = 10;
idx=crossvalind('kfold',N,10);
for kk = 1:K
    % re-order categories
    catVec = randperm(nCategories);
    for jj = catVec
        s=scatter(reZ(idx==kk,jj),imZ(idx==kk,jj),'o');
        s.MarkerFaceAlpha   = alpha;
        s.MarkerEdgeAlpha   = 0;
        s.SizeData          = markerSize(idx==kk,jj);
        s.MarkerEdgeColor   = colors(jj,:);
        s.MarkerFaceColor   = colors(jj,:);
    end
end
% for jj = 1:nCategories
%     s=scatter(reZ(:,jj),imZ(:,jj),'o');
%     s.MarkerFaceAlpha   = alpha;
%     s.MarkerEdgeAlpha   = 0;
%     s.SizeData          = markerSize(:,jj);
%     s.MarkerEdgeColor   = colors(jj,:);
%     s.MarkerFaceColor   = colors(jj,:);
% end

% resulting mean vectors
if meanVecs
    for jj = 1:nCategories
        plot([0 real(zM(jj))],[0 imag(zM(jj))],'linewidth',5,'color',colors(jj,:)*0.9)
    end
end
if magText
    xpos = xLims(2)*0.75;
    xstr = round(xpos*100)/100;
    t=text(xpos,0,num2str(xstr));
    t.FontSize=18;
    t.HorizontalAlignment='center';
    t.VerticalAlignment='top';
end

if connect
    for ii=1:N
        plot(reZ(ii,:),imZ(ii,:),'color',[0.6 0.6 0.6])
    end
end

axis off

return