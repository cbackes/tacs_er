function han = PolarPlot(th,rho,opts)
% make polar scatter plot for multiple categories

if any(size(th)==1)
    th = th(:);
    rho = rho(:);
    nCategories=1;
else
    nCategories     = size(th,2);
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

% marker sizes
if ~isfield(opts,'markerSize')
    markerSize = ones(1,nCategories)*100;
else
    markerSize = opts.markerSize;
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
    xpos = xLims(2)*0.75; xstr = round(xpos*100)/100;
    t=text(xpos,0,num2str(xstr));
    t.FontSize=18;
    t.HorizontalAlignment='center';
    t.VerticalAlignment='top';
end

% scatter of dots by category
for jj = 1:nCategories
    s=scatter(reZ(:,jj),imZ(:,jj),'o');
    s.MarkerFaceAlpha   = 0.7;
    s.MarkerEdgeAlpha   = 0.7;
    s.SizeData          = markerSize(:,jj);
    s.MarkerEdgeColor   = colors(jj,:);
    s.MarkerFaceColor   = colors(jj,:);    
end

% resulting mean vectors
if meanVecs
    for jj = 1:nCategories
        plot([0 real(zM(jj))],[0 imag(zM(jj))],'linewidth',5,'color',colors(jj,:)*0.8)
    end
end

if connect
    for ii=1:N
        plot(reZ(ii,:),imZ(ii,:),'color',[0.6 0.6 0.6])
    end
end

axis off

return