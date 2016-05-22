function han = xyScatter(x,y,opts)

nCategories     = size(x,2);

if nargin <3
    opts=[];
end
if  ~isfield(opts,'colors')
    colors = unifrnd(0.4,0.8,[nCategories 3]);
else 
    colors = opts.colors;
end

if ~isfield(opts,'markerSize')
    markerSize = ones(1,nCategories)*100;
else
    markerSize = opts.markerSize;
end

han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[400 400],'paperposition',[0 0 400 400])
set(gcf,'position',[75,75,250,250])
hold on

maxX = max(x(:)); minX = min(x(:));
maxY = max(y(:)); minY = min(y(:));

dX = abs(maxX-minX);
dY = abs(maxY-minY);
xLims = [minX-dX/10 maxX+dX/10];
xlim(xLims)
yLims = [minY-dY/10 maxY+dY/10];
ylim(yLims)
% lCol = [0.9 0.9 0.9];
% plot(xLims,[0 0],'color',lCol,'linewidth',2);
% plot([0 0],yLims,'color',lCol,'linewidth',2);
%

for jj = 1:nCategories
    
    s=scatter(x(:,jj),y(:,jj), 'o');
    s.MarkerFaceAlpha   = 0.7;
    s.MarkerEdgeAlpha   = 0.7;
    s.SizeData          = markerSize(:,jj);
    s.MarkerEdgeColor   = colors(jj,:);
    s.MarkerFaceColor   = colors(jj,:);
    
    if isfield(opts,'polyfitN')
        n = opts.polyfitN;
        P = polyfit(x(:,jj),y(:,jj),n);
        xx = linspace(min(x(:,jj)),max(x(:,jj)),100);
        yy = zeros(size(xx));
        for ii = 1:(n+1)
            yy = yy+P(ii)*xx.^(n-ii+1);
        end
        plot(xx,yy,'linewidth',5,'color',colors(jj,:)*0.8);
    end
end

if isfield(opts,'xlabel')
    xlabel(opts.xlabel);
end
if isfield(opts,'ylabel')
    ylabel(opts.ylabel);
end
if isfield(opts,'xytext')
    xy = opts.xytext;
    text(xy(1),xy(2), opts.text , 'fontsize', 18)
end

set(gca,'LineWidth',2,'FontSize',20)