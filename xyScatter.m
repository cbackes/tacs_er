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
    markersize = opts.markersize;
end

han=figure(); clf;
set(gcf,'paperpositionmode','auto','color','white')
set(gcf,'paperUnits','points','papersize',[800 800],'paperposition',[0 0 800 800])
set(gcf,'position',[110,100,500,500])
hold on

dX = abs(max(x)-min(x));
xLims = [min(x)-dX/10 max(x)+dX/10];
xlim(xLims)
yLims = [min(y)-dX/10 max(y)+dX/10];
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
        P = polyfit(x,y,n);
        xx = linspace(min(x),max(x),100);
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