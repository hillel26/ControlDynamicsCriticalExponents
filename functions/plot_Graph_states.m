function plot_Graph_states(G,h,color,zlims,forced)
% G is graph struct 
% G.Nodes.Activity is column vector of a node activity

figure(h);
hh = plot(G,'NodeColor',color,'NodeLabel',{},'MarkerSize',5,'EdgeAlpha',0.1,'LineWidth',1);
z = G.Nodes.Activity;
hh.ZData = z;
x = hh.XData; y = hh.YData;
if ~isempty(forced); highlight(hh,forced,'NodeColor',[0 50 100]/255); end
xlim([min(x),max(x)]); ylim([min(y),max(y)]); zlim(zlims)
set(gca,'XTick',[],'YTick',[],'ZTick',[],'Visible','on','FontSize',28,'linewidth',2,'box','on');
colormap(flipud(jet))
% hc.FontSize = 20;
% hc.Title.String = 'Activity';
% caxis(clim)
set(gcf,'color','w')
view(55,15)
