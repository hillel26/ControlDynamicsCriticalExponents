function plot_Graph_Shells(G,s,h,clim)
% G is graph struct 
% s - source of bfs
% G.Nodes.Activity is column vector of a node activity
if isempty(s); s = randi(G.numnodes); elseif length(s)>1; s=s(1);end 

d = distances(G,s);
% if the graph is not connected
% G = rmnode(G,find(isinf(d))); 
% d(isinf(d))=[];
[d,idx] = sort(d);
G = reordernodes(G,idx);

r = d;
th = 1.5*(0:length(d)-1);
x = r.*cos(th);
y = r.*sin(th);
z = G.Nodes.Activity;

figure(h);
hh = plot(G,'XData',x,'YData',y,'ZData',z,'NodeCData',z,'NodeLabel',{},'MarkerSize',5,'EdgeAlpha',0.3,'LineWidth',1.5);
% highlight(hh,1,'MarkerSize',13,'Marker');
% hold on; scatter3(x(1),y(1),z(1),100,'k','filled'); hold off;
xlim([min(x),max(x)]); ylim([min(y),max(y)]); zlim(clim)
set(gca,'XTick',[],'YTick',[],'Visible','on','FontSize',28,'linewidth',2);
colormap(flipud(jet))
% hc = colorbar;
% hc.FontSize = 20;
% hc.Title.String = 'Activity';
caxis(clim)
set(gcf,'color','w')
view(55,15)
