function setDefaultColorbarsPrintDifferent(ax1,ax2,polarizationLimits,...
    mapRB,titleImage,sz)
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

% Give each one its own colormap
caxis(polarizationLimits)
colormap(ax2,mapRB)
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.88 .11 .065 .815],'FontSize',16);
cb1.Label.String = titleImage;%'Group Velocity'
cb1.Label.FontSize = sz;
cb1.Label.FontWeight = 'bold';
cb1.FontSize = sz-2;
% set(cb1,'XTickLabel',{'8','4','0'}, 'XTick', 0:4*10^-3:8*10^-3)
% cb1.Alpha=0;

% cb2 = colorbar(ax2,'Position',[.88 .11 .065 .815],'FontSize',16);
% cb2.Label.String = titlePolarization;%'R'
% cb2.Label.FontSize = sz;
% cb2.Label.FontWeight = 'bold';
% cb2.FontSize = sz-2;
% set(cb2, 'YAxisLocation','left')
x1=ax1.Position;
x=get(cb1,'Position');
x(3)=0.02;
set(cb1,'Position',x)
set(ax1,'position',x1)


