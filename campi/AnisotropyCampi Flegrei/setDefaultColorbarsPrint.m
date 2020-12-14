function setDefaultColorbarsPrint(ax1,ax2,polarizationLimits,mapRB,...
    titleImage,titlePolarization,sz)
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

cb2 = colorbar(ax2,'Position',[.88 .11 .065 .815],'FontSize',16);
cb2.Label.String = titlePolarization;%'R'
cb2.Label.FontSize = sz;
cb2.Label.FontWeight = 'bold';
cb2.FontSize = sz-2;
set(cb2, 'YAxisLocation','left')

ax1.Position = ax1.Position - [0.05 0 0 0];
ax2.Position = ax2.Position - [0.05 0 0 0];

