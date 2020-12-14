function setDefaultColorbars(ax1,ax2,polarizationLimits,mapRB,...
    titleImage,titlePolarization)
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
cb1 = colorbar(ax1,'Position',[.05 .11 .065 .815],'FontSize',16);
cb1.Label.String = titleImage;%'Group Velocity'
cb1.Label.FontSize = 20;
cb1.Label.FontWeight = 'bold';

cb2 = colorbar(ax2,'Position',[.88 .11 .065 .815],'FontSize',16);
cb2.Label.String = titlePolarization;%'R'
cb2.Label.FontSize = 20;
cb2.Label.FontWeight = 'bold';



