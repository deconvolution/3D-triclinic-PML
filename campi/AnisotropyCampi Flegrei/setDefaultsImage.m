function setDefaultsImage(xLimits,yLimits,velocityLimits,axNumber,mapRB)
set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',12);

caxis(velocityLimits)
xlim(xLimits)
ylim(yLimits)
xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',20)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',20)
colormap(axNumber,mapRB)
shading interp



