function setDefaultsImagePrint(xLimits,yLimits,velocityLimits,axNumber,mapRB,sz)
set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',sz);

caxis(velocityLimits)
xlim(xLimits)
ylim(yLimits)
xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',sz)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',sz)
colormap(axNumber,mapRB)
shading interp



