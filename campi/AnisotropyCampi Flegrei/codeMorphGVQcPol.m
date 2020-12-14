%Script to image interferometric and polarization parameters together
clear; clc; close all
warning off all
load S11_13.mat
S83_84                          =   load('Seismicity_UTM_1983_1984.rtf');
% S05_16=load('Seismicity_UTM_2005_2016.txt');
load polarization.mat
Q                               =   load('QcH.txt');
%% 2 s vs 0.2-1 Hz

xLimits                         =   [420000 432000];
yLimits                         =   [4518000 4524000];
velocityLimits                  =   [0 1.5];
polarizationLimits              =   [0.2 0.6];

% Velocity measurements
X                               =   load('2P.txt');

R                               =   6; %R polarization parameter
tresholdR                       =   0.25;
Az                              =   4; %Az polarization parameter

conditionalX                    =   (X(:,1)<xLimits(1) |...
    X(:,1)>xLimits(2) |...
    X(:,2)<yLimits(1) | X(:,2)>yLimits(2));
X(conditionalX,:)=[];

conditionalP1                   =   (P011(:,1)<xLimits(1) |...
    P011(:,1)>xLimits(2) |...
    P011(:,2)<yLimits(1) | P011(:,2)>yLimits(2));

P011(conditionalP1,:)           =   [];
cVect                           =   X(:,4);

PVect                           =   P011(:,R);
PAz                             =   P011(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

interp                          =   100;

sz                              =   24;

[xi,yi]                         =...
    meshgrid(xLimits(1):interp:xLimits(2), yLimits(1):interp:yLimits(2));
zi                              =   griddata(X(:,1),X(:,2),cVect,xi,yi);

C                               =   shaperead('COSs.shp');
F                               =   shaperead('FAULTS.shp');

figure('Name',...
    'Group Velocity vs Residual Length and Azimuth, 2 s vs 0.2-1 Hz',...
    'NumberTitle','off','Position',[10 10 1300 560]);
ax1                             =   axes;
pcolor(xi,yi,zi);
chosenCMap                      =   flipud(inferno);

setDefaultsImagePrint(xLimits,yLimits,velocityLimits,ax1,chosenCMap,sz)

hold on

plot(S11_13(:,2),S11_13(:,3),'d','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[1 1 1],'MarkerSize',8,'LineWidth',1)

mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

ax2 = axes;

scatter(ax2,P011(:,1), P011(:,2), 300, PVect, 'Filled',...
    'MarkerEdgeColor','k','LineWidth',2);

hold on
conditionP = P011(:,R)>tresholdR;

quiver(P011(conditionP,1), P011(conditionP,2),...
    cX(conditionP),sX(conditionP),.35,'-k','LineWidth',2)
setDefaultColorbarsPrintDifferent(ax1,ax2,polarizationLimits,chosenCMap,...
    'Group Velocity',sz)
hold off

print('GV_vs_P_2s_02_1Hz','-dtiff','-r300');

%% 0.2-1 s vs 1-5 Hz
conditionalP2                   =   (P15(:,1)<xLimits(1) |...
    P15(:,1)>xLimits(2)| P15(:,2)<yLimits(1) | P15(:,2)>yLimits(2));

P15(conditionalP2,:)            =   [];
PVect                           =   P15(:,R);
PAz                             =   P15(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));

figure('Name',...
    'Group Velocity vs Residual Length and Azimuth, 0.9 s vs 0.2-1 Hz',...
    'NumberTitle','off','Position',[10 10 1300 560]);
ax1                             =   axes;

pcolor(xi,yi,zi);
setDefaultsImagePrint(xLimits,yLimits,velocityLimits,ax1,chosenCMap,sz)

hold on

plot(S11_13(:,2),S11_13(:,3),'d','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[1 1 1],'MarkerSize',8,'LineWidth',1)

mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

ax2                             =   axes;

scatter(ax2,P15(:,1), P15(:,2), 300, PVect, 'Filled',...
    'MarkerEdgeColor','k','LineWidth',2);


hold on

conditionP = P15(:,R)>tresholdR;

quiver(P15(conditionP,1), P15(conditionP,2),...
    cX(conditionP),sX(conditionP),.35,'-k','LineWidth',2)

setDefaultColorbarsPrintDifferent(ax1,ax2,polarizationLimits,chosenCMap,...
    'Group Velocity',sz)

hold off
print('GV_vs_P_09s_1_5Hz','-dtiff','-r300');

%% Q 2-4 Hz vs 1-5 Hz
xLimitsQ                        =   [xLimits(1)+100 xLimits(2)-100];
yLimitsQ                        =   [yLimits(1)+100 yLimits(2)-100];
QLimits                         =   [0 0.008];

[xi,yi]                         =...
    meshgrid(xLimitsQ(1):interp:xLimitsQ(2),...
    yLimitsQ(1):interp:yLimitsQ(2));


conditionalQ                    =   (Q(:,1)<xLimits(1) |...
    Q(:,1)>xLimits(2) | Q(:,2)<yLimits(1) | Q(:,2)>yLimits(2));
Q(conditionalQ,:)               =   [];
cVectQ                          =   Q(:,4);

S                               =   S83_84(804:986,:);
depth                           =	S(:,4)>-5000 & S(:,4)<-2200;

chosenCMapQ                     =   inferno;


[xiQ,yiQ]                       =...
    meshgrid(xLimitsQ(1):interp:xLimitsQ(2),...
    yLimitsQ(1):interp:yLimitsQ(2));
ziQ                             =...
    griddata(Q(:,1),Q(:,2),cVectQ,xi,yi,'natural');
zi                              =   interp2(xiQ,yiQ,ziQ,xi,yi);
zi                              =   inpaintn(zi);

figure('Name',...
    'Qc vs Residual Length and Azimuth, 2-4 Hz vs 1-5 Hz',...
    'NumberTitle','off','Position',[10 10 1300 560]);

ax1                             =   axes;
pcolor(xi,yi,zi);
hold on

plot(S(depth,2),S(depth,3),'d','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[1 1 1],'MarkerSize',8,'LineWidth',1)

setDefaultsImagePrint(xLimits,yLimits,QLimits,ax1,chosenCMapQ,sz)

hold on

mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

ax2 = axes;

scatter(ax2,P15(:,1), P15(:,2), 300, PVect, 'Filled',...
    'MarkerEdgeColor','w','LineWidth',2);

hold on
conditionP = P15(:,R)>tresholdR;

quiver(P15(conditionP,1), P15(conditionP,2),...
    cX(conditionP),sX(conditionP),.35,'-w','LineWidth',2)


setDefaultColorbarsPrintDifferent(ax1,ax2,polarizationLimits,chosenCMap,...
    'Coda Attenuation',sz)

hold off

print('A_vs_P_3Hz_1_5Hz','-dtiff','-r300');

%% Geomorphology
load polarization.mat
xLimits                         =   [420000 432000];
yLimits                         =   [4516000 4528000];

figure('Name','Geomorphology','NumberTitle','off',...
    'Position',[10 10 1300 560]);

subplot(1,2,1)
mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

hold on
P=P011;
PAz                             =   P(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));
conditionP = P(:,R)>tresholdR;
scatter(P(:,1), P(:,2), 200, P(:,R), 'Filled',...
    'MarkerEdgeColor','k','LineWidth',2);
colormap(chosenCMap)
caxis(polarizationLimits)
hold on
quiver(P(conditionP,1), P(conditionP,2),...
    cX(conditionP),sX(conditionP),.6,'-k','LineWidth',2)

ax=gca;

set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',14);

ax.Position = ax.Position - [0.07 0 0 0];

subplot(1,2,2)
mapshow(C,'FaceAlpha',0);
mapshow(F,'Color', 'black')

hold on
P=P15;
PAz                             =   P(:,Az)-90;
cX                              =   cos(deg2rad(PAz));
sX                              =   -sin(deg2rad(PAz));
conditionP = P(:,R)>tresholdR;
scatter(P(:,1), P(:,2), 200, P(:,R), 'Filled',...
    'MarkerEdgeColor','k','LineWidth',2);
colormap(chosenCMap)
caxis(polarizationLimits)
hold on
quiver(P(conditionP,1), P(conditionP,2),...
    cX(conditionP),sX(conditionP),.6,'-k','LineWidth',2)
set(gca,'XTick',xLimits(1):4000:xLimits(2));
set(gca,'YTick',yLimits(1):2000:yLimits(2));
set(gca,'FontSize',14);

ax=gca;
ax.Position = ax.Position - [0.07 0 0 0];


cb1=colorbar('Position',[.9 .11 .065 .815],'FontSize',16);
cb1.Label.String = 'R';%'Group Velocity'
cb1.Label.FontSize = sz;
cb1.Label.FontWeight = 'bold';
cb1.FontSize = sz-2;

x1=get(gca,'Position');
x=get(cb1,'Position');
x(3)=0.02;
set(cb1,'Position',x)
set(gca,'position',x1)

hold off

print('Morphology','-dtiff','-r300');
