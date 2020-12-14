%% Anisotropic, viscoelastic SH-wave propagation at Campi Flegrei
% 20x20 km, 1 Hz
%
% O(2,4) finite-difference scheme
% dx=3*10^3 (m/s)/ 10 (Hz) -> 20 m
% ddz=20000 m / 20 m = 1000 nodes
% ddx=20000 m / 20 m = 1000 nodes
% dt=(2/pi)* 20 (m)/ (3*10^3) (m/s)->0.004 s
clear
close all
clc
%% INPUTS
%isotropic, effective anisotropic or scattering
chooseModel                         =   1; 

%velocity fluctuations
rmsVelocityFluctuations             =   0.07;
correlationLength                   =   2;

% V0=sqrt(4.5*10^9/2500)/1000; %Heap if setting from rock physics
% volcano

if chooseModel ==  1
    V0                              =   3;
else
    V0_SN                           =   2;
    V0_WE                           =   4;
    csi_SN                              =   V0_SN*(1+csi);
    csi_WE                              =   V0_WE*(1+csi);

end

spatialSampling                     =   1/100;
nodesX                              =   1000;
nodesZ                              =   1000;

% Material prop. - Homogeneous (1), effective medium (2), scattering (3)

rhoAverage                          =   2500;
rhoScattering                       =   2500;
rho                                 =   rhoAverage*ones(nodesX,nodesZ);
    
Q2Average                           =   100;
Q4Average                           =   100;
Q2Scattering                        =   100;
Q4Scattering                        =   100;

%% Definition of the velocity fluctuations
spatialWavenumber                   =   zeros(nodesZ);
spatialWavenumber(1,:)              =...
    spatialSampling:spatialSampling:nodesZ/100;
spatialWavenumber(:,1)              =...
    spatialSampling:spatialSampling:nodesZ/100;

%dummy
m1                                  =   zeros(nodesZ,1);
for i = 2:nodesZ
    m1(i)=spatialSampling*i;
    spatialWavenumber(2:end,i)      =...
        sqrt(m1(i).^2+spatialWavenumber(2:end,1).^2);
end

%Power spectral density function - exponential
Pexp                                =...
    4*pi*rmsVelocityFluctuations^2*correlationLength^2./...
    (1+correlationLength^2*spatialWavenumber.^2).^2;
phase                               =...
    2*pi.*rand(size(spatialWavenumber));
ePhase                              =   (exp(1i*phase));

%velocity fluctuations
csi                                 =...
    1/4/pi^2*abs(fft2(sqrt(Pexp).*ePhase));
csi                                 =   csi(1:nodesX,1:nodesZ);
csi                                 =   csi-mean(mean(csi));

% for both isotropic and anisotropic cases
csi_h                               =   V0*(1+csi);

%% Material properties - stiffness coefficients c44, c66,c46
if chooseModel == 1
    
    c44                       =   rhoAverage.*(csi_h*1000).^2;
    c66                       =   rhoAverage.*(csi_h*1000).^2;
    c46                       =   zeros(nodesX,nodesZ);
     
elseif chooseModel==2
    
    c44Average                      =   rhoAverage.*(csi_WE*1000).^2;
    c66Average                      =   rhoAverage.*(csi_SN*1000).^2;
    c46Average                      =   -1000000000*ones(nodesX,nodesZ);
%     c46=-0.5*sqrt(c44.*c66);

elseif chooseModel==3

    c44Average                      =   rhoAverage.*(csi_WE*1000).^2;
    c44Scattering                   =   rhoScattering.*(csi_WE*1000).^2;
    c66Average                      =   rhoAverage.*(csi_SN*1000).^2;
    c66Scattering                   =   rhoScattering.*(csi_SN*1000).^2;
    c46Average                      =   -0.5*sqrt(c44Scattering*c66Scattering);
    c46Scattering                   =   0.5*sqrt(c44Scattering*c66Scattering);

end

fRelaxationTimes                    =   10;% To set where Q acts

tau=1/(2*pi*fRelaxationTimes);

ts2Upper=(tau/Q2Average)*(sqrt(Q2Average*Q2Average+1)-1);
te2Upper=(tau/Q2Average)*(sqrt(Q2Average*Q2Average+1)+1);
phi2Upper=1/te2Upper-1/ts2Upper;
ts4Upper=(tau/Q4Average)*(sqrt(Q4Average*Q4Average+1)-1);
te4Upper=(tau/Q4Average)*(sqrt(Q4Average*Q4Average+1)+1);
phi4Upper=1/te4Upper-1/ts4Upper;


ts2=ts2Upper*ones(nodesX,nodesZ);
phi2=phi2Upper*ones(nodesX,nodesZ);
ts4=ts4Upper*ones(nodesX,nodesZ);
phi4=phi4Upper*ones(nodesX,nodesZ);

if chooseModel == 3
    c44=c44Upper*ones(nodesX,nodesZ);
    c44(interx:nxt3,interz:nzt3)=c44Scattering*ones(nxt2,nzt2);
    c66=c66Upper*ones(nxt,nzt);
    c66(interx:nxt3,interz:nzt3)=c66Scattering*ones(nxt2,nzt2);
    c46=c46Upper*ones(nxt,nzt);
    c46(interx:nxt3,interz:nzt3)=c46Lower*ones(nxt2,nzt2);
    rho=rhoUpper*ones(nxt,nzt);
    rho(interx:nxt3,interz:nzt3)=rhoLower*ones(nxt2,nzt2);
    ts2(interx:nxt3,interz:nzt3)=ts2Lower*ones(nxt2,nzt2);
    phi2(interx:nxt3,interz:nzt3)=phi2Lower*ones(nxt2,nzt2);
    ts4(interx:nxt3,interz:nzt3)=ts4Lower*ones(nxt2,nzt2);
    phi4(interx:nxt3,interz:nzt3)=phi4Lower*ones(nxt2,nzt2);
    %     ts2(R1)=ts2Lower;
    %     phi2(R1)=phi2Lower;
    %     ts4(R1)=ts4Lower;
    %     phi4(R1)=phi4Lower;
    %update field variables
    ts2Lower=(tau/Q2Scattering)*(sqrt(Q2Scattering*Q2Scattering+1)-1);
    te2Lower=(tau/Q2Scattering)*(sqrt(Q2Scattering*Q2Scattering+1)+1);
    phi2Lower=1/te2Lower-1/ts2Lower;
    ts4Lower=(tau/Q4Scattering)*(sqrt(Q4Scattering*Q4Scattering+1)-1);
    te4Lower=(tau/Q4Scattering)*(sqrt(Q4Scattering*Q4Scattering+1)+1);
    phi4Lower=1/te4Lower-1/ts4Lower;

end

%% Spatial inputs
dx                                  =   20;
dz                                  =   20;
origin                              =   [418100 4512100];

% vectors
X                                   =   origin(1) + (dx:dx:nodesX*dx);
Y                                   =   origin(2) + (dz:dz:nodesZ*dz);

% 2D space
[XX,YY]                             =   ndgrid(X,Y);

% Time and frequency inputs - see calculations at the start
dt                                  =   .0005; %s
nstep                               =   16000; %8 s
fSource                             =   1;% Dominant frequency

% Vector for times
Time                                =   0:dt:(nstep-1)*dt;
nframe                              =   100; % Visualizes every

% set anisotropic zone, anisotropy and anelasticity
interx=nodesX/2 - 0/dx;
nxt3=nodesX/2 - 0/dx;
nxt2=length(interx:nxt3);

interz=nodesZ/2-0/dz;
nzt3=nodesZ/2 - 0/dz;
nzt2=length(interz:nzt3);

% x0=nxt/2+600/dx;
% z0=nzt/2-4000/dz;
% xx=-nxt/2+1:nxt/2;
% zz=-nzt/2+1:nzt/2;
% R=zeros(nxt,nzt);
% for i=1:length(xx)
%     for j=1:length(zz)
%         R(i,j)=sqrt((xx(i)-x0)^2+(zz(j)-z0)^2);
%     end
% end
% rim=0/dx;
% R1=R>rim;
chosenCMap                      =   flipud(inferno);

figure
hn = surf(XX,YY,csi);
hn.EdgeColor = 'none';
hn.FaceColor = 'interp';
% oldcmap = colormap;
colormap(chosenCMap);
c = colorbar;
c.Label.String = 'Velocity fluctuations (km/s)';
view(0,90)
hold on

C = shaperead('COSs.shp');
mapshow(C,'FaceAlpha',0);

Faults= shaperead('FAULTS.shp');
mapshow(Faults,'Color', 'black')

axis equal
xlim([418000 434000])
ylim([4514000 4528000])

xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',20)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',20)
hold off


% absorbing boundaries parameters
r=1;
nab=120;
sab1(1:nab,1)=r.^(nab:-1:1);
sabtop=repmat(sab1,1,nodesZ);
sabbottom=flipud(repmat(sab1,1,nodesZ));
sabright=flipud(repmat(sab1,1,nodesX))';
sableft=repmat(sab1,1,nodesX)';

%% source/rec location, frequency of relaxation peaks and  energy
% 151007 0910 50.68 40-49.50  14-09.02  01.53    2.5 16  76   .5  .11   .2 
correlationLength=40+49.50/60;  b=14+9.02/60;
[Es,Ns,~]=deg2utm(correlationLength,b);
ix=floor((Es-origin(1))/dx);
iz=floor((Ns-origin(2))/dz);

% ix=floor((428377-origin(1))/dx);
% iz=floor((4520295-origin(2))/dz);

% GET cordinates stations
filename = 'stazioni_151007.txt';
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
stazioni=A.data;
namest=A.textdata;
correlationLength=stazioni(:,1)+stazioni(:,2)/60;
b=stazioni(:,3)+stazioni(:,4)/60;
latlong=[correlationLength,b];
[Er,Nr,Zone]=deg2utm(latlong(:,1),latlong(:,2));
rx = floor(Er-origin(1))/dx;
rz = floor(Nr-origin(2))/dz;


rhoSource=2500;
M=2.7;
% Me = (2/3) (log ES - 4.4)
E = 10^(4.4 + 1.5 * M)/fSource^2/rhoSource*dt^3; %Choy and Boatwright, 1995
Amp=sqrt(E);
% 1.58997E-06 (V/counts) / 1 * 1500 (V/m/s)
visibility='on';
rate=nstep/nframe;

% field variable
% u2: displacement
% e4 and e6: strain components
% s4 and s6: stress components
% e23 and e12: memory variables
u2=zeros(nodesX,nodesZ);
u=zeros(nodesX,nodesZ);
s4=zeros(nodesX,nodesZ);
s6=zeros(nodesX,nodesZ);
e23=zeros(nodesX,nodesZ);
e12=zeros(nodesX,nodesZ);

AmpTime1 = zeros(nstep,1);
AmpTime2 = zeros(nstep,1);
AmpTime3 = zeros(nstep,1);
AmpTime4 = zeros(nstep,1);
AmpTime5 = zeros(nstep,1);
AmpTime6 = zeros(nstep,1);
AmpTime7 = zeros(nstep,1);
AmpTime8 = zeros(nstep,1);
AmpTime9 = zeros(nstep,1);

% Set plot wavefield
h = figure('Name','Wavefield evolution',...
    'NumberTitle','off','Position',[1,200,900,800]);
h1=surf(XX,YY,u2);
colormap(chosenCMap)
hold on

if chooseModel==3
    plot([origin(1)+interx*dx origin(1)+interx*dx],...
        [origin(2)+interz*dz origin(2)+nzt3*dz],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+nxt3*dx origin(1)+nxt3*dx],...
        [origin(2)+interz*dz origin(2)+nzt3*dz],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+interx*dx origin(1)+nxt3*dx],...
        [origin(2)+interz*dz origin(2)+interz*dz],'LineWidth',...
        2,'Color',[.5 .5 .5]);
    plot([origin(1)+interx*dx origin(1)+nxt3*dx],...
        [origin(2)+nzt3*dz origin(2)+nzt3*dz],'LineWidth',...
        2,'Color',[.5 .5 .5]);  
end

light('Position',[0 0 0],'Style','infinite');
view(0,90)
lR=200;

for i=1:length(rx)
    rectangle('Position',[origin(1)+rx(i)*dx origin(2)+rz(i)*dz lR lR],...
        'FaceColor',[0 0 0])
    text(origin(1)+rx(i)*dx,origin(2)+rz(i)*dz+2*lR,namest{i},...
        'FontSize',12,'Color',[1 1 1])
end

viscircles([origin(1)+ix*dx origin(2)+iz*dz],400,'Color',[0 0 0])
% viscircles([origin(1)+x0*dx origin(2)+z0*dz],4000,'Color',[1 1 1])

mapshow(C,'FaceAlpha',0);
mapshow(Faults,'Color', 'black')

xlabel('WE (UTM/WGS84)','FontWeight','bold','FontSize',20)
ylabel('SN (UTM/WGS84)','FontWeight','bold','FontSize',20)

hold off

set(h1,'xdatasource','XX',...
    'ydatasource','YY',...
    'zdatasource','u2')
set(h1,'FaceLighting','phong','FaceColor','interp',...
    'AmbientStrength',1,'SpecularStrength',1,...
    'DiffuseStrength',0.9,'SpecularColorReflectance',1,...
    'SpecularExponent',10,'FaceAlpha',1, 'linestyle','none')

axis equal
xticks(418000:4000:434000)
xticklabels({'418000','422000','426000','430000','434000'})
yticks(4514000:4000:4526000)
yticklabels({'4514000','4518000','4522000','4526000'})
xlim([418000 434000])
ylim([4514000 4528000])
xlabel('WE [UTM/WGS84]'); ylabel('SN [UTM/WGS84]')
%colormap(jet)
colorbar

% Set plot wavefield
ha = figure('Name','Wavefield evolution',...
    'NumberTitle','off','Position',[1000,200,900,800]);

subplot(3,3,1)
% Set plot amplitudes
h2 = plot(Time,AmpTime1,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{1})
set(h2,'xdatasource','Time','ydatasource','AmpTime1','visible',visibility)

subplot(3,3,2)
% Set plot amplitudes
h3 = plot(Time,AmpTime2,'r','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{2})
set(h3,'xdatasource','Time','ydatasource','AmpTime2','visible',visibility)

subplot(3,3,3)
% Set plot amplitudes
h4 = plot(Time,AmpTime3,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{3})
set(h4,'xdatasource','Time','ydatasource','AmpTime3','visible',visibility)

subplot(3,3,4)
% Set plot amplitudes
h5 = plot(Time,AmpTime4,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{4})
set(h5,'xdatasource','Time','ydatasource','AmpTime4','visible',visibility)

subplot(3,3,5)
% Set plot amplitudes
h6 = plot(Time,AmpTime5,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{5})
set(h6,'xdatasource','Time','ydatasource','AmpTime5','visible',visibility)

subplot(3,3,6)
% Set plot amplitudes
h7 = plot(Time,AmpTime6,'r','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{6})
set(h7,'xdatasource','Time','ydatasource','AmpTime6','visible',visibility)

subplot(3,3,7)
% Set plot amplitudes
h8 = plot(Time,AmpTime7,'b','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{7})
set(h8,'xdatasource','Time','ydatasource','AmpTime7','visible',visibility)

subplot(3,3,8)
% Set plot amplitudes
h9 = plot(Time,AmpTime8,'r','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{8})
set(h9,'xdatasource','Time','ydatasource','AmpTime8','visible',visibility)

subplot(3,3,9)
% Set plot amplitudes
h10 = plot(Time,AmpTime9,'r','LineWidth',2);
xlabel('Time (s)','FontSize',12,'FontWeight','normal','Color','k')
ylabel('Amplitude ','FontSize',12,'FontWeight','normal','Color','k')
title(namest{9})
set(h10,'xdatasource','Time','ydatasource','AmpTime9','visible',visibility)


%% MODEL


% source wavelet - modified after errata
dt2=dt/2;
[f,nw2]=morlet1(fSource,dt2);
nw=nw2/2;
% figure
% plot(f)
% grid on
% title('Morlet Wavelet')

% Finite Differences weights
x1=9/(8*dx);
x2=-1/(24*dx);
z1=9/(8*dz);
z2=-1/(24*dz);

t=0;
index=0;
kx=3:nodesX-2;
kz=3:nodesZ-2;

%% TIME STEPPING
tic
for n=1:nstep
    if n==1000
        pause
    end
    t=t+dt;
    %
    % Absorbing boundaries
    
    % Horizontal stripes
    u2(1:nab,:)=u2(1:nab,:).*sabtop;
    u2(nodesX-nab+1:nodesX,:)=u2(nodesX-nab+1:nodesX,:).*sabbottom;
    
    % Vertical stripes
    u2(:,1:nab)=u2(:,1:nab).*sableft;
    u2(:,nodesZ-nab+1:nodesZ)=u2(:,nodesZ-nab+1:nodesZ).*sabright;
    
    %strains
    % Eqs. 9.24 and 9.53
    %i-3/2 -> i-2
    %i-1/2 -> i-1
    %i+1/2 -> i
    %i+3/2 -> i+1
    e4=z1*(u2(kx,kz)-u2(kx,kz-1))+z2*(u2(kx,kz+1)-u2(kx,kz-2));
    e6=x1*(u2(kx,kz)-u2(kx-1,kz))+x2*(u2(kx+1,kz)-u2(kx-2,kz));
    %memory variables
    f1=2*ts2(kx,kz)-dt;
    f2=2*ts2(kx,kz)+dt;
    ee=e23(kx,kz);
    
    %Eqs. (4.232)4 and 9.54
    e23(kx,kz)=(2*dt*ts2(kx,kz).*phi2(kx,kz).*e4+f1.*e23(kx,kz))./f2;
    e23(kx,kz)=0.5*(e23(kx,kz)+ee);
    f1=2*ts4(kx,kz)-dt;
    f2=2*ts4(kx,kz)+dt;
    ee=e12(kx,kz);
    
    %Eqs. (4.149)6 and 9.54
    e12(kx,kz)=(2*dt*ts4(kx,kz).*phi4(kx,kz).*e6+f1.*e12(kx,kz))./f2;
    e12(kx,kz)=0.5*(e12(kx,kz)+ee);
    % stresses
    %Eq. 4.2333
    s4(kx,kz)=c44(kx,kz).*(e4+e23(kx,kz))+c46(kx,kz).*e6;
    s6(kx,kz)=c66(kx,kz).*(e6+e12(kx,kz))+c46(kx,kz).*e4;
    
    %strains
    % Eqs. 9.24
    %i-3/2 -> i-2
    %i-1/2 -> i-1
    %i+1/2 -> i
    %i+3/2 -> i+1
    ds4=z1*(s4(kx,kz+1)-s4(kx,kz))+z2*(s4(kx,kz+2)-s4(kx,kz-1));
    ds6=x1*(s6(kx+1,kz)-s6(kx,kz))+x2*(s6(kx+2,kz)-s6(kx-1,kz));
    
    %acceleration
    acc=(ds4+ds6)./rho(kx,kz);
    
    %source - corrected errata
    if n<=nw
        n1=n*2-1;
        source=Amp*f(n1);
    elseif n>nw
        source =0;
    end
    
    % Euler equation
    % Eqs. (1.46)1 and 9.52
    
    u(kx,kz)=2*u2(kx,kz)-u(kx,kz)+dt*dt*acc;
    u(ix,iz)=u(ix,iz)+source;
    
    
    % Update displacement
    uu=u2(kx,kz);
    u2(kx,kz)=u(kx,kz);
    u(kx,kz)=uu;
    
    %--------------------
    AmpTime1(n)   =   u2(floor(rx(1)),floor(rz(1)));%/sqrt(((ix-rx(1))^2+(iz-rz(1))^2)*10)/2;
    AmpTime2(n)   =   u2(floor(rx(2)),floor(rz(2)));%/sqrt(((ix-rx(2))^2+(iz-rz(2))^2)*10)/2;
    AmpTime3(n)   =   u2(floor(rx(3)),floor(rz(3)));%/sqrt(((ix-rx(3))^2+(iz-rz(3))^2)*10)/2;
    AmpTime4(n)   =   u2(floor(rx(4)),floor(rz(4)));%/sqrt(((ix-rx(4))^2+(iz-rz(4))^2)*10)/2;
    AmpTime5(n)   =   u2(floor(rx(5)),floor(rz(5)));%/sqrt(((ix-rx(5))^2+(iz-rz(5))^2)*10)/2;
    AmpTime6(n)   =   u2(floor(rx(6)),floor(rz(6)));%/sqrt(((ix-rx(6))^2+(iz-rz(6))^2)*10)/2;
    AmpTime7(n)   =   u2(floor(rx(7)),floor(rz(7)));%/sqrt(((ix-rx(7))^2+(iz-rz(7))^2)*10)/2;
    AmpTime8(n)   =   u2(floor(rx(8)),floor(rz(8)));%/sqrt(((ix-rx(8))^2+(iz-rz(8))^2)*10)/2;
    AmpTime9(n)   =   u2(floor(rx(9)),floor(rz(9)));%/sqrt(((ix-rx(9))^2+(iz-rz(9))^2)*10)/2;
    
    
    % Plot
    if isequal(mod(n,nframe),0)
        index=index+1;
        refreshdata(h1)
        refreshdata(h2)
        refreshdata(h3)
        refreshdata(h4)
        refreshdata(h5)
        refreshdata(h6)
        refreshdata(h7)
        refreshdata(h8)
        refreshdata(h9)
        refreshdata(h10)
        drawnow
        F(index) = getframe(gcf);
        
    end
    
end
toc

if chooseModel==1
    label='Isotropic_1Hz.avi';
    label1='Isotropic_1Hz.fig';
elseif chooseModel==2
    label='Anisotropic_1Hz.avi';
    label1='Anisotropic_1Hz.fig';
elseif chooseModel==3
    label='Anisotropic_scattering.avi';
    label1='Anisotropic_scattering.fig';
end

video = VideoWriter(label,'Motion JPEG AVI' );

video.FrameRate=3;

open(video);
writeVideo(video,F);
close(video);
saveas(h,label1)
saveas(ha,strcat(num2str(fSource),'_Modelled'))

%%
fSource=1;
l11=50.68;
l22=58.68;

nmarker0='151007/*N.SAC';

list                        = dir(nmarker0);
filenames                   = {list.name}'; %create cell array of file names
filedir                     = {list.folder}'; %create cell array of file folder
listasac                    = strcat(filedir,'/',filenames);
lls                         = length(listasac);
[tempis,sismaew,SAChdr]          =   fget_sac(listasac{1}); % Imports SAC files
lsis=length(sismaew);

title1=cell(lls,1);
title2=cell(lls,1);

% GET cordinates stations
filename = 'stazioni_151007.txt';
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
stazioni=A.data;
namest=A.textdata;

nmarker=nmarker0(9);
% for i=1:lls
%     title1(i)=strcat(namest(i), nmarker);
% end


% sisma=zeros(lls,lsis);
% figure('Name','Seismograms',...
%     'NumberTitle','off','Position',[100,200,1600,800]);

% Yl=3e-05;
% for i=1:lls
%     [tempis,sisma1,SAChdr]=fget_sac(listasac{i}); % Imports SAC files
%     srate                       =   1/SAChdr.times.delta; %sampling f0uency
%     sisma2  =   detrend(sisma1/SAChdr.data.scale,1);%remove trend
%     sisma(i,:) = sisma2;
%     subplot(3,3,i)
%     plot(tempis,sisma(i,:),'k','LineWidth',2);
%     title(title1{i})
%     xlim([l11,l22]);
%     
%     %     ylim([-Yl,Yl]);
%     
%     xlabel('Time (s)'); ylabel('Meters')
% end

tu                      =   tukeywin(lsis,0.05);%avoid windowing
window=20;

nmarker2=[nmarker0(9) '- 1 Hz'];
for i=1:lls
    title2(i)=strcat(namest(i), nmarker2);
end

% onset=1;
% Yl1=1e-05;
% Yl2=1e-05;

hb=figure('Name','Filtered - 1 Hz',...
    'NumberTitle','off','Position',[1,200,900,800]);


for i=1:lls
    [tempis,sisma1,SAChdr]=fget_sac(listasac{i}); % Imports SAC files
    srate                       =   1/SAChdr.times.delta; %sampling f0uency
    Wn                      =   ([fSource-fSource/3 fSource+fSource/3]/srate*2); %frequency band
    [z,p,k]                 =   butter(4,Wn,'bandpass'); %butter filter
    [sos,g]                 =   zp2sos(z,p,k); % Convert to SOS form
    sisma2  =   detrend(sisma1/SAChdr.data.scale,1);%remove trend
    sisma(i,:) = sisma2;
    
    tsisma         =   tu'.*sisma(i,:);
    fsisma         =   filtfilt(sos,g,tsisma);% filtered waveform
    subplot(3,3,i)
    switch i
        case {1, 3, 4, 5, 7}
            plot(tempis-50,fsisma,'b','LineWidth',2);
        case {2, 6, 8, 9}
            plot(tempis-50,fsisma,'r','LineWidth',2);
    end
    title(title2{i})
    xlim([l11,l22]);
    xlabel('Time (s)'); ylabel('Meters')
end

saveas(hb,strcat(num2str(fSource),'_Data'))
