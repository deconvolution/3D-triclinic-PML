%% Dimensions
clear all;
close all;

dx=10;
dy=10;
dz=10;

nx=101;
ny=101;
nz=101;

nt=1000;
ns=nt;
dt=10^-3;
%% grid specifying 3D coordinates
[Y,X,Z]=meshgrid((1:ny)*dy,(1:nx)*dx,(1:nz)*dz);
%% initialize density and stiffness
C.rho=zeros(nx,ny,nz);
C.C11=zeros(nx,ny,nz);
C.C12=C.C11;
C.C13=C.C11;
C.C14=C.C11;
C.C15=C.C11;
C.C16=C.C11;
C.C22=C.C11;
C.C23=C.C11;
C.C24=C.C11;
C.C25=C.C11;
C.C26=C.C11;
C.C33=C.C11;
C.C34=C.C11;
C.C35=C.C11;
C.C36=C.C11;
C.C44=C.C11;
C.C45=C.C11;
C.C46=C.C11;
C.C55=C.C11;
C.C56=C.C11;
C.C66=C.C11;
%% transfer stiffness, density and viscosity
lambda=10^9;
mu=10^9;
C.C11(:)=lambda+2*mu;
C.C12(:)=lambda;
C.C13(:)=lambda;
C.C22(:)=lambda+2*mu;
C.C23(:)=lambda;
C.C33(:)=lambda+2*mu;
C.C44(:)=mu;
C.C55(:)=mu;
C.C66(:)=mu;
C.rho(:)=1000;
C.lambda2=0;
C.mu2=0;
%% Source and source signals
% source locations
s1=[fix(nx/2)];
s2=[fix(ny/2)];
s3=[75];

% source frequency
freq=6;

% magnitude
M=2;

% source signal generation
singles=rickerWave(freq,dt,ns,M);

% assign signal to coresponding source
src1=zeros(nt,1);
src2=src1;
src3=src1;

src1=1*[singles];
src2=1*[singles];
src3=1*[singles];

% souce type, 'P' for P-wave source, 'S' for S-wave source, 'D' for
% directional source
source_type=['D'];
%% Receiver
% receiver location
r1=[21:3:101-20];
r2=[ones(size(r1))*fix(ny/2)];
r3=[ones(size(r1))*31];
%% PML
% PML layers
lp=20;

% theoretical reflection coeficcient
Rc=.0001;

% PML power, normally 2
nPML=2;
%% plot
% point interval in time steps
plot_interval=50;

% figure path
p2=mfilename('fullpath');
path=[p2 '/'];

% a slice function is used
% view angle
view_angle=[];
% x slice
x2=s1;
% y slice
y2=s2;
% z slice
z2=s3;
%% solving wavefield
[R1,R2,R3,v1,v2,v3,E]=triclinic_3D(dt,dx,dy,dz,nt,nx,ny,nz,C,...
    s1,s2,s3,src1,src2,src3,source_type, ...
    r1,r2,r3, ...
    lp,nPML,Rc, ...
    X,Y,Z,...
    x2,y2,z2,plot_interval,view_angle,...
    path);
%% plot energy
figure
plot(dt:dt:dt*nt,E);
xlabel('t [s]');
ylabel('E [J]');
print(gcf,[path 'E.png'],'-dpng','-r200');
%% recording saving
[rec,simu_info]=rec_conversion(s1,s2,s3,r1,r2,r3,R1,R2,R3,dt,dx,dy,dz,nt,nx,ny,nz,lp,ones(1,length(s1))*freq);
save([path 'rec.mat'],'rec');
save([path 'simu_info.mat'],'rec');
figure('name','rec v3 [m/s]');
imagesc(1:length(r3),dt:dt:dt*nt,R3');
xlabel('Nr');
ylabel('t [s]');
title('v3 [m/s]');
%% make gif
sources=[path '/pic/'];
delaytime=.2;
filename='animation';
gifmaker(filename,delaytime,sources);