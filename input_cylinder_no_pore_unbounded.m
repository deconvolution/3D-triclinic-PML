%% Introduction
%{
This is the code to simulate triclinic viscoelastic acoustic-elastic
velocity field. The viscoelastic model belongs to Kelvin-Voigt type. PML is
implemented. This script provides all neccessary parameters to the solver
'triclinic_3D.m'.

Finite difference method:
1st order in time,
2nd order in space.

Coordinate convention:
In 3D array: 1-x, 2-y, 3-z.
In 3D plot: 1-y, 2-x, 3-z.
%}
%% Dimensions
clear all;
close all;

dx=.1;
dy=.1;
dz=.1;

nx=101;
ny=101;
nz=101;

nt=3000;
ns=nt;
dt=10^-5;
%% construct rock sample
S=ones(nx,ny,nz);

% form cylinder
%{
cx=50;
cy=50;
tt=zeros(nx,ny);
[tty,ttx]=meshgrid(1:ny,1:nx);
IND=find((ttx-cx).^2+(tty-cy).^2<30^2);
tt(IND)=1;
S(:,:,11:80)=repmat(tt,[1,1,80-11+1]);
%}

% form sphere pore
cx=35;
cy=50;
cz=50;
[tty,ttx,ttz]=meshgrid(1:ny,1:nx,1:nz);
IND=find((ttx-cx).^2+(tty-cy).^2+(ttz-cz).^2<10^2);
S(IND)=0;
S(:)=1;
IND_solid=find(S==1);
IND_fluid=find(S==0);
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
mu=2.94*10^3*(1490)^2;
lambda=2.94*10^3*(2860)^2-2*mu;
C.C11(IND_solid)=lambda+2*mu;
C.C12(IND_solid)=lambda;
C.C13(IND_solid)=lambda;
C.C22(IND_solid)=lambda+2*mu;
C.C23(IND_solid)=lambda;
C.C33(IND_solid)=lambda+2*mu;
C.C44(IND_solid)=mu;
C.C55(IND_solid)=mu;
C.C66(IND_solid)=mu;
C.rho(IND_solid)=2.94*10^3;


lambdaf=1225*340^2;
C.rho(IND_fluid)=1225;
C.C11(IND_fluid)=lambdaf;
C.C12(IND_fluid)=lambdaf;
C.C13(IND_fluid)=lambdaf;
C.C22(IND_fluid)=lambdaf;
C.C23(IND_fluid)=lambdaf;
C.C33(IND_fluid)=lambdaf;

C.lambda2=0;
C.mu2=0;
%% Source and source signals
% source locations
s1=[fix(nx/2)];
s2=[fix(ny/2)];
s3=[75];

% source frequency
freq=600;

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
source_type=['S'];
%% Receiver
% receiver location
r1=1:3:101;
r2=ones(size(r1))*fix(ny/2);
r3=ones(size(r1))*21;
%% PML
% PML layers
lp=20;

% theoretical reflection coeficcient
Rc=.001;

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