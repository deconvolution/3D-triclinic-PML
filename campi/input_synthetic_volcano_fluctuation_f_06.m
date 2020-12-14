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
close all;
clear all;
%%
M=table2array(readtable('./modvPS.txt'));

X0=unique(M(:,1));
Y0=unique(M(:,2));
Z0=unique(M(:,3));

nx=length(X0);
ny=length(Y0);
nz=length(Z0);

X=permute(reshape(M(:,1),[nz,ny,nx]),[3,2,1]);
Y=permute(reshape(M(:,2),[nz,ny,nx]),[3,2,1]);
Z=permute(reshape(M(:,3),[nz,ny,nx]),[3,2,1]);

vp=permute(reshape(M(:,4),[nz,ny,nx]),[3,2,1]);
vs=permute(reshape(M(:,5),[nz,ny,nx]),[3,2,1]);

tit={'WE','SN','Al','vp','vs','vp/vs'};
vtkwrite('campi_vp.vtk','structured_grid',X,Y,Z,'scalars','campi_vp',vp);
vtkwrite('campi_vs.vtk','structured_grid',X,Y,Z,'scalars','campi_vs',vs);
%% Dimensions
dx=100;
dy=100;
dz=100;

nt=1000;
ns=nt;
dt=10^-2;
%% compute stiffness
mu=1000*vs.^2;
lambda=1000*vp.^2-2*mu;
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
C.rho=ones(nx,ny,nz)*1000;
C.C11=lambda+2*mu;
C.C12=lambda;
C.C13=lambda;
C.C14=0;
C.C15=0;
C.C16=0;
C.C22=lambda+2*mu;
C.C23=lambda;
C.C24=0;
C.C25=0;
C.C26=0;
C.C33=lambda+2*mu;
C.C34=0;
C.C35=0;
C.C36=0;
C.C44=mu;
C.C45=0;
C.C46=0;
C.C55=mu;
C.C56=0;
C.C66=mu;
C.lambda2=0;
C.mu2=0;
%% Source and source signals
% source locations
s1=[fix(nx/2)];
s2=[fix(ny/2)];
s3=130;

% source frequency
freq=.6;

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
% receiver location, find topography
r1=[1:5:nx,ones(1,length(1:5:ny))*fix(nx/2)];
r2=[ones(size(1:5:nx))*fix(ny/2),1:5:ny];
r3=ones(size(r1));
for i=1:length(r1)
    tt=C.C44(r1(i),r2(i),:);
    tt2=find(tt~=0);
    tt3=min(tt2)+2;
    r3(i)=tt3;
end
%% PML
% PML layers
lp=15;

% theoretical reflection coeficcient
Rc=.0001;

% PML power, normally 2
nPML=2;
%% plot
% point interval in time steps
plot_interval=20;

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
[R1,R2,R3,v1,v2,v3,E]=triclinic_3D_themas(dt,dx,dy,dz,nt,nx,ny,nz,C,...
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
%% plot seismogram
figure;
tt=find(r1==96&r2==96);
ax=plot(dt:dt:dt*nt,R1(tt(1),:),'red');
hold on;
tt=find(r1==21&r2==96);
ax2=plot(dt:dt:dt*nt,R1(tt(1),:),'blue');
hold on;
tt=find(r1==171&r2==96);
ax3=plot(dt:dt:dt*nt,R1(tt(1),:),'green');
hold on;
tt=find(r1==96&r2==21);
ax4=plot(dt:dt:dt*nt,R1(tt(1),:),'black');
hold on;
tt=find(r1==96&r2==171);
ax5=plot(dt:dt:dt*nt,R1(tt(1),:),'cyan');
hold on;
legend([ax,ax2,ax3,ax4,ax5], ...
    'center','x-','x+','y-','y+');
xlabel('t [s]');
ylabel('v_1 [m/s]');

figure;
tt=find(r1==96&r2==96);
ax=plot(dt:dt:dt*nt,R2(tt(1),:),'red');
hold on;
tt=find(r1==21&r2==96);
ax2=plot(dt:dt:dt*nt,R2(tt(1),:),'blue');
hold on;
tt=find(r1==171&r2==96);
ax3=plot(dt:dt:dt*nt,R2(tt(1),:),'green');
hold on;
tt=find(r1==96&r2==21);
ax4=plot(dt:dt:dt*nt,R2(tt(1),:),'black');
hold on;
tt=find(r1==96&r2==171);
ax5=plot(dt:dt:dt*nt,R2(tt(1),:),'cyan');
hold on;
legend([ax,ax2,ax3,ax4,ax5], ...
    'center','x-','x+','y-','y+');
xlabel('t [s]');
ylabel('v_2 [m/s]');

figure;
tt=find(r1==96&r2==96);
ax=plot(dt:dt:dt*nt,R3(tt(1),:),'red');
hold on;
tt=find(r1==21&r2==96);
ax2=plot(dt:dt:dt*nt,R3(tt(1),:),'blue');
hold on;
tt=find(r1==171&r2==96);
ax3=plot(dt:dt:dt*nt,R3(tt(1),:),'green');
hold on;
tt=find(r1==96&r2==21);
ax4=plot(dt:dt:dt*nt,R3(tt(1),:),'black');
hold on;
tt=find(r1==96&r2==171);
ax5=plot(dt:dt:dt*nt,R3(tt(1),:),'cyan');
hold on;
legend([ax,ax2,ax3,ax4,ax5], ...
    'center','x-','x+','y-','y+');
xlabel('t [s]');
ylabel('v_3 [m/s]');