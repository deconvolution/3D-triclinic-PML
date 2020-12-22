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
load('./campi_model.mat');
nx=campi_model.nx;
ny=campi_model.ny;
nz=campi_model.nz;

X=campi_model.X;
Y=campi_model.Y;
Z=campi_model.Z;
%% Dimensions
dx=campi_model.dx;
dy=campi_model.dy;
dz=campi_model.dz;

nt=6000;
ns=nt;
dt=10^-2/4;
%% initialize density and stiffness
C=campi_model.C;
air_boundary=18;
C.C11(:,:,1:air_boundary)=0;
C.C12(:,:,1:air_boundary)=0;
C.C13(:,:,1:air_boundary)=0;
C.C22(:,:,1:air_boundary)=0;
C.C23(:,:,1:air_boundary)=0;
C.C33(:,:,1:air_boundary)=0;
C.C44(:,:,1:air_boundary)=0;
C.C55(:,:,1:air_boundary)=0;
C.C66(:,:,1:air_boundary)=0;

C.C11=C.C11*2/3;
C.C55=C.C55*2/3;
C.C66=C.C66*2/3;
%% Source and source signals
% source locations
s1=campi_model.s1;
s2=campi_model.s2;
s3=campi_model.s3;

% source frequency
freq=1;

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
r1=campi_model.r1;
r2=campi_model.r2;
r3=campi_model.r3;
%% PML
% PML layers
lp=15;

% theoretical reflection coeficcient
Rc=.0001;

% PML power, normally 2
nPML=2;

%
PMLxp=1;
PMLxm=1;
PMLyp=1;
PMLym=1;
PMLzp=1;
PMLzm=1;
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
    lp,nPML,Rc,PMLxp,PMLxm,PMLyp,PMLym,PMLzp,PMLzm, ...
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
[rec,simu_info]=rec_conversion(X,Y,Z,s1,s2,s3,r1,r2,r3,R1,R2,R3,dt,dx,dy,dz,nt,nx,ny,nz,lp,ones(1,length(s1))*freq);
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
for l3=1:size(R3,1)
    subplot(3,3,l3)
    plot(dt:dt:dt*nt,R1(l3,:));
    xlabel('t [s]');
    ylabel('v_1 [m/s]');
end
print(gcf,[path '/rec_v1' '.png'],'-dpng','-r200');

figure;
for l3=1:size(R3,1)
    subplot(3,3,l3)
    plot(dt:dt:dt*nt,R2(l3,:));
    xlabel('t [s]');
    ylabel('v_2 [m/s]');
end
print(gcf,[path '/rec_v2' '.png'],'-dpng','-r200');

figure;
for l3=1:size(R3,1)
    subplot(3,3,l3)
    plot(dt:dt:dt*nt,R3(l3,:));
    xlabel('t [s]');
    ylabel('v_3 [m/s]');
end
print(gcf,[path '/rec_v3' '.png'],'-dpng','-r200');

figure;
for l3=1:size(R3,1)
    subplot(3,3,l3)
    plot(dt:dt:dt*nt,R1(l3,:)./R3(l3,:));
    xlabel('t [s]');
    ylabel('v_1/v_3');
end
print(gcf,[path '/rec_v1_v3' '.png'],'-dpng','-r200');

figure;
for l3=1:size(R3,1)
    subplot(3,3,l3)
    plot(dt:dt:dt*nt,R2(l3,:)./R3(l3,:));
    xlabel('t [s]');
    ylabel('v_2/v_3');
end
print(gcf,[path '/rec_v2_v3' '.png'],'-dpng','-r200');