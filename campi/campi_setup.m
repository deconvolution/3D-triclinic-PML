%% model set up for event
%% Dimensions
close all;
clear all;
%%
M=table2array(readtable('./modvPS.txt'));

X0=unique(M(:,1));
Y0=flip(unique(M(:,2)),1);
Z0=unique(M(:,3));

nx=length(X0);
ny=length(Y0);
nz=length(Z0);


X=flip(permute(reshape(M(:,1),[nz,ny,nx]),[3,2,1]),2);
Y=flip(permute(reshape(M(:,2),[nz,ny,nx]),[3,2,1]),2);
Z=flip(permute(reshape(M(:,3),[nz,ny,nx]),[3,2,1]),2);

vp=flip(permute(reshape(M(:,4),[nz,ny,nx]),[3,2,1]),2);
vs=flip(permute(reshape(M(:,5),[nz,ny,nx]),[3,2,1]),2);

tit={'WE','SN','Al','vp','vs','vp/vs'};
vtkwrite('campi_vp.vtk','structured_grid',X,Y,Z,'scalars','campi_vp',vp);
vtkwrite('campi_vs.vtk','structured_grid',X,Y,Z,'scalars','campi_vs',vs);
%% Dimensions
dx=100;
dy=100;
dz=100;
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
%% receiver location
filename = 'stazioni_151007.txt';
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename,delimiterIn,headerlinesIn);
stazioni=A.data;
namest=A.textdata;
coor=stazioni(:,1)+stazioni(:,2)/60;
b=stazioni(:,3)+stazioni(:,4)/60;
latlong=[coor,b];
[Er,Nr,Zone]=deg2utm(latlong(:,1),latlong(:,2));
r1=fix((Er'-min(X0))/dx);
r2=fix((Nr'-min(Y0))/dy);
r3=fix((max(Z0)-stazioni(:,5)'*1000)/dz);
%% Source and source signals
[Es,Ns,~]=deg2utm(40+49.50/60,14+9.02/60);
s1=fix((Es'-min(X0))/dx);
s2=fix((Ns'-min(Y0))/dy);
s3=fix((max(Z0)-(-1.53)*1000)/dz);
%% save model setup
campi_model.X=X;
campi_model.Y=Y;
campi_model.Z=Z;
campi_model.nx=nx;
campi_model.ny=ny;
campi_model.nz=nz;
campi_model.dx=dx;
campi_model.dy=dy;
campi_model.dz=dz;
campi_model.C=C;
campi_model.s1=s1;
campi_model.s2=s2;
campi_model.s3=s3;
campi_model.r1=r1;
campi_model.r2=r2;
campi_model.r3=r3;
save('campi_model.mat','campi_model');
