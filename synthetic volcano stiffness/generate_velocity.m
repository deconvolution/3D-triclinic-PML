% All of the parameters are set SI units in this script
%% read modulus
clear all;
close all;

load('ShearMod.mat');

X=Data.X;
Y=Data.Y;
Z=Data.Z;
G=Data.G;

[nx,ny,nz]=size(G);

X_min=min(X,[],'all');
X_max=max(X,[],'all');

Y_min=min(Y,[],'all');
Y_max=max(Y,[],'all');

Z_min=min(Z,[],'all');
Z_max=max(Z,[],'all');

dx=(X_max-X_min)/(nx)*1000;
dy=(Y_max-Y_min)/(ny)*1000;
dz=(Z_max-Z_min)/(nz)*1000;

air_ind=find(G==0);
nu=.2;
G=flip(G,3);
X=flip(X,3);
Y=flip(Y,3);
Z=flip(Z,3);
tt=G(:,:,1:55);
tt(tt>0 & tt<2*10^10)=0;
G(:,:,1:55)=tt;
%% calculate lambda and mu
lambda=2*G.*nu./(1-2*nu);
mu=G;
rho=ones(nx,ny,nz)*10^3;
vtkwrite('G.vtk','structured_grid',X,Y,Z,'scalars','G',G);
%%
vp=sqrt((lambda+2*mu)./rho);
vp(isnan(vp))=0;
vs=sqrt(mu./rho);
%% add velocity fluctuations
% 3 layers
%{
N=-.2;
% layer 1
layer=46:53;
kappa=200;
a=1000;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vp(:,:,layer)=vp(:,:,layer)+tt2;
kappa=100;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vs(:,:,layer)=vs(:,:,layer)+tt2;
%}
%{
% layer 2
layer=28:35;
kappa=300;
a=1500;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vp(:,:,layer)=vp(:,:,layer)+tt2;
kappa=150;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vs(:,:,layer)=vs(:,:,layer)+tt2;


% layer 3
layer=40:65;
kappa=400;
a=2000;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vp(:,:,layer)=vp(:,:,layer)+tt2;
kappa=200;
tt2=random_media_creation(nx,ny,length(layer),dx,dy,dz,kappa,a,N);
vs(:,:,layer)=vs(:,:,layer)+tt2;
%}
%% add vel gradient
z1=40;
z2=154;
g=.0262;
vp=add_vel_gradient(nx,ny,nz,dx,dy,dz,vp,z1,z2,g);
g=.0262/1.72;
vs=add_vel_gradient(nx,ny,nz,dx,dy,dz,vs,z1,z2,g);
figure
subplot(1,2,1)
plot(reshape(vp(fix(nx/2),fix(ny/2),:),[nz,1]),(1:nz)*dz);
xlabel('vp [m/s]');
ylabel(['z*' num2str(dz) '[m]']);
set(gca,'ydir','reverse');
subplot(1,2,2)
plot(reshape(vs(fix(nx/2),fix(ny/2),:),[nz,1]),(1:nz)*dz);
xlabel('vs [m/s]');
ylabel(['z*' num2str(dz) '[m]']);
set(gca,'ydir','reverse');
%% transform to C
C=zeros([6,6,size(vp)]);
C(1,1,:,:,:)=rho.*vp.^2;
C(1,2,:,:,:)=rho.*vp.^2-2*rho.*vs.^2;
C(1,3,:,:,:)=rho.*vp.^2-2*rho.*vs.^2;
C(2,2,:,:,:)=rho.*vp.^2;
C(2,3,:,:,:)=rho.*vp.^2-2*rho.*vs.^2;
C(3,3,:,:,:)=rho.*vp.^2;
C(4,4,:,:,:)=rho.*vs.^2;
C(5,5,:,:,:)=rho.*vs.^2;
C(6,6,:,:,:)=rho.*vs.^2;
%%
mean(sqrt(reshape(C(1,1,:,:,:),[nx,ny,nz])./rho),'all')
%%
lambda=reshape(C(1,3,:,:,:),[nx,ny,nz]);
G=reshape(C(4,4,:,:,:),[nx,ny,nz]);
%% fracture 1
% 1 set fracture above the magma body
% fracture normal parallel to 1 axis
f_x=10:110;
f_y=10:110;
f_z=41:45;
e=zeros(nx,ny,nz);
e(f_x,f_y,f_z)=.1;
theta=0;
C=add_fracture(lambda,G,e,theta);
%% fracture 2
% 2 sets of fractures
% fracture normal parallel to 1 axis
% fracture normal parallel to pi/6 relative to 1 axis (right hand rule)
f_x2=10:110;
f_y2=10:110;
f_z2=56:60;
e=zeros(nx,ny,nz);
e(f_x2,f_y2,f_z2)=.1;
theta2=pi/6;
C2=add_fracture2(lambda,G,e,theta,theta2);
%%
C(:,:,f_x2,f_y2,f_z2)=C2(:,:,f_x2,f_y2,f_z2);
%%
figure('name','C');
subplot(1,2,1)
slice(reshape(C(1,1,:,:,:),[nx,ny,nz]),fix(nx/2),fix(ny/2),fix(nz/2));
set(gca,'zdir','reverse');
xlabel(['y*' num2str(dy) '[m]']);
ylabel(['x*' num2str(dx) '[m]']);
zlabel(['z*' num2str(dz) '[m]']);
colorbar;
title('C11 [Pa]');
subplot(1,2,2)
slice(reshape(C(1,6,:,:,:),[nx,ny,nz]),40,40,fix(nz/2));
set(gca,'zdir','reverse');
xlabel(['y*' num2str(dy) '[m]']);
ylabel(['x*' num2str(dx) '[m]']);
zlabel(['z*' num2str(dz) '[m]']);
title('C16 [Pa]');
colorbar;
shg;
hold off;
%% assign air
for i=1:3
    for j=i:3
        tt=C(i,j,:,:,:);
        tt(tt==0)=rho(tt==0)*340^2;
        C(i,j,:,:,:)=tt;
    end
end
%%
c.X=X;
c.Y=Y;
c.Z=Z;
c.dx=dx;
c.dy=dy;
c.dz=dz;
c.C11=reshape(C(1,1,:,:,:),[nx,ny,nz]);
c.C12=reshape(C(1,2,:,:,:),[nx,ny,nz]);
c.C13=reshape(C(1,3,:,:,:),[nx,ny,nz]);
c.C14=reshape(C(1,4,:,:,:),[nx,ny,nz]);
c.C15=reshape(C(1,5,:,:,:),[nx,ny,nz]);
c.C16=reshape(C(1,6,:,:,:),[nx,ny,nz]);
c.C22=reshape(C(2,2,:,:,:),[nx,ny,nz]);
c.C23=reshape(C(2,3,:,:,:),[nx,ny,nz]);
c.C24=reshape(C(2,4,:,:,:),[nx,ny,nz]);
c.C25=reshape(C(2,5,:,:,:),[nx,ny,nz]);
c.C26=reshape(C(2,6,:,:,:),[nx,ny,nz]);
c.C33=reshape(C(3,3,:,:,:),[nx,ny,nz]);
c.C34=reshape(C(3,4,:,:,:),[nx,ny,nz]);
c.C35=reshape(C(3,5,:,:,:),[nx,ny,nz]);
c.C36=reshape(C(3,6,:,:,:),[nx,ny,nz]);
c.C44=reshape(C(4,4,:,:,:),[nx,ny,nz]);
c.C45=reshape(C(4,5,:,:,:),[nx,ny,nz]);
c.C46=reshape(C(4,6,:,:,:),[nx,ny,nz]);
c.C55=reshape(C(5,5,:,:,:),[nx,ny,nz]);
c.C56=reshape(C(5,6,:,:,:),[nx,ny,nz]);
c.C66=reshape(C(6,6,:,:,:),[nx,ny,nz]);
c.rho=rho;
save('c.mat','c');
%%
plot(reshape(c.C33(40,40,:),[nz,1]))
%%
c.C11(find(c.C11<10^10 & c.C11>1000*340^2))
%%
figure;
imagesc(reshape(G(100,:,:),[ny,nz]));