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
%%
dx=100;
dy=100;
dz=100;
%% PML input
lp=20;
lpn=2;
Rc=10^-3;
%% compute C0
lambda=reshape(M3(:,:,:,4).^2,[1,1,nx0,ny0,nz0]);
mu=reshape(M3(:,:,:,5).^2,[1,1,nx0,ny0,nz0]);
C0=zeros(6,6,nx0,ny0,nz0);
C0(1,1,:,:,:)=lambda;
C0(2,2,:,:,:)=lambda;
C0(3,3,:,:,:)=lambda;
C0(4,4,:,:,:)=mu;
C0(5,5,:,:,:)=mu;
C0(6,6,:,:,:)=mu;
C0(1,2,:,:,:)=lambda-2*mu;
C0(1,3,:,:,:)=lambda-2*mu;
C0(2,3,:,:,:)=lambda-2*mu;
%% dimensions
dx=12;
dy=12;
dz=12;

nx2=round(nx0/dx*dx0);
ny2=round(ny0/dy*dy0);
nz2=round(nz0/dz*dz0);

C=zeros(6,6,nx2+2*lp+2,ny2+2*lp+2,nz2+2*lp+2);
[~,~,nx,ny,nz]=size(C);

for i=1:6
    for j=1:6
        C(i,j,lp+2:nx-lp-1,lp+2:ny-lp-1,lp+2:nz-lp-1)=reshape(imresize3(reshape(C0(i,j,:,:,:),[nx0,ny0,nz0]),[nx2,ny2,nz2]),[1,1,nx2,ny2,nz2]);
    end
end

C(:,:,1:lp+1,:,:)=repmat(C(:,:,lp+2,:,:),[1,1,lp+1,1,1]);
C(:,:,nx-lp:nx,:,:)=repmat(C(:,:,nx-lp-1,:,:),[1,1,lp+1,1,1]);

C(:,:,:,1:lp+1,:)=repmat(C(:,:,:,lp+2,:),[1,1,1,lp+1,1]);
C(:,:,:,ny-lp:ny,:)=repmat(C(:,:,:,ny-lp-1,:),[1,1,1,lp+1,1]);

C(:,:,:,:,1:lp+1)=repmat(C(:,:,:,:,lp+2),[1,1,1,1,lp+1]);
C(:,:,:,:,nz-lp:nz)=repmat(C(:,:,:,:,nz-lp-1),[1,1,1,1,lp+1]);
%%
figure
pcolor3(reshape(C(3,3,:,:,:),[nx,ny,nz]));
xlabel('x');
ylabel('y');
zlabel('z');
colorbar;
title('C33');
set(gca,'zdir','reverse');
shg;
%% time step
dt=10^-3; % [s]
nt=8000; % Amount of time steps
ns=nt;
%% Define viscoelastic parameters
theta=0;
tau=5/8*theta;

Eta=zeros(6,6,nx,ny,nz);
%{
lambda2=lambda*theta+2/3*mu*(theta-tau);
mu2=mu*tau;
Eta(1,1,:,:,:)=lambda2+2*mu2;
Eta(2,2,:,:,:)=Eta(1,1,:,:,:);
Eta(3,3,:,:,:)=Eta(1,1,:,:,:);
Eta(1,2,:,:,:)=lambda2;
Eta(1,3,:,:,:)=Eta(1,2,:,:,:);
Eta(2,3,:,:,:)=Eta(1,2,:,:,:);
Eta(4,4,:,:,:)=mu2;
Eta(5,5,:,:,:)=Eta(4,4,:,:,:);
Eta(6,6,:,:,:)=Eta(4,4,:,:,:);
%}
%% density
rho=ones(nx,ny,nz);
%% Source and source signals
M=2.7;
sx=25;
sy=25;
sz=200;
sn=length(sx);
freq=5;
singles=rickerWave(freq,dt,ns,M);
srcx=zeros(nt,1);
srcy=srcx;
srcz=srcx;
srcx=1*singles;
srcy=1*singles;
srcz=1*singles;
%%

rx=[22,30,40];
ry=[20,30,40];
rz=[30,40,40];
rt=[2,3];

rx=[22:100:1000,1000];
ry=ones(size(rx))*90;
rz=ones(size(rx));
rt=[100,1000,2000,3000,4000,5000,6000,7000,8000];
huge_model=1;

tol=10^-5;

[Rx,Ry,Rz,Rux,Ruy,Ruz,ux,uy,uz]=solver2(dt,dx,dy,dz,nt,nx,ny,nz,huge_model,sx,sy,sz,rt,srcx,srcy,srcz,rx,ry,rz,lp,C,Eta,rho,lpn,Rc,tol);
%%
tt=permute(Rux,[2,1,3,4]);
lim2=.01*[min(ux(:)),max(ux(:))];

col=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0;
    1,1,1;
    0,0,0];
rt2=dt:dt:dt*nt;
t3=zeros(length(rx),nt);
t4=t3;
for i=1:size(t3,1)
    t3(i,:)=real(reshape(ux(rx(i),ry(i),rz(i),:),[1,nt]));
    t4(i,:)=real(reshape(uz(rx(i),ry(i),rz(i),:),[1,nt]));
end

for l2=1:length(rt)
    figure(3)
    %set(gcf,'Visible','on');
    %set(gcf,'position',[0,0,1500,600]);
    
    pcolor3(tt(:,:,:,l2));
    colorbar;
    set(gca,'zdir','reverse');
    xlabel(['x*' num2str(dx) '[m]']);
    ylabel(['y*' num2str(dx) '[m]']);
    zlabel(['z*' num2str(dx) '[m]']);
    title(['t=' num2str(dt*l2) 's']);
    
    for i=1:size(sx,2)
        hold on;
        ax2=plot3(sx(i),sy(i),sz(i),'o','color','red');
    end
    for i=1:size(rx,2)
        hold on;
        ax3=plot3(rx(i),ry(i),rz(i),'.','color',col(i,:),'markersize',20);
    end
    
    legend([ax2,ax3],'source','receiver','Location',[0.5,0.03,0.005,0.005],'orientation','horizontal');
    hold off;
    shg;
    % print(gcf,['.\marmousi2\' num2str(l) '.png'],'-dpng','-r100');
end
%%
figure('name','seismogram');
for i=1:length(rx)
    plot(rt2,Rx(:,i),'color',col(i,:));
    hold on;
end
xlabel('t [s]');
ylabel('ux [m]');
xlim([0,rt2(end)]);
ylim([min(Rx(:)),max(Rx(:))]);
