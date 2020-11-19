function u2=random_media_creation(nx,ny,nz,dx,dy,dz,kappa,a,N)
rng(1)
W=randn(nx,ny,nz);

% sampling interval
ksx=1/dx;
ksy=1/dy;
ksz=1/dz;

nkx=500;
nky=500;
nkz=500;

kx=ksx*((-nkx/2):(nkx/2))/nkx;
ky=ksx*((-nky/2):(nky/2))/nky;
kz=ksx*((-nkz/2):(nkz/2))/nkz;

kx=kx(2:end);
ky=ky(2:end);
kz=kz(2:end);

Lx=nx;
Ly=ny;
Lz=nz;

FW=fftn(W,[nkx,nky,nkz]);
%% spectral filter
FF=R2(kx,ky,kz,kappa,a,N);
Fv=FF.*FW;
%% spacial domain of FF
F=ifftn(FF);
%% inverse transform of velocity in frequency domain
u2=real(ifftn(Fv,[nkx,nky,nkz]));
u2=u2(1:nx,1:ny,1:nz);
end