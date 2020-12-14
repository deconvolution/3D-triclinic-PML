function [f,nw]=morlet1(fb,dt)
% wb=2*pi*fb;
t0=6/5/fb;
% Dw=0.5*wb;
nw=2*t0/dt;
% n=1:nw;
% t=(n-1)*dt;
% lb = -4;
% ub = 4;
lb = -2;
ub = 2;
[f,~] = morlet(lb,ub,nw);

% plot(xval,f)
% grid on
% title('Morlet Wavelet')
% D=t-t0;
% f=exp(-Dw*Dw*D.*D/4).*cos(wb*D);