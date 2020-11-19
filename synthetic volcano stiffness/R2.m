function FF=R2(kx,ky,kz,kappa,a,N)
[kx2,ky2,kz2]=meshgrid(kx,ky,kz);
FF=zeros(size(kx2));
k=sqrt(kx2.^2+ky2.^2+kz2.^2);

%%
%{
% von karman
d=3;
FF=kappa.^2*(a.^-2+k.^2).^(-d/4-N/2);
%}

% gaussian
FF=kappa.*exp(-a.^2.*k.^2/8);
end