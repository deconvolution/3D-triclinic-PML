function C2=add_fracture(lambda,G,e,theta)
nu=lambda./2./(lambda+G);
[nx,ny,nz]=size(G);
C=zeros(6,6,nx,ny,nz);
C(1,1,:,:,:)=lambda+2*G;
C(2,2,:,:,:)=lambda+2*G;
C(3,3,:,:,:)=lambda+2*G;
C(1,2,:,:,:)=lambda;
C(1,3,:,:,:)=lambda;
C(2,3,:,:,:)=lambda;
C(4,4,:,:,:)=G;
C(5,5,:,:,:)=G;
C(6,6,:,:,:)=G;
C=fill_lower_triangle(C);

S0=inv_C_S(C);
%%
E=G.*(3*lambda+2.*G)./(lambda+G);
gammab=G./(E.*nu./(1+nu)./(1-2*nu)+E./(1+nu));
Et=16./3./(3-2*gammab).*e;
En=4./3./gammab./(1-gammab).*e;
Zn=En./(E.*nu./(1+nu)./(1-2*nu)+E./2./(1+nu));
Zt=Et./(E./2./(1+nu));
Sf0=zeros(6,6,nx,ny,nz);
Sf=Sf0;
Sf0(5,5,:,:,:)=Zt;
Sf0(6,6,:,:,:)=Zt;
Sf0(1,1,:,:,:)=Zn;
Sf0(isnan(Sf0))=0;
%% Bond transform
M=bond_transform(theta,nx,ny,nz);
for i=1:nx
    for j=1:ny
        for k=1:nz
            Sf(:,:,i,j,k)=M(:,:,i,j,k)*Sf0(:,:,i,j,k)*M(:,:,i,j,k)';
        end
    end
end
%%
S=S0+Sf;
%%
C2=inv_C_S(S);
C2(isnan(C2))=0;
end