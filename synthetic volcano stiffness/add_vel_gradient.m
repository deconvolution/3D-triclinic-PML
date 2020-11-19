function v2=add_vel_gradient(nx,ny,nz,dx,dy,dz,v,z1,z2,g)
G3=zeros(1,nz);
G=zeros(nx,ny,nz);
for k=1:length(z1)
    nz2=z2(k)-z1(k)+1;
    G3(z1(k):z2(k))=g*dz*((1:nz2));
    G(:,:,z1(k):z2(k))=repmat(reshape(G3(z1(k):z2(k)),[1,1,z2(k)-z1(k)+1]),[nx,ny]);
end
v2=v+G;
end