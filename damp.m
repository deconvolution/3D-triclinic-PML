function [C1,C2,C3,C7,C8,C12,C16,C19,C21,lambda2,mu2]=damp(C1,C2,C3,C7,C8,C12,C16,C19,C21,lambda2,mu2,l_damp,Qs)
% l_damp
if length(Qs)==1
    Qs=inf;end
if length(Qs)==1
    Qs=inf;
end
%%
C1=C1;
C7=C7;
C8=C8;
C12=C12;
C16=C16;
C19=C19;
C21=C21;
C2=C2;
C3=C3;
lambda2=lambda2;
mu2=mu2;
[nx,ny,nz]=size(C1);
%% create cube for corners
cube=zeros(l_damp,l_damp,l_damp);
for i=1:l_damp
    cube(i,:,:)=Qs(l_damp-i+1);
    cube(:,i,:)=Qs(l_damp-i+1);
    cube(:,:,i)=Qs(l_damp-i+1);
end
%% create plane
plane=zeros(l_damp,l_damp);
for i=1:l_damp
    plane(i,:)=Qs(l_damp-i+1);
    plane(:,i)=Qs(l_damp-i+1);
end
%% corner
cube2=zeros(nx,ny,nz);
cube2(1:l_damp,1:l_damp,1:l_damp)=flip(flip(flip(cube,1),2),3);
cube2(end-l_damp+1:end,end-l_damp+1:end,end-l_damp+1:end)=cube;
cube2(end-l_damp+1:end,1:l_damp,1:l_damp)=flip(flip(cube,2),3);
cube2(1:l_damp,end-l_damp+1:end,1:l_damp)=flip(flip(cube,1),3);
cube2(1:l_damp,1:l_damp,end-l_damp+1:end)=flip(flip(cube,1),2);
cube2(1:l_damp,end-l_damp+1:end,end-l_damp+1:end)=flip(cube,1);
cube2(end-l_damp+1:end,1:l_damp,end-l_damp+1:end)=flip(cube,2);
cube2(end-l_damp+1:end,end-l_damp+1:end,1:l_damp)=flip(cube,3);
%%
cube2(1:l_damp,l_damp+1:end-l_damp,l_damp+1:end-l_damp)=repmat(reshape((Qs),[l_damp,1,1]),[1,ny-2*l_damp,nz-2*l_damp]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(Qs),[l_damp,1,1]),[1,ny-2*l_damp,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape((Qs),[1,l_damp,1]),[nx-2*l_damp,1,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(flip(Qs),[1,l_damp,1]),[nx-2*l_damp,1,nz-2*l_damp]);
cube2(l_damp+1:end-l_damp,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape((Qs),[1,1,l_damp]),[nx-2*l_damp,ny-2*l_damp,1]);
cube2(l_damp+1:end-l_damp,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(flip(Qs),[1,1,l_damp]),[nx-2*l_damp,ny-2*l_damp,1]);
%
cube2(l_damp+1:end-l_damp,1:l_damp,1:l_damp)=repmat(reshape(flip(flip(plane,1),2),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,end-l_damp+1:end)=repmat(reshape(plane,[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,1:l_damp,end-l_damp+1:end)=repmat(reshape(flip(plane,1),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);
cube2(l_damp+1:end-l_damp,end-l_damp+1:end,1:l_damp)=repmat(reshape(flip(plane,2),[1,l_damp,l_damp]),[nx-2*l_damp,1,1]);

cube2(1:l_damp,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape(flip(flip(plane,1),2),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(plane,[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(1:l_damp,l_damp+1:end-l_damp,end-l_damp+1:end)=repmat(reshape(flip(plane,1),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);
cube2(end-l_damp+1:end,l_damp+1:end-l_damp,1:l_damp)=repmat(reshape(flip(plane,2),[l_damp,1,l_damp]),[1,ny-2*l_damp,1]);

cube2(1:l_damp,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(flip(plane,1),2),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(end-l_damp+1:end,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(plane,[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(1:l_damp,end-l_damp+1:end,l_damp+1:end-l_damp)=repmat(reshape(flip(plane,1),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
cube2(end-l_damp+1:end,1:l_damp,l_damp+1:end-l_damp)=repmat(reshape(flip(plane,2),[l_damp,l_damp,1]),[1,1,nz-2*l_damp]);
%%
mu2(:,:,l_damp+1:end)=C16(:,:,l_damp+1:end).*cube2(:,:,l_damp+1:end);
lambda2(:,:,l_damp+1:end)=C12(:,:,l_damp+1:end).*cube2(:,:,l_damp+1:end);
end