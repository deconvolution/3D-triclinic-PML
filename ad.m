function u2=ad(u,n,i)
switch i
    case 1
        u2=zeros(size(u,1)+n,size(u,2),size(u,3));
        u2(1:end-n,:,:)=u;
    case 2
        u2=zeros(size(u,1),size(u,2)+n,size(u,3));
        u2(:,1:end-n,:)=u;
    case 3
        u2=zeros(size(u,1),size(u,2),size(u,3)+n);
        u2(:,:,1:end-n)=u;
    case -1
        u2=zeros(size(u,1)+n,size(u,2),size(u,3));
        u2(n+1:end,:,:)=u;
    case -2
        u2=zeros(size(u,1),size(u,2)+n,size(u,3));
        u2(:,n+1:end,:)=u;
    case -3
        u2=zeros(size(u,1),size(u,2),size(u,3)+n);
        u2(:,:,n+1:end)=u;
end