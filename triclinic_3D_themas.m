function [R1,R2,R3,v1,v2,v3,E]=triclinic_3D(dt,dx,dy,dz,nt,nx,ny,nz,C,...
    s1,s2,s3,src1,src2,src3,source_type, ...
    r1,r2,r3, ...
    lp,nPML,Rc, ...
    X,Y,Z,...
    x2,y2,z2,plot_interval,view_angle,...
    path)
%% initialize folder for saving
tic;
PVDFileName =[path 'timeinfo'];
if  exist([PVDFileName,'.pvd'],'file')
    delete([PVDFileName,'.pvd'])
end
if ~exist(path,'dir')
    mkdir(path);
end

if ~exist([path 'vts/'],'dir')
    mkdir([path 'vts/']);
end

if ~exist([path 'pic/'],'dir')
    mkdir([path 'pic/']);
end

n_picture=1;
%% vtk write source and receiver locations
if length(s1)==1
    vtkwrite([path 'source.vtk'],'polydata','lines',min(X(:))+[s1,s1]'*dx/1000,min(Y(:))+[s2,s2]'*dy/1000,max(Z(:))-[s3,s3]'*dz/1000);
else
    % change for volcano
    vtkwrite([path 'source.vtk'],'polydata','lines',min(X(:))+s1'*dx/1000,min(Y(:))+s2'*dy/1000,max(Z(:))-s3'*dz/1000);
end

if length(r3)==1
    vtkwrite([path 'receiver.vtk'],'polydata','lines',min(X(:))+[r1,r1]'*dx/1000,min(Y(:))+[r2,r2]'*dy/1000,max(Z(:))-[r3,r3]'*dz/1000);
else
    % change for volcano
    vtkwrite([path 'receiver.vtk'],'polydata','lines',min(Y(:))+r2'*dy/1000,min(X(:))+r1'*dx/1000,max(Z(:))-r3'*dz/1000);
end
%% initialize parameters
format shortg;
% energy
E=zeros(nt,1);

% receiver index
INDr=sub2ind([nx,ny,nz],r1,r2,r3)+2*nx*ny*nz;

% receiver
R1=zeros(length(r3),nt);
R2=R1;
R3=R1;
% velocity
v1=zeros(nx,ny,nz,3);
v2=v1;
v3=v1;
p=v1;

% sigmas
sigmas11=zeros(nx,ny,nz);
sigmas12=sigmas11;
sigmas13=sigmas11;
sigmas22=sigmas11;
sigmas23=sigmas11;
sigmas33=sigmas11;
p=sigmas11;
%% PML coefficient
tt=zeros(nx,ny,nz,3);
tt(:,:,:,1)=C.C11;
tt(:,:,:,2)=C.C22;
tt(:,:,:,3)=C.C33;

% maximum velocity
vmax=sqrt(reshape(max(tt,[],4),[nx,ny,nz])./C.rho);

% PML coefficient
beta1=zeros(nx,ny,nz);
tt=beta1;
tt(2:lp+1,:,:)=repmat(reshape((abs((1:lp)-lp-1)/lp).^nPML,[lp,1,1]),[1,ny,nz]);
tt(nx-lp:nx-1,:,:)=repmat(reshape(((abs(nx-lp-(nx-lp+1:nx)))/lp).^nPML,[lp,1,1]),[1,ny,nz]);
beta1=vmax*(nPML+1)*log(1/Rc)/2/lp/dx.*tt;

beta2=zeros(nx,ny,nz);
tt=beta2;
tt(:,2:lp+1,:)=repmat(reshape((abs((1:lp)-lp-1)/lp).^nPML,[1,lp,1]),[nx,1,nz]);
tt(:,ny-lp:ny-1,:)=repmat(reshape(((abs(nx-lp-(nx-lp+1:nx)))/lp).^nPML,[1,lp,1]),[nx,1,nz]);
beta2=vmax*(nPML+1)*log(1/Rc)/2/lp/dy.*tt;

beta3=zeros(nx,ny,nz);
tt=beta3;
tt(:,:,2:lp+1)=repmat(reshape((abs((1:lp)-lp-1)/lp).^nPML,[1,1,lp]),[nx,ny,1]);
tt(:,:,nz-lp:nz-1)=repmat(reshape(((abs(nz-lp-(nz-lp+1:nz)))/lp).^nPML,[1,1,lp]),[nx,ny,1]);
beta3=vmax*(nPML+1)*log(1/Rc)/2/lp/dz.*tt;

beta1(1,:,:)=beta1(2,:,:);
beta1(end,:,:)=beta1(end-1,:,:);

beta2(:,1,:)=beta2(:,2,:);
beta2(:,end,:)=beta2(:,end-1,:);

beta3(:,:,1)=beta3(:,:,2);
beta3(:,:,end)=beta3(:,:,end-1);

% 3D PML coefficient
IND=unique(find(beta1.*beta2.*beta3~=0));
IND2=unique(find(beta1.*beta2+beta2.*beta3+beta3.*beta1));
IND3=setdiff(IND2,IND);
beta=beta1+beta2+beta3;
beta(IND)=beta(IND)/3;
beta(IND3)=beta(IND3)/2;
clear vmax beta01 beta02 beta03 tt beta1 beta2 beta3 IND IND2 IND3
%% auxiliary variable
v1_t=sigmas11;
v2_t=sigmas11;
v3_t=sigmas11;
v1_x=sigmas11;
v2_x=sigmas11;
v3_x=sigmas11;
v1_y=sigmas11;
v2_y=sigmas11;
v3_y=sigmas11;
v1_z=sigmas11;
v2_z=sigmas11;
v3_z=sigmas11;

% plot limit
lim2=zeros(1,2);
%% for l=2 source allocation
l=1;
for is=1:length(s3)
    switch source_type(is)
        case 'D'
            
            i=s1(is);
            j=s2(is);
            k=s3(is);
            v1(i,j,k,3)=v1(i,j,k,3)+dt./C.rho(i,j,k)*src1(l,is);
            v2(i,j,k,3)=v2(i,j,k,3)+dt./C.rho(i,j,k)*src2(l,is);
            v3(i,j,k,3)=v3(i,j,k,3)+dt./C.rho(i,j,k)*src3(l,is);
            
        case 'P'
            
            i=s1(is);
            j=s2(is);
            k=s3(is);
            
            S1=zeros(3,3,3);
            S2=S1;
            S3=S1;
            
            S1(2,2,2)=src1(l,is);
            S2(2,2,2)=src2(l,is);
            S3(2,2,2)=src3(l,is);
            
            v1(i-1:i+1,j-1:j+1,k-1:k+1,3)=v1(i-1:i+1,j-1:j+1,k-1:k+1,3) ...
                +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S1,1,1),1,1)/dx;
            v2(i-1:i+1,j-1:j+1,k-1:k+1,3)=v2(i-1:i+1,j-1:j+1,k-1:k+1,3)...
                +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S2,1,2),1,2)/dy;
            v3(i-1:i+1,j-1:j+1,k-1:k+1,3)=v3(i-1:i+1,j-1:j+1,k-1:k+1,3)...
                +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S3,1,3),1,3)/dz;
            
            
        case 'S'
            
            i=s1(is);
            j=s2(is);
            k=s3(is);
            S1=zeros(3,3,3);
            S2=S1;
            S3=S1;
            
            S1(2,2,2)=src1(l,is);
            S2(2,2,2)=src2(l,is);
            S3(2,2,2)=src3(l,is);
            
            v1(i-1:i+1,j-1:j+1,k-1:k+1,3)=v1(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S3,1,2),1,2)/dy-diff(ad(S2,1,3),1,3)/dz);
            v2(i-1:i+1,j-1:j+1,k-1:k+1,3)=v2(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S1,1,3),1,3)/dz-diff(ad(S3,1,1),1,1)/dx);
            v3(i-1:i+1,j-1:j+1,k-1:k+1,3)=v3(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S2,1,1),1,1)/dx-diff(ad(S1,1,2),1,2)/dy);
    end
end
%% for l>=3
for l=2:nt-1
    %% shift v in time
    for l2=1:2
        v1(:,:,:,l2)=v1(:,:,:,l2+1);
        v2(:,:,:,l2)=v2(:,:,:,l2+1);
        v3(:,:,:,l2)=v3(:,:,:,l2+1);
    end
    %% compute sigmas
    l4=2;
    v1_x=1/dx*D(v1(:,:,:,l4),1);
    v2_x=1/dx*D(v2(:,:,:,l4),-1);
    v3_x=1/dx*D(v3(:,:,:,l4),-1);
    v1_y=1/dy*D(v1(:,:,:,l4),-2);
    v2_y=1/dy*D(v2(:,:,:,l4),2);
    v3_y=1/dy*D(v3(:,:,:,l4),-2);
    v1_z=1/dz*D(v1(:,:,:,l4),-3);
    v2_z=1/dz*D(v2(:,:,:,l4),-3);
    v3_z=1/dz*D(v3(:,:,:,l4),3);
    
    % add viscous terms
    v1_t=1/dt*(v1(:,:,:,2)-v1(:,:,:,1));
    v2_t=1/dt*(v2(:,:,:,2)-v2(:,:,:,1));
    v3_t=1/dt*(v3(:,:,:,2)-v3(:,:,:,1));
    
    sigmas11=1/3*dt*(4*C.mu2/dx.*D(v1_t,1) ...
        -2*C.mu2/dy.*D(v2_t,2) ...
        -2*C.mu2/dz.*D(v3_t,3) ...
        +(2*C.C11-C.C12-C.C13).*v1_x ...
        +(2*C.C16-C.C26-C.C36).*v1_y ...
        +(2*C.C15-C.C25-C.C35).*v1_z ...
        +(2*C.C16-C.C26-C.C36).*v2_x ...
        +(2*C.C12-C.C22-C.C23).*v2_y ...
        +(2*C.C14-C.C24-C.C34).*v2_z ...
        +(2*C.C15-C.C25-C.C35).*v3_x ...
        +(2*C.C14-C.C24-C.C34).*v3_y ...
        +(2*C.C13-C.C23-C.C33).*v3_z) ...
        +sigmas11...
        -dt*beta.*sigmas11;
    
    sigmas22=1/3*dt*(-2*C.mu2/dx.*D(v1_t,1)...
        +4*C.mu2/dy.*D(v2_t,2) ...
        -2*C.mu2/dz.*D(v3_t,3) ...
        +(-C.C11+2*C.C12-C.C13).*v1_x ...
        +(-C.C16+2*C.C26-C.C36).*v1_y ...
        +(-C.C15+2*C.C25-C.C35).*v1_z ...
        +(-C.C16+2*C.C26-C.C36).*v2_x ...
        +(-C.C12+2*C.C22-C.C23).*v2_y ...
        +(-C.C14+2*C.C24-C.C34).*v2_z ...
        +(-C.C15+2*C.C25-C.C35).*v3_x ...
        +(-C.C14+2*C.C24-C.C34).*v3_y ...
        +(-C.C13+2*C.C23-C.C33).*v3_z) ...
        +sigmas22...
        -dt*beta.*sigmas22;
    
    sigmas33=1/3*dt*(-2*C.mu2/dx.*D(v1_t,1)...
        -2*C.mu2/dy.*D(v2_t,2) ...
        +4*C.mu2/dz.*D(v3_t,3) ...
        +(-1*C.C11-C.C12+2*C.C13).*v1_x ...
        +(-1*C.C16-C.C26+2*C.C36).*v1_y ...
        +(-1*C.C15-C.C25+2*C.C35).*v1_z ...
        +(-1*C.C16-C.C26+2*C.C36).*v2_x ...
        +(-1*C.C12-C.C22+2*C.C23).*v2_y ...
        +(-1*C.C14-C.C24+2*C.C34).*v2_z ...
        +(-1*C.C15-C.C25+2*C.C35).*v3_x ...
        +(-1*C.C14-C.C24+2*C.C34).*v3_y ...
        +(-1*C.C13-C.C23+2*C.C33).*v3_z) ...
        +sigmas33...
        -dt*beta.*sigmas33;
    
    sigmas12=dt*(C.mu2/dx.*D(v2_t,-1)+C.mu2/dy.*D(v1_t,-2) ...
        +C.C16.*(v1_x)+C.C56.*(v3_x)+C.C66.*(v2_x)+C.C26.*(v2_y)+C.C46.*(v3_y)+C.C66.*(v1_y)+C.C36.*(v3_z)+C.C46.*(v2_z)+C.C56.*(v1_z))...
        +sigmas12...
        -dt*beta.*sigmas12;
    
    sigmas13=dt*(C.mu2/dx.*D(v3_t,-1)+C.mu2/dz.*D(v1_t,-3) ...
        +C.C15.*(v1_x)+C.C56.*(v2_x)+C.C55.*(v3_x)+C.C25.*(v2_y)+C.C45.*(v3_y)+C.C56.*(v1_y)+C.C35.*(v3_z)+C.C45.*(v2_z)+C.C55.*(v1_z)) ...
        +sigmas13...
        -dt*beta.*sigmas13;
    
    sigmas23=dt*(C.mu2/dy.*D(v3_t,-2)+C.mu2/dz.*D(v1_t,-3) ...
        +C.C14.*(v1_x)+C.C46.*(v2_x)+C.C45.*(v3_x)+C.C24.*(v2_y)+C.C46.*(v1_y)+C.C44.*(v3_y)+C.C34.*(v3_z)+C.C45.*(v1_z)+C.C44.*(v2_z)) ...
        +sigmas23...
        -dt*beta.*sigmas23;
    %% compute p
    p=-1/3*dt*((2*C.lambda2+2*C.mu2)/dx.*D(v1_t,1)...
        +(2*C.lambda2+2*C.mu2)/dy.*D(v2_t,2) ...
        +(2*C.lambda2+2*C.mu2)/dz.*D(v3_t,3) ...
        +(C.C11+C.C12+C.C13).*v1_x ...
        +(C.C16+C.C26+C.C36).*v1_y ...
        +(C.C15+C.C25+C.C35).*v1_z ...
        +(C.C16+C.C26+C.C36).*v2_x ...
        +(C.C12+C.C22+C.C23).*v2_y ...
        +(C.C14+C.C24+C.C34).*v2_z ...
        +(C.C15+C.C25+C.C35).*v3_x ...
        +(C.C14+C.C24+C.C34).*v3_y ...
        +(C.C13+C.C23+C.C33).*v3_z) ...
        +p...
        -dt*beta.*p;
    %% compute v
    v1(:,:,:,3)=dt./C.rho.*(1/dx*D(sigmas11-p,-1)+1/dy*D(sigmas12,2)+1/dz*D(sigmas13,3))+v1(:,:,:,2)...
        -dt*beta.*v1(:,:,:,2);
    
    v2(:,:,:,3)=dt./C.rho.*(1/dx*D(sigmas12,1)+1/dy*D(sigmas22-p,-2)+1/dz*D(sigmas23,3))+v2(:,:,:,2)...
        -dt*beta.*v2(:,:,:,2);
    
    v3(:,:,:,3)=dt./C.rho.*(1/dx*D(sigmas13,1)+1/dy*D(sigmas23,2)+1/dz*D(sigmas33-p,-3))+v3(:,:,:,2)...
        -dt*beta.*v3(:,:,:,3);
    
    % allocate source
    for is=1:length(s3)
        switch source_type(is)
            case 'D'
                
                i=s1(is);
                j=s2(is);
                k=s3(is);
                v1(i,j,k,3)=v1(i,j,k,3)+dt./C.rho(i,j,k)*src1(l,is);
                v2(i,j,k,3)=v2(i,j,k,3)+dt./C.rho(i,j,k)*src2(l,is);
                v3(i,j,k,3)=v3(i,j,k,3)+dt./C.rho(i,j,k)*src3(l,is);
                
            case 'P'
                
                i=s1(is);
                j=s2(is);
                k=s3(is);
                
                S1=zeros(3,3,3);
                S2=S1;
                S3=S1;
                
                S1(2,2,2)=src1(l,is);
                S2(2,2,2)=src2(l,is);
                S3(2,2,2)=src3(l,is);
                
                v1(i-1:i+1,j-1:j+1,k-1:k+1,3)=v1(i-1:i+1,j-1:j+1,k-1:k+1,3) ...
                    +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S1,1,1),1,1)/dx;
                v2(i-1:i+1,j-1:j+1,k-1:k+1,3)=v2(i-1:i+1,j-1:j+1,k-1:k+1,3)...
                    +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S2,1,2),1,2)/dy;
                v3(i-1:i+1,j-1:j+1,k-1:k+1,3)=v3(i-1:i+1,j-1:j+1,k-1:k+1,3)...
                    +1/3*dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*diff(ad(S3,1,3),1,3)/dz;
                
                
            case 'S'
                
                i=s1(is);
                j=s2(is);
                k=s3(is);
                S1=zeros(3,3,3);
                S2=S1;
                S3=S1;
                
                S1(2,2,2)=src1(l,is);
                S2(2,2,2)=src2(l,is);
                S3(2,2,2)=src3(l,is);
                
                v1(i-1:i+1,j-1:j+1,k-1:k+1,3)=v1(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S3,1,2),1,2)/dy-diff(ad(S2,1,3),1,3)/dz);
                v2(i-1:i+1,j-1:j+1,k-1:k+1,3)=v2(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S1,1,3),1,3)/dz-diff(ad(S3,1,1),1,1)/dx);
                v3(i-1:i+1,j-1:j+1,k-1:k+1,3)=v3(i-1:i+1,j-1:j+1,k-1:k+1,3)+dt./C.rho(i-1:i+1,j-1:j+1,k-1:k+1).*(diff(ad(S2,1,1),1,1)/dx-diff(ad(S1,1,2),1,2)/dy);
        end
    end
    
    % fixed-boundary condition
    v1(1,:,:,3)=0;
    v1(end,:,:,3)=0;
    v1(:,1,:,3)=0;
    v1(:,end,:,3)=0;
    v1(:,:,1,3)=0;
    v1(:,:,end,3)=0;
    
    v2(1,:,:,3)=0;
    v2(end,:,:,3)=0;
    v2(:,1,:,3)=0;
    v2(:,end,:,3)=0;
    v2(:,:,1,3)=0;
    v2(:,:,end,3)=0;
    
    v3(1,:,:,3)=0;
    v3(end,:,:,3)=0;
    v3(:,1,:,3)=0;
    v3(:,end,:,3)=0;
    v3(:,:,1,3)=0;
    v3(:,:,end,3)=0;
    
    % calculate energy
    E(l+1)=.5*sum(C.rho.*dx*dy*dz.*(v1.^2+dy*v2.^2+v3.^2),'all');
    %% assign recordings
    R1(:,l+1)=v1(INDr);
    R2(:,l+1)=v2(INDr);
    R3(:,l+1)=v3(INDr);
    %% plot
    % visualize v3
    if mod(l,plot_interval)==0 || l==nt
        
        figure('visible','off');
        set(gcf,'position',[80,80,1000,800]);
        % need change later
        ax=slice(X,Y,Z,v3(:,:,:,3),y2*dy,x2*dx,z2*dz);
        if length(view_angle)==2
            view(view_angle);
        end
        set(ax,'edgecolor','none');
        tt=caxis;
        lim2(1)=min(lim2(1),tt(1));
        lim2(2)=max(lim2(2),tt(2));
        caxis(.04*lim2);
        colorbar;
        xlabel('y [m]');
        ylabel('x [m]');
        zlabel('z [m]');
        set(gca,'zdir','reverse');
        title({['v3 [m/s]'],['t=' num2str(dt*(l+1)) 's']});
        hold on;
        %{
        for is=1:length(s3)
            hold on;
            ax2=plot3(s2*dy,s1*dx,s3*dz,'v','color','red');
        end
        for ir=1:length(r3)
            hold on;
            ax3=plot3(r2*dy,r1*dx,r3*dz,'^','color','black');
        end
        
        legend([ax2,ax3],...
            'source','receiver',...
            'Location',[0.5,0.02,0.005,0.002],'orientation','horizontal');
        %}
        print(gcf,[path 'pic/' num2str(n_picture) '.png'],'-dpng','-r200')
        %% save vtk
        c11=C.C11;
        filename = [path 'vts/velocity_model_Timestep', num2str(n_picture)];
        
        data.X=X;
        data.Y=Y;
        data.Z=Z;
        data.SCAL=c11;
        data.SCAL_name='C11';
        data.VEC_X=v1(:,:,:,3);
        data.VEC_Y=v2(:,:,:,3);
        data.VEC_Z=v3(:,:,:,3);
        data.VEC_name='velocity';
        mat3D2vts(filename,data);
        UpdatePVD_File(PVDFileName, ['./vts/velocity_model_Timestep' num2str(n_picture) '.vts'], l*dt);
        
        n_picture=n_picture+1;
        
    end
    %% timing
    fprintf('\ntime step=%d/%d',l+1,nt);
    fprintf('\n    epalsed time=%.2fs',toc);
    fprintf('\n    n_picture=%d',n_picture);
    d=clock;
    fprintf('\n    current time=%d %d %d %d %d %.0d  ',d(1),d(2),d(3),d(4),d(5),d(6));
end
%% tell where the files are saved
fprintf(['\npictures, vts, vts files are stored in \t' path]);

end