function [T,T2]=rec_conversion(s1,s2,s3,r1,r2,r3,R1,R2,R3,dt,dx,dy,dz,nt,nx,ny,nz,lp,freq)
% T is recording file
% T2 is simulation detail
%%
T=table([r1';r1';r1'],[r2';r2';r2'],[r3';r3';r3'],[ones(length(r3),1);ones(length(r3),1)*2;ones(length(r3),1)*3],[R1;...
    R2;R3]);
T.Properties.VariableNames={'x','y','z','v','time'};
T1=table(dt,nt,0,0,dt:dt:(dt*nt));
T1.Properties.VariableNames={'x','y','z','v','time'};
T=[T1;T];
%%
T2=table([nt;zeros(length(s3)-1,1)],[nx;zeros(length(s3)-1,1)],[ny;zeros(length(s3)-1,1)],...
    [nz;zeros(length(s3)-1,1)],[lp;zeros(length(s3)-1,1)],[dt;zeros(length(s3)-1,1)],[dx;zeros(length(s3)-1,1)],...
    [dy;zeros(length(s3)-1,1)],[dz;zeros(length(s3)-1,1)],freq',s1',s2',s3');
T2.Properties.VariableNames={'nt','nx','ny','nz','lp','dt','dx','dy','dz','freq','s1','s2','s3'};
end
