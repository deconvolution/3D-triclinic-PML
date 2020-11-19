% load('820604.3878,9686534.895,-9380, rec_S_f_04.mat');
close all;
fprintf('UTM easting range');
disp([min(rec_S_f_04.easting(4:end)),max(rec_S_f_04.easting(4:end))])
fprintf('\n');
fprintf('UTM northing range');
disp([min(rec_S_f_04.northing(4:end)),max(rec_S_f_04.northing(4:end))])
%%
x=8.2*10^5;
y=9.68*10^6;
rec=rec_S_f_04;
%% extract characteristics
nx=rec.easting(1);
ny=rec.northing(1);
nz=rec.elevation(1);
l_damp=rec.easting(2);
t=rec.time(3,:);
%% v
% v1
[~,tt]=min((rec.easting(4:(nx-2*l_damp)*(ny-2*l_damp)+3)-x).^2+(rec.northing(4:(nx-2*l_damp)*(ny-2*l_damp)+3)-y).^2);
tt=tt+3;
% v2
[~,tt2]=min((rec.easting((nx-2*l_damp)*(ny-2*l_damp)+4:(nx-2*l_damp)*(ny-2*l_damp)*2+3)-x).^2+(rec.northing((nx-2*l_damp)*(ny-2*l_damp)+4:(nx-2*l_damp)*(ny-2*l_damp)*2+3)-y).^2);
tt2=tt2+(nx-2*l_damp)*(ny-2*l_damp)+3;
% v3
[~,tt3]=min((rec.easting((nx-2*l_damp)*(ny-2*l_damp)*2+4:(nx-2*l_damp)*(ny-2*l_damp)*3+3)-x).^2+(rec.northing((nx-2*l_damp)*(ny-2*l_damp)*2+4:(nx-2*l_damp)*(ny-2*l_damp)*3+3)-y).^2);
tt3=tt3+2*(nx-2*l_damp)*(ny-2*l_damp)+3;
%%
u1=cumtrapz((rec.time(tt,:)))*(t(2)-t(1));
u2=cumtrapz((rec.time(tt2,:)))*(t(2)-t(1));
u3=cumtrapz((rec.time(tt3,:)))*(t(2)-t(1));
%% plot
figure
subplot(3,1,1)
plot(t,u1);
xlabel('time [s]');
ylabel('u1 [m/s]');
title({['x=' num2str(x) 'm'],['y=' num2str(y) 'm'],['z=' num2str(rec.elevation(tt)) 'm']});
subplot(3,1,2)
plot(t,u2);
xlabel('time [s]');
ylabel('u2 [m/s]');
subplot(3,1,3)
plot(t,u3);
xlabel('time [s]');
ylabel('u3 [m/s]');
