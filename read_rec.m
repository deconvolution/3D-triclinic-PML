% load('820604.3878,9686534.895,-9380, rec_S_f_04.mat');
close all;
%%
x=8.2*10^5;
y=9.68*10^6;
rec=rec_S_f_04;
%%
nr=2;
nr2=3;
nr3=4;
u1=cumtrapz((rec.time(nr,:)))*(dt);
u2=cumtrapz((rec.time(nr2,:)))*(dt);
u3=cumtrapz((rec.time(nr3,:)))*(dt);
t=table2array(rec(1,4:end));
%% plot
figure
subplot(3,1,1)
plot(t(1:end-1),u1);
xlabel('time [s]');
ylabel('u1 [m/s]');
title({['x=' num2str(rec.x(nr)) 'm'],['y=' num2str(rec.y(nr)) 'm'],['z=' num2str(rec.z(nr)) 'm']});
subplot(3,1,2)
plot(t(1:end-1),u2);
xlabel('time [s]');
ylabel('u2 [m/s]');
subplot(3,1,3)
plot(t(1:end-1),u3);
xlabel('time [s]');
ylabel('u3 [m/s]');
