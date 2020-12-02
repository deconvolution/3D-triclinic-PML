close all
clear all
load('v.mat');
%% UTM range
% km
a=[min(v.x(:)),max(v.x(:))]
b=[min(v.y(:)),max(v.y(:))]
(a(2)-a(1))/2
(b(2)-b(1))/2
%%
dx=v.x(1,2,1)-v.x(1,1,1);
dy=v.y(2,1,1)-v.y(1,1,1);
dz=v.z(1,1,2)-v.z(1,1,1);
%%
LAT=[-4.96,-.44];
LON=[33.74,38.26];
Z=[min(v.z(:)),max(v.z(:))];
[nx,ny,nz]=size(v.vp);
dlat=(LAT(2)-LAT(1))/nx;
dlon=(LON(2)-LON(1))/ny;
dz=(Z(2)-Z(1))/nz;
[lon,lat,alt]=meshgrid(LON(1)+dlon:dlon:LON(2),LAT(1)+dlat:dlat:LAT(2),Z(1)+dz:dz:Z(2));
%%
v2.lat=lat;
v2.lon=lon;
v2.z=v.z;
v2.easting=v.x;
v2.northing=v.y;
v2.vp=v.vp;
v2.vs=v.vs;
v2.dlat=dlat;
v2.dlon=dlon;
v2.dz=dz;
v2.deasting=dx;
v2.dnorthing=dy;
save('v2.mat');