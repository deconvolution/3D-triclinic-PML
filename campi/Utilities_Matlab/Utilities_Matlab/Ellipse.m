dkm=111;
tc=80; %lapse time
v=4; %S-wave velocity
a=v*tc/2/dkm;%principal horizontal axis in degrees
x1=26.5; y1=45.5; z1=50/dkm;% coordinates of the source in km
nf=1; % 0 (near) or 1 (far) station

if nf==0
    x2= 26.6; y2=45.8; z2=0;% coordinates of the near station in km
    hy=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);%hypocentral distance for near station in km
elseif nf==1
    x2= 25; y2=44.5; z2=0; % coordinates of the far station in km
    hy=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);%hypocentral distance for far station in km
end
b=sqrt(a^2-(hy/2)^2); %secondary horizontal axis
ecc=sqrt(1-b^2/a^2); %eccentricity of the ellipse

% Plot the ellipse in red + source and station (Sato, 1978). I have picked from the map
% in Fig. 1
a1 = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
b1 = a*sqrt(1-ecc^2);
t = linspace(0,2*pi);
X = a*cos(t);
Y = b*sin(t);
w = atan2(y2-y1,x2-x1);
x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
plot(x,y,'r-')
hold on
plot(x1,y1,'*')
plot(x2,y2,'^')
hold off
%axis equal

%Average depth of volume (in km) of the medium from which coda wave generation would occur
avdepth=130; % average depth of the seismicity
avv=sqrt(b^2+avdepth/dkm)*dkm; %average depth of the volume