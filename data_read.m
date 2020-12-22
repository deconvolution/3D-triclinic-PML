clear all;
close all;
nmarker0='./151007/*N.sac';

list                        = dir(nmarker0);
filenames                   = {list.name}'; %create cell array of file names
filedir                     = {list.folder}'; %create cell array of file folder
listasac                    = strcat(filedir,'/',filenames);
lls                         = length(listasac);
tt          =   readsac(listasac{1}); % Imports SAC files
%%
fSource=2;
tt=readsac(listasac{1});
fs=1/tt.DELTA;
dt=tt.DELTA;
[b,a] = butter(4,([fSource-fSource/3 fSource+fSource/3]/fs*2),'bandpass');

for i=1:lls
    tt=readsac(listasac{i}); % Imports SAC files
    u=filter(b,a,detrend(tt.DATA1,1));
    
    subplot(3,3,i)
    
    plot(dt:dt:dt*length(tt.DATA1),u,'b','LineWidth',2);
    
    xlabel('t [s]'); ylabel('u [m]')
end
%%
fc = 300;
fs = 1000;



dataIn = randn(1000,1);
dataOut = filter(b,a,dataIn);