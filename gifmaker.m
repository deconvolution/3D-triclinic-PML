function gifmaker(filename,delaytime,sources)
%% read sources
[source_path,~]=fileparts(sources);
source_info=struct2table(dir(fullfile(source_path, '*.png')));

idx = length(source_info.name);
filename = strcat([source_path '/'],filename,'.gif');
%% sort
tt=regexp(source_info.name,'\d*','match');
tt2=str2double(cat(1,tt{:}));
[~,ii]=sort(tt2);
source_info.name=source_info.name(ii);
%%
for i = 1:idx
    myimg = imread(strcat([source_path '/'],source_info.name{i}));
    [A,map] = rgb2ind(myimg,256);
    
    if i==1
        imwrite(A,map,filename,'LoopCount',Inf,'DelayTime',delaytime);
    else
        imwrite(A,map,filename,'WriteMode','append','DelayTime',delaytime);
    end
    fprintf('\n write to gif time step=%d/%d',i,idx);
    d=clock;
    fprintf('\n    current time=%d %d %d %d %d %.0d',d(1),d(2),d(3),d(4),d(5),d(6));
end