function [] = UpdatePVD_File(PVDFileName, filename, time)
% This creates a PVD file (if it does not exist yet), and adds a filename
%


% figure out if the file exist
fid = fopen([PVDFileName,'.pvd'],'r');

if fid>0
    
    % determine the # of lines in the file
    numlines=0;
    while ~feof(fid)
        fgetl(fid);
        numlines = numlines+1;
    end
    fclose(fid);
    
    % move to the end of the file
    fid   = fopen([PVDFileName,'.pvd'],'r');
    fid_w = fopen(['temp.pvd'],'w');
    for i=1:numlines-2
        line = fgetl(fid);
        fprintf(fid_w, [line,'\n']);
    end
    fclose(fid);
    
else
    fid_w = fopen(['temp.pvd'],'w');
    
    
    %   create a new file & write header
    fprintf(fid_w,'<?xml version="1.0"?> \n');
    fprintf(fid_w,'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian"> \n');
    fprintf(fid_w,'<Collection>\n');
    
end


% Write new timestep
line = ['	<DataSet timestep="', num2str(time,'%1.6e'), '" file="',filename,'"/> \n']; 
fprintf(fid_w, line);


% write bottom
fprintf(fid_w,'</Collection>\n');
fprintf(fid_w,'</VTKFile>\n');


% rename file
fclose(fid_w);

% rename the file
movefile('temp.pvd',[PVDFileName,'.pvd'])





