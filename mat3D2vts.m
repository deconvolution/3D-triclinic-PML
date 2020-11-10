function mat3D2vts(filename,data)
%
% transform matlab data (scalar or vecotr to pv format)
% Tobias Baumann, Mainz, 2013
%

% Transform the arrays in vector-format, and change them to the appropriate
% precision. This example file assumes that data comes in single precision
% Double precision requires changing in ths code below (Float32 into Float64)


if isfield(data,'ASCII')
    ASCII       = data.ASCII;
else
    ASCII       = false;                       % ASCII or BINARY?
end

if isfield(data,'Lat')
    % we have data on Latitude/Longitude/Depth format
    data.R = 6371 + data.Depth;         % this assumes that depth is given in km, and the depth is negative (-1000 is at 1000 km depth)
    if any(data.R(:))<0
        error('Are you sure that you gave depth in units of kilometer [and depth as negative numbers]?')
    end
    
    [data.X,data.Y,data.Z] = LonLatR2xyz(data.Lon,data.Lat,data.R);
    
end


if(isfield(data,'X') && isfield(data,'Y') && isfield(data,'Z'))
    Points          = [single(data.X(:)),single(data.Y(:)),single(data.Z(:))];
end


% scalar data exists ?
if(isfield(data,'SCAL'))
    scalar          = data.SCAL(:);
else
    scalar          = [];
end

if(isfield(data,'SCAL1'))
    scalar1          = data.SCAL1(:);
else
    scalar1          = [];
end

if(isfield(data,'SCAL2'))
    scalar2          = data.SCAL2(:);
else
    scalar2          = [];
end

if(isfield(data,'SCAL_name'))
    fieldname_scalar   = data.SCAL_name;
else
    fieldname_scalar   = 'scalarfield_data';
end

if(isfield(data,'SCAL1_name'))
    fieldname_scalar1   = data.SCAL1_name;
else
    fieldname_scalar1   = 'scalarfield1_data';
end

if(isfield(data,'SCAL2_name'))
    fieldname_scalar2   = data.SCAL2_name;
else
    fieldname_scalar2   = 'scalarfield2_data';
end

% vector data exists ?
if(isfield(data,'VEC_X') && isfield(data,'VEC_Y') && isfield(data,'VEC_Z'))
    vector          = [single(data.VEC_X(:)),single(data.VEC_Y(:)),single(data.VEC_Z(:))];
else
    vector          = [];
end
if(isfield(data,'VEC_name'))
    fieldname_vector   = data.VEC_name;
else
    fieldname_vector   = 'vectorfield_data';
end




%==========================================================================
% Definitions and initialization
sizeof_Float32  =   4;
sizeof_Float64  =   4;
sizeof_UInt32   =   4;
Offset          =   0;      % Initial offset

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fname_vtk       = [filename '.vts'];
fid             = fopen(fname_vtk,'w','b');           % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'<?xml version="1.0"?> \n');
fprintf(fid,'<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian" >\n');
fprintf(fid,'  <StructuredGrid  WholeExtent="%i %i %i %i %i %i">\n', [0 size(data.X,1)-1 0 size(data.X,2)-1 0 size(data.X,3)-1]);
fprintf(fid,'  <Piece Extent="%i %i %i %i %i %i">\n', [0 size(data.X,1)-1 0 size(data.X,2)-1 0 size(data.X,3)-1]);

%--------------------------------------------------------------------------
% Add point-wise data
%--------------------------------------------------------------------------
% if isempty(fieldname_scalar1)
%     fprintf(fid,'    <PointData Scalars=" %s " Vectors=" %s"  >\n',fieldname_scalar,fieldname_vector);
% else
fprintf(fid,'    <PointData Scalars=" %s %s %s " Vectors=" %s"  >\n',fieldname_scalar,fieldname_scalar1, fieldname_scalar2, fieldname_vector);
% end


% SCALAR ----------------
if ~isempty(scalar)
    if ASCII
        % ASCII:
        if isinteger(data.SCAL)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="ascii">\n',fieldname_scalar);
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="ascii">\n',fieldname_scalar);
        end
        
        for i=1:length(scalar)
            fprintf(fid,'        %g \n',single(scalar(i)));
        end
    else
        % BINARY:
        if isinteger(data.SCAL)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar, int32(Offset));
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar, int32(Offset));
        end
        
        Offset = Offset + length(scalar(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
end
% -----------------------

% SCALAR1 ----------------
if ~isempty(scalar1)
    if ASCII
        % ASCII:
        if isinteger(data.SCAL1)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="ascii">\n',fieldname_scalar1);
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="ascii">\n',fieldname_scalar1);
        end
        
        for i=1:length(scalar1)
            fprintf(fid,'        %g \n',single(scalar1(i)));
        end
    else
        % BINARY:
        if isinteger(data.SCAL1)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar1, int32(Offset));
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar1, int32(Offset));
        end
        
        Offset = Offset + length(scalar1(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
end
% -----------------------

% SCALAR2 ----------------
if ~isempty(scalar2)
    if ASCII
        % ASCII:
        if isinteger(data.SCAL2)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="ascii">\n',fieldname_scalar2);
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="ascii">\n',fieldname_scalar2);
        end
        
        for i=1:length(scalar2)
            fprintf(fid,'        %g \n',single(scalar2(i)));
        end
    else
        % BINARY:
        if isinteger(data.SCAL1)
            fprintf(fid,'      <DataArray type="Int32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar2, int32(Offset));
        else
            fprintf(fid,'      <DataArray type="Float32" Name="%s" format="appended" offset="%i">\n',fieldname_scalar2, int32(Offset));
        end
        
        Offset = Offset + length(scalar2(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
end
% -----------------------

% VECTOR ----------------
if ~isempty(vector)
    if ASCII
        % ASCII:
        fprintf(fid,'      <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="ascii">\n',fieldname_vector);
        for i=1:length(T)
            fprintf(fid,'   %g %g %g \n',single(vector(i,:)));
        end
    else
        % BINARY:
        fprintf(fid,'      <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%i">\n',fieldname_vector,int32(Offset));
        Offset = Offset + length(vector(:))*sizeof_Float32 + 1*sizeof_UInt32;
    end
    fprintf(fid,'      </DataArray>\n');
end
% -----------------------


fprintf(fid,'    </PointData>\n');
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------

fprintf(fid,'    <Celldata>\n');
fprintf(fid,'    </Celldata>\n');


%--------------------------------------------------------------------------
% Add coordinates of structured grid
%--------------------------------------------------------------------------
fprintf(fid,'    <Points>\n');

% ASCII
if ASCII
    if isinteger(data.SCAL)
        fprintf(fid,'      <DataArray type="Int32" Name="Array" NumberOfComponents="3" format="ascii">\n');
        for i=1:size(Points,1)
            fprintf(fid,' %g %g %g \n',[int32(Points(i,:))]);
        end
    else
        fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="ascii">\n');
        for i=1:size(Points,1)
            fprintf(fid,' %g %g %g \n',[Points(i,:)]);
        end
    end
    
    
else
    fprintf(fid,'      <DataArray type="Float32" Name="Array" NumberOfComponents="3" format="appended" offset="%i" >\n',int32(Offset));
    
end
fprintf(fid,'      </DataArray>\n');
fprintf(fid,'    </Points>\n');
%--------------------------------------------------------------------------

fprintf(fid,'  </Piece> \n');
fprintf(fid,'  </StructuredGrid> \n');


if ~ASCII
    % Append binary data in raw format: the order in which data arrays are
    % added should be the same as how they are defined above
    fprintf(fid,'  <AppendedData encoding="raw"> \n');
    fprintf(fid,'_');
    
    if ~isempty(scalar)
        fwrite(fid,uint32(length(scalar)*sizeof_Float32),'uint32');
        if isinteger(scalar)
            % Add scalar data in int32 format
            fwrite(fid,int32(scalar).'      ,   'int32');
        else
            % Add scalar data in binary format (convert to single precision
            % to save space)
            fwrite(fid,single(scalar).'      ,   'float32');
        end
        
    end
    
    if ~isempty(scalar1)
        fwrite(fid,uint32(length(scalar1)*sizeof_Float32),'uint32');
        if isinteger(scalar1)
            % Add scalar data in int32 format
            fwrite(fid,int32(scalar1).'      ,   'int32');
        else
            % Add scalar data in binary format (convert to single precision
            % to save space)
            fwrite(fid,single(scalar1).'      ,   'float32');
        end
        
    end
    if ~isempty(scalar2)
        fwrite(fid,uint32(length(scalar2)*sizeof_Float32),'uint32');
        if isinteger(scalar2)
            % Add scalar data in int32 format
            fwrite(fid,int32(scalar2).'      ,   'int32');
        else
            % Add scalar data in binary format (convert to single precision
            % to save space)
            fwrite(fid,single(scalar2).'      ,   'float32');
        end
        
    end
    
    if ~isempty(vector)
        % Add vector data in binary format
        fwrite(fid,uint32(length(vector(:))*sizeof_Float32),'uint32');
        fwrite(fid,single(vector).' ,   'float32');
    end
    
    % Add Coordinates in binary format
    fwrite(fid,uint32(length(Points(:))*sizeof_Float32),'uint32');
    fwrite(fid,single(Points).' ,   'float32');
    
    fprintf(fid,'  </AppendedData> \n');
end


fprintf(fid,'</VTKFile>\n');
fclose(fid);


if ~ASCII
    disp(['Created Binary XML-VTK output file ',fname_vtk])
else
    disp(['Created ASCII XML-VTK output file ',fname_vtk])
end

end