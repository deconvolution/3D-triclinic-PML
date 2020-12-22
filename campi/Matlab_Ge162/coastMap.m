function coastMap(varargin)
% Basic coastline mapping of the world.
% usage: >coastMap
%
% Options:
%	coastMap('state')    (uses built in mapping toolbox)
%	coastMap('boundary') (uses boundary files in ./Data/)
%	coastMap('plate') (uses boundary files in ./Data/)
%	
%	To use GMT coastline in datalib:
%	coastMap('GMT_c')  crude 
%	coastMap('GMT_l')  low 
%	coastMap('GMT_i')  intermediate
%	coastMap('GMT_h')  high  
%	coastMap('GMT_f')  full 
%
%	coastMap('BoundingBox',bbox,'GMT_i')
%		where bbox = [lon1 lon2 ; lat1 lat2];
%		(faster query if only a subset is used)
%		(w/ GMT datasets, must come before 'GMT_X')
% 
% can set custom viewing range with:
% axis([minlon,maxlon, minlat,maxlat])
%
% % add new points to the map with:
% hh = plot(lons,lats,'b.');
% set(hh,'markersize',16);
%
% % Great circle, see track1 or track2
% % Small circle, see scircle1 or scircle2
% [lat,lon] = track2(37.39,-122.02,49.55,8.67)
% plot(lon,lat,'b');
%
% Also see "worldmap" for built in matlab functions
%   useful for proper projections
%
% Written by Daniel Bowden, dbowden@caltech.edu
%
% Simple coastline data acquired from Paul Earle in 2011
% Geologic boundaries from Fan-Chi Lin 2014
% GMT coastline data acquired from:
%  http://www.soest.hawaii.edu/pwessel/gshhg/
%  protected under GNU license, see data directory for more info


% This will load and plot from coastline.data if no other 
% coastline files are used
use_simple_text = 1;


hold on

%fprintf('Total number of inputs = %d\n',nargin);
nVarargs = length(varargin);
k = 1;
while k <= nVarargs
   fprintf('   %s\n', varargin{k})
   if(strcmp(varargin{k},'BoundingBox'))
	k = k+1;
	bbox = varargin{k}
   elseif(length(strfind(varargin{k},'GMT')))
   	if(strcmp(varargin{k},'GMT_c'))
	   if(exist('bbox'))
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/c/GSHHS_c_L1.shp','BoundingBox',bbox);
	   else
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/c/GSHHS_c_L1.shp');
	   end
	elseif(strcmp(varargin{k},'GMT_l'))
	   if(exist('bbox'))
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/l/GSHHS_l_L1.shp','BoundingBox',bbox);
	   else
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/l/GSHHS_l_L1.shp');
	   end
	elseif(strcmp(varargin{k},'GMT_i'))
	   if(exist('bbox'))
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/i/GSHHS_i_L1.shp','BoundingBox',bbox);
	   else
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/i/GSHHS_i_L1.shp');
	   end
	elseif(strcmp(varargin{k},'GMT_h'))
	   if(exist('bbox'))
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/h/GSHHS_h_L1.shp','BoundingBox',bbox);
	   else
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/h/GSHHS_h_L1.shp');
	   end
	elseif(strcmp(varargin{k},'GMT_f'))
	   if(exist('bbox'))
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/f/GSHHS_f_L1.shp','BoundingBox',bbox);
	   else
		S1 = shaperead('/home/datalib/GMT_Coastlines_Rivers/GMT_GSHHG_2.3.4/GSHHS_shp/f/GSHHS_f_L1.shp');
	   end
	else
		disp(['ERROR! Not recognized: ',varargin{k}])
	end
	use_simple_text = 0;  % Replace use of coastline.data with GMT data
	plot([S1.X],[S1.Y],'k');
	hold on
	set(gca,'dataaspectratio',[1,0.8,1]);

	%axesm('MapProjection','mercator','MapLatLimit',[31 36],'MapLonLimit',[-127 -116])
	%lat = [S1.Y];
	%lon = [S1.X];
	%geoshow(lat,lon);

   elseif(strcmp(varargin{k},'state'))
   	fprintf('    --plotting state boundaries\n')
	% Assumes map toolbox on path, for "usastatelo"
	states = shaperead('usastatelo', 'UseGeoCoords', true);
	for ii = 1:length(states)
		plot(states(ii).Lon, states(ii).Lat,'k','linewidth',1)
	end

   elseif(strcmp(varargin{k},'boundary'))

	fprintf('    --plotting geologic boundaries\n')
	boundary_files = [{'wus_province_II.dat'},{'Grenville.dat'},{'BD.20.dat'},{'rift.dat'}];
	for j = 1:length(boundary_files)
		this_file = ['/home/datalib/Boundaries/',boundary_files{j}];
		fid = fopen(this_file);	
		tline = fgetl(fid);
		clat = [];
		clon = [];
		while(ischar(tline))
			if(length(strfind(tline,'>'))>0)
				clon = [clon, NaN];
				clat = [clat, NaN];
			else
				lonlat = strread(tline,'%f',2);
				clon = [clon, lonlat(1)];
				clat = [clat, lonlat(2)];
			end
			tline = fgetl(fid);
		end
		fclose(fid);
		clon( clon >= 999 ) = NaN;
		clat( clat >= 999 ) = NaN;
		clon(abs(diff(clon))>180) = NaN;

		ishift = find(clon>180);
		clon( ishift ) =  clon( ishift ) - 360;

		hh = plot(clon,clat,'k','linewidth',1);
	end
   elseif(strcmp(varargin{k},'fault'))
	fprintf('    --plotting faults\n')
	this_file = ['/home/datalib/Faults/Holocene_LatestPleistocene.gmt'];
	fid = fopen(this_file);	
	tline = fgetl(fid);
	clat = [];
	clon = [];
	while(ischar(tline))
		if(length(strfind(tline,'>'))>0)
			clon = [clon, NaN];
			clat = [clat, NaN];
		else
			%disp(tline)
			lonlat = strread(tline,'%f',2);
			clon = [clon, lonlat(1)];
			clat = [clat, lonlat(2)];
		end
		tline = fgetl(fid);
	end
	fclose(fid);
	clon( clon >= 999 ) = NaN;
	clat( clat >= 999 ) = NaN;
	clon(abs(diff(clon))>180) = NaN;

	ishift = find(clon>180);
	clon( ishift ) =  clon( ishift ) - 360;

	hh = plot(clon,clat,'r','linewidth',1);
   elseif(strcmp(varargin{k},'plate'))
	fprintf('    --plotting plate boundaries\n')
	this_file = ['/home/datalib/Boundaries/platebound.gmt'];
	fid = fopen(this_file);	
	tline = fgetl(fid);
	clat = [];
	clon = [];
	while(ischar(tline))
		if(length(strfind(tline,'>'))>0)
			clon = [clon, NaN];
			clat = [clat, NaN];
		else
			lonlat = strread(tline,'%f',2);
			clon = [clon, lonlat(1)];
			clat = [clat, lonlat(2)];
		end
		tline = fgetl(fid);
	end
	fclose(fid);
	clon( clon >= 999 ) = NaN;
	clat( clat >= 999 ) = NaN;
	clon(abs(diff(clon))>180) = NaN;

	ishift = find(clon>180);
	clon( ishift ) =  clon( ishift ) - 360;

	hh = plot(clon,clat,'r','linewidth',1);
	%end
   elseif(strcmp(varargin{k},'ocean_white'))
	fprintf('    --ocean_white\n')
	fprintf('      (must have bbox and GMT coasts\n')

	k = k+1;
	in_ocean = varargin{k};

	% expand to cover completely
	bbox2 = bbox;
	bbox2(1,:) = bbox2(1,:) - 0.5;
	bbox2(2,:) = bbox2(2,:) + 0.5;
	
	% find coastline points outside
	inS1 = find(inpolygon([S1.X],[S1.Y],[bbox2(:,1)],[bbox2(:,2)]));
	in_ocean
	patch([S1.X(inS1) in_ocean(1)],[S1.Y(inS1) in_ocean(2)],'w')

	

	
   end
   k = k+1;
end

if(use_simple_text)
	load /home/datalib/Boundaries/coastline.data
	coastline(coastline == 99999) = NaN;
	clat = coastline(:,2);
	clon = coastline(:,1);
	clon(abs(diff(clon))>180) = NaN;
	hh = plot(clon,clat,'k','linewidth',1);
	hold on
	axcl = axis;
	set(gca,'dataaspectratio',[1,0.8,1]);
end
axis([-180 180 -80 80])
if(exist('bbox'))
	axis([ bbox(1,1) bbox(1,2) bbox(2,1) bbox(2,2)])
end

