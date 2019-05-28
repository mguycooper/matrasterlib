function [ h,ax,c ] = rastercontour( Z,R,varargin )
%RASTERCONTOUR rastercontour(Z,R,varargin) project and display spatially
%referenced raster Z associated with map/geographic raster reference object
%or valid referencing vector/matrix as a 2-d labeled contour plot using
%in-built Matlab functions mapshow.m or geoshow.m. Default behavior sets
%properties ('DisplayType','contour'), ('ShowText','on'), ('Box','on'),
%('TickDir','out'), ('LineWidth',1.5), axis image, and adds a colorbar. The
%function accepts any input to mapshow.m or geoshow.m and returns the
%graphics, axis, and colorbar handle objects so the user can make ex-post
%adjustments as desired.

%   This function is a wrapper for the in-built Matlab functions geoshow.m
%   and mapshow.m (Copyright 2016 The MathWorks, Inc.). The function is
%   designed to display spatially referenced 2-d gridded data as a labeled
%   contour plot by default, whereas geoshow and mapshow are more general,
%   and are designed to display vector or raster data, requiring additional
%   pre- or post-processing to optimize for either. The function determines
%   on-the-fly whether planar (projected) or geographic raster data is
%   passed to it, and inherits the functionality of either geoshow.m or
%   mapshow.m. The function makes reasonable default assumptions about user
%   preferences including setting 'DisplayType' to 'contour', 'ShowText' to
%   'on', adding a colorbar, and a few other graphics settings. All of
%   these can be changed, and the function accepts any argument to mapshow
%   or geoshow, e.g. the user can modify the property 'DisplayType' from
%   the default 'surface' to alternate values such as 'mesh', 'texturemap',
%   or 'contour', or the property 'FaceColor' (shading) to alternate values
%   such as 'flat', 'faceted', or 'interp'. The functions rastersurf.m,
%   rastermesh.m, rastercontour.m, and rastertexture.m are alias' for
%   rastershow.m with alternate displaytypes. Also see rastersurf.m
%   rastersurf3.m, and rastercontour3.m.

%   Author: Matt Cooper, guycooper@ucla.edu, May 2019
%   Citation: Matthew Cooper (2019). matrasterlib
%   (https://www.github.com/mguycooper/matrasterlib), GitHub. Retrieved MMM
%   DD, YYYY

%   Syntax

%   RASTERCONTOUR(...) displays a regular data grid Z and appends a
%   colorbar to the current axes in the default (right) location. By
%   default, DISPLAYTYPE is set to 'contour' and the axis is set image. Z
%   is a 2-D array of class double. If Z is projected in map (planar)
%   coordinates, R must be a referencing matrix or a map raster reference
%   object that relates the subscripts of Z to map coordinates. If Z is
%   geographic data grid, R can be a geographic raster reference object, a
%   referencing vector, or a referencing matrix.
%
%   If R is a geographic raster reference object, its RasterSize property
%   must be consistent with size(Z).
%
%   If R is a referencing vector, it must be a 1-by-3 with elements:
%
%     [cells/degree northern_latitude_limit western_longitude_limit]
%
%   If R is a referencing matrix, it must be 3-by-2 and transform raster
%   row and column indices to/from geographic coordinates according to:
% 
%                     [lon lat] = [row col 1] * R.
%
%   If R is a referencing matrix, it must define a (non-rotational,
%   non-skewed) relationship in which each column of the data grid falls
%   along a meridian and each row falls along a parallel. 
%
%   H = RASTERCONTOUR(...) returns a handle to a MATLAB graphics object
%
%   [H,AX] = RASTERCONTOUR(...) returns a handle to a MATLAB graphics
%   object and a MATLAB axes object
%
%   [H,AX,C] = RASTERCONTOUR(...) returns a handle to a MATLAB graphics
%   object, a MATLAB axes object, and a MATLAB colorbar object

%   RASTERCONTOUR(..., Name, Value) specifies name-value pairs that set
%   MATLAB graphics properties. Parameter names can be abbreviated and are
%   case-insensitive. Refer to the MATLAB Graphics documentation on surface
%   for a complete description of these properties and their values.

%   See also geoshow, mapshow, rastersurf, rastersurf3, rastercontour3

%% Check inputs

% confirm mapping toolbox is installed
assert( license('test','map_toolbox')==1, ...
                        'rastercontour requires Matlab''s Mapping Toolbox.')

% If varargin{2} is not a MapCells or GeographicCellsReference, first
% try converting using refvecToGeoRasterReference, which will fail if
% R is not a vector. Then try refmatToGeoRasterReference, which will
% generally fail because the values in R will produce lat/lon values that
% are inconsistent. Finally try refmatToMapRasterReference

if (~isa(R,'map.rasterref.MapCellsReference')) && ...
        (~isa(R,'map.rasterref.GeographicCellsReference'))
    try 
        R = refvecToGeoRasterReference(R,size(Z));
    catch
        try
            R = refmatToGeoRasterReference(R,size(Z));
        catch
            try 
                R = refmatToMapRasterReference(R,size(Z));
            catch
                error(['Expected input argument 2, R, to be a referencing ' ...
                    'matrix or a map raster reference object that relates ' ...
                    'the subscripts of Z to map coordinates. If Z is ' ...
                    'geographic data grid, R can be a geographic raster ' ...
                    'reference object, a referencing vector, or a ' ...
                    'referencing matrix.']);
            end
        end
    end
end     

% confirm R is either a MapCells or GeographicCellsReference objects. Note,
% this is redundant with try catch block above, but keeping just in case
validateattributes(R, ...
                        {'map.rasterref.MapCellsReference', ...
                        'map.rasterref.GeographicCellsReference'}, ...
                        {'scalar'}, 'rastercontour', 'R', 2)
                    
% confirm Z is a numeric or logical grid of size R.RasterSize
validateattributes(Z,   {'numeric', 'logical'}, ...
                        {'size', R.RasterSize}, 'rastercontour', 'Z', 1)
                    
%% apply the function 

if strcmp(R.CoordinateSystemType,'planar')
    h               =   mapshow(Z,R,'DisplayType','contour', ...
                            'ShowText','on', varargin{:});
elseif strcmp(R.CoordinateSystemType,'geographic')
    h               =   geoshow(Z,R,'DisplayType','contour', ...
                            'ShowText','on', varargin{:});
end


shading flat
axis image
hold on

ax                  =   gca; 
ax.Box              =   'on';
ax.TickDir          =   'out';
ax.LineWidth        =   1.5;
% ax.DataAspectRatio  =   [diff(ax.XLim) diff(ax.XLim) diff(ax.ZLim)];

% add colorbar using default location
c                   =   colorbar;
% uncomment for narrow colorbar
% axpos               =   ax.Position;
% cpos                =   c.Position;
% cpos(3)             =   0.5*cpos(3); % make colorbar half width
% c.Position          =   cpos;
% ax.Position         =   axpos;

%% arrange the outputs

switch nargout
    case 0
    case 1
        varargout{1} = h; 
    
    case 2
        varargout{1} = h; 
        varargout{2} = ax; 
        
    case 3         
        varargout{1} = h; 
        varargout{2} = ax; 
        varargout{3} = c; 
        
    otherwise
        error('Unrecognized number of outputs.') 
end

%% clean up

if nargout==0
    clear h ax c
end

end

