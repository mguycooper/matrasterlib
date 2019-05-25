function [ h,ax,c ] = geoshowsurface( Z,R,varargin )
%GEOSHOWSURFACE alias for GEOSHOW with default "DisplayType" set to
%"surface", adds a colorbar, sets axis tight, box on, tickdir out, and
%linewidth = 2

%   This function is an alias for the in-built Matlab function
%   geoshow.m (Copyright 2016 The MathWorks, Inc.), designed to map 
%   2-d gridded data by default with a 'Surface' display type and a
%   colorbar, and some improved graphics the user can adjust ex-post

%   Author: Matt Cooper, guycooper@ucla.edu, May 2019

%   Syntax

%   GEOSHOWSURFACE(...) displays a regular data grid Z and appends a
%   colorbar to the current axes in the default (right) location. By
%   default, DISPLAYTYPE is set to 'surface' and the axis is set tight. Z
%   is a 2-D array of class double. R can be a geographic raster reference
%   object, a referencing vector, or a referencing matrix.
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
%   H = GEOSHOWSURFACE(...) returns a handle to a MATLAB graphics object
%
%   [H,AX] = GEOSHOWSURFACE(...) returns a handle to a MATLAB graphics
%   object and a MATLAB axes object
%
%   [H,AX,C] = GEOSHOWSURFACE(...) returns a handle to a MATLAB graphics
%   object, a MATLAB axes object, and a MATLAB colorbar object

%   GEOSHOWSURFACE(..., Name, Value) specifies name-value pairs that set
%   MATLAB graphics properties. Parameter names can be abbreviated and are
%   case-insensitive. Refer to the MATLAB Graphics documentation on surface
%   for a complete description of these properties and their values.

%   See also geoshow, mapshow

h = geoshow(Z,R,'DisplayType','surface',varargin{:}); ax = gca;
axis tight
c = colorbar;
ax.Box = 'on';
ax.TickDir = 'out';
ax.LineWidth = 2;

% uncomment for narrow colorbar
% axpos = get(gca,'Position');
% cpos = get(c,'Position');
% cpos(3) = 0.5*cpos(3);
% set(c,'Position',cpos);
% set(gca,'Position',axpos);

end

