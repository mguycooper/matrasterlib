function [ h,ax,c ] = mapshowsurface( Z,R,varargin )
%MAPSHOWSURFACE alias for MAPSHOW with default "DisplayType" set to
%"surface", adds a colorbar, sets axis tight, box on, tickdir out, and
%linewidth = 2

%   This function is an alias for the in-built Matlab function
%   mapshow.m (Copyright 2016 The MathWorks, Inc.), designed to map 
%   2-d gridded data by default with a 'Surface' display type, a
%   colorbar, and some improved graphics the user can adjust ex-post

%   Author: Matt Cooper, guycooper@ucla.edu, May 2019

%   Syntax

%   MAPSHOWSURFACE(...) displays a regular data grid Z and appends a colorbar to
%   the current axes in the default (right) location.  Z is a 2-D array of
%   class double and R is a referencing matrix or a map raster reference
%   object that relates the subscripts of Z to map coordinates. DISPLAYTYPE
%   is set to 'surface' and the axis is set tight.

%   H = MAPSHOWSURFACE(...) returns a handle to a MATLAB graphics object
%   [H,C] = MAPSHOWSURFACE(...) returns a handle to a MATLAB graphics
%   object and a MATLAB colorbar object

%   MAPSHOWSURFACE(..., Name, Value) specifies name-value pairs that set
%   MATLAB graphics properties. Parameter names can be abbreviated and are
%   case-insensitive. Refer to the MATLAB Graphics documentation on surface
%   for a complete description of these properties and their values.

%   See also geoshow, mapshow

h = mapshow(Z,R,'DisplayType','surface',varargin{:}); ax = gca;
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

