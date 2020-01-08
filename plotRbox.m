function [h,coords] = plotRbox(R,varargin)
%PLOTRBOX plotRbox plots a box corresponding to the geographic limits of
%the map/geographic raster reference object R. 

%   Author: Matt Cooper, guycooper@ucla.edu, May 2019
%   Citation: Matthew Cooper (2019). matrasterlib
%   (https://www.github.com/mguycooper/matrasterlib), GitHub. Retrieved MMM
%   DD, YYYY

% confirm mapping toolbox is installed
assert(license('test','map_toolbox')==1, ...
                        'rasterinterp requires Matlab''s Mapping Toolbox.')
                    
% confirm R is either MapCells or GeographicCellsReference objects
validateattributes(R, ...
                        {'map.rasterref.MapCellsReference', ...
                        'map.rasterref.GeographicCellsReference', ...
                        'map.rasterref.MapPostingsReference', ...
                        'map.rasterref.GeographicPostingsReference'}, ...
                        {'scalar'}, 'rasterinterp', 'R', 2)
                    
% determine if R is planar or geographic and call the appropriate function
if strcmp(R.CoordinateSystemType,'planar')
    [h,coords]      =   mapplotRbox(R,varargin{:});
elseif strcmp(R.CoordinateSystemType,'geographic')
    [h,coords]      =   geoplotRbox(R,varargin{:});
end

%% apply the appropriate function

    function [h,coords] = mapplotRbox(R,varargin)
        
        uplefty     =   R.YWorldLimits(2);
        upleftx     =   R.XWorldLimits(1);
        uprightx    =   R.XWorldLimits(2);
        uprighty    =   R.YWorldLimits(2);
        lowlefty    =   R.YWorldLimits(1);
        lowleftx    =   R.XWorldLimits(1);
        lowrighty   =   R.YWorldLimits(1);
        lowrightx   =   R.XWorldLimits(2);
        coords.y    =   [uprighty uplefty lowlefty lowrighty uprighty];
        coords.x    =   [uprightx upleftx lowleftx lowrightx uprightx];
        h           =   plot(coords.x,coords.y,varargin{:});
    end

    function [h,coords] = geoplotRbox(R,varargin)
        
        uplefty     =   R.LatitudeLimits(2);
        upleftx     =   R.LongitudeLimits(1);
        uprightx    =   R.LongitudeLimits(2);
        uprighty    =   R.LatitudeLimits(2);
        lowlefty    =   R.LatitudeLimits(1);
        lowleftx    =   R.LongitudeLimits(1);
        lowrighty   =   R.LatitudeLimits(1);
        lowrightx   =   R.LongitudeLimits(2);
        coords.y    =   [uprighty uplefty lowlefty lowrighty uprighty];
        coords.x    =   [uprightx upleftx lowleftx lowrightx uprightx];
        h           =   plot(coords.x,coords.y,varargin{:});
    end
end

