function R = rasterref(X,Y,rasterInterpretation)

%RASTERREF R = rasterref(X,Y,rasterInterpretation) constructs spatial
%referencing object R from 2-d grids X,Y that define the E-W and N-S
%coordinates of each grid cell. The default interpretation is 'Postings'
%i.e. that the X,Y coordinate pairs represent the centroids of each grid
%cell. User can specify 'Cells' if the X,Y coordinate pairs represent the
%edges of the grid cells. X and Y can be geographic or projected (planar)
%coordinates.

%   This function exists because many geospatial datasets are provided with
%   vectors or grids of latitude/longitude and/or planar coordinate values
%   that represent the centroid or cell edges of each grid-cell, but many
%   Matlab Mapping Toolbox functions require the spatial referencing object
%   R as input. This function creates the object R from the
%   latitude/longitude or x/y coordinate vectors or grids. 

%   Author: Matt Cooper, guycooper@ucla.edu, June 2019 Citation: Matthew
%   Cooper (2019). matrasterlib: a software library for processing
%   satellite and climate model data in Matlab©
%   (https://www.github.com/mguycooper/matrasterlib), GitHub. Retrieved MMM
%   DD, YYYY.

%   Syntax

%   R = rasterref(X,Y)
%   R = rasterref(LON,LAT)
%   R = rasterref(__,rasterInterpretation)

%   Description

%   R = rasterref(X,Y) constructs spatial referencing object R from 2-d
%   grids of E-W (longitude) and N-S (latitude) coordinates X,Y. X and Y
%   are 2-d numeric matrices (grids) of equal size oriented E-W and N-S
%   such that the linear (row,col) index (1,1) is the northwest corner grid
%   cell and the coordinate pair X(1,1),Y(1,1) is (by default) intepreted
%   as the centroid of the northwest corner grid cell. This is equivalent
%   to a 'Postings' cell interpretation. X defines the E-W (longitude or
%   projected coordinate) value of every grid cell centroid and Y defines
%   the N-S (latitude or projected coordinate) value of every grid cell
%   centroid. The default 'Postings' assumption is consistent with an
%   interpretation of raster grid cell values representing the centroid of
%   each grid cell, for example as would be provided in a netcdf or h5
%   scientific raster dataset or if using [X,Y]=meshgrid(x,y) where x and y
%   are vectors that span min(x):x_cell_extent:max(x) and
%   min(y):y_cell_extent:max(y) where x and y define the coordinates of
%   data values (not cell edges). X and Y can be geographic or planar
%   coordinate systems.

%   R = rasterref(X,Y,'Cells') applies the 'Cells' intepretation instead of
%   'Postings', which is consistent with an interpretation that X,Y define
%   the E-W and N-S coordinates of the grid cell edges, for example as
%   is commonly the case for image-based data (e.g. MODIS satellite
%   imagery).

%   See also rasterref, georefcells, maprefcells, meshgrid

%   Examples - forthcoming

%% Check inputs
                    
% confirm there are at least two and no more than 3 inputs
narginchk(2,3)

% if no rasterInterpretation is given, assign 'postings', otherwise make
% sure either 'cells' or 'postings' is given
if nargin == 2
    rasterInterpretation = 'postings';
elseif nargin == 3
    assert(strcmp(rasterInterpretation,'cells') | ...
            strcmp(rasterInterpretation,'postings'), ...
            ['Input argument 3, rasterInterpretation, must be either ', ...
            '"cells" or "postings"']);    
end

% confirm X and Y are 2d numeric grids of equal size
validateattributes(X, {'numeric'}, {'2d','size',size(Y)}, 'rasterref', 'X', 1)
validateattributes(Y, {'numeric'}, {'2d','size',size(X)}, 'rasterref', 'Y', 2)

% get an estimate of the grid resolution (i.e. the grid cell size)
xres                =   diff(X(1,:));
yres                =   diff(Y(:,1));

% check that the gridding is uniform. 
assert(all(xres(1) == xres), ...
            'Input argument 1, X, must have uniform grid spacing');
assert(all(yres(1) == yres), ...
            'Input argument 2, Y, must have uniform grid spacing');

% confirm X is oriented E-W and Y is N-S
assert(X(1,1)==min(X(:)) & X(1,2)>X(1,1), ...
                'Input argument 1, X, must be oriented E-W');
assert(Y(1,1)==max(Y(:)) & Y(1,1)>Y(2,1), ...
                'Input argument 2, Y, must be oriented N-S');

% determine if R is planar or geographic and call the appropriate function
tf                  =   islatlon(Y(1,1),X(1,1));
if tf
    R               =   rasterrefgeo(X,Y,rasterInterpretation);
else
    R               =   rasterrefmap(X,Y,rasterInterpretation);
end

%% apply the appropriate function

    function R = rasterrefmap(X,Y,rasterInterpretation)
        
        % NOTE: i convert to double because I have encountered X,Y grids
        % provided by netcdf files that are stored as type uint, but it
        % might be better to convert X and Y to single or double and
        % confirm what is most compatible with the functions called
        xlims       =   double([min(X(:)) max(X(:))]);
        ylims       =   double([min(Y(:)) max(Y(:))]);
        rasterSize  =   size(X);
        if strcmp(rasterInterpretation,'cells')
            R       =   maprefcells(xlims,ylims,rasterSize, ...
                            'ColumnsStartFrom','north', ...
                            'RowsStartFrom', 'west');
        elseif strcmp(rasterInterpretation,'postings')
            R       =   maprefpostings(xlims,ylims,rasterSize, ...
                            'ColumnsStartFrom','north', ...
                            'RowsStartFrom', 'west');
        end
    end

    function R = rasterrefgeo(X,Y,rasterInterpretation)
        
        % 'columnstartfrom','south' is default and corresponds to S-N
        % oriented grid as often provided by netcdf and h5 but for sanity 
        % I require this function accept N-S oriented grids i.e. index 
        % (1,1) is NW corner, consequently set 'columnstartfrom','north'
        
        xlims       =   double([min(X(:)) max(X(:))]);
        ylims       =   double([min(Y(:)) max(Y(:))]);
        rasterSize  =   size(X);
        
        if strcmp(rasterInterpretation,'cells')
            R       =   georefcells(ylims,xlims,rasterSize, ...
                            'ColumnsStartFrom','north', ...
                            'RowsStartFrom', 'west');
        elseif strcmp(rasterInterpretation,'postings')
            R       =   georefpostings(ylims,xlims,rasterSize, ...
                            'ColumnsStartFrom','north', ...
                            'RowsStartFrom', 'west');
        end
                        
        % note: R.CellExtentInLongitude should equal:
        % (R.LongitudeLimits(2)-R.LongitudeLimits(1))/R.RasterSize(2)
        % but unless pre-processing is performed on standard netcdf or h5
        % grid to adjust for edge vs centroid, a 'cells' interpretation
        % gets it wrong
    end
        
end

% Notes 

% for reference, the N-S E-W could be handled here but the built-in error
% message handling is less informative than the custom message option with
% using assert

% % confirm X and Y are 2d numeric grids of equal size oriented N-S and E-W
% validateattributes(X(1,:),  {'numeric'}, ... 
%                             {'2d','size',size(Y(1,:)),'increasing'}, ... 
%                             'R2grid', 'X', 1)
% validateattributes(Y(1,:),  {'numeric'}, ...
%                             {'2d','size',size(X(1,:)),'decreasing'}, ...
%                             'R2grid', 'Y', 2)

