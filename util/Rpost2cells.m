function Rcells = Rpost2cells(Rpost)
%RPOST2CELLS Rpost2cells converts the geospatial referencing object R from
%type 'postings' to type 'cells'
%   Geographic data referenced to a regular grid spacing can be of type
%   'postings', corresponding to an interpretation that the
%   latitude/longitude or planar x/y coordinates represent the location of
%   the data i.e. the centroid of each grid cell, or of type 'cell',
%   corresponding to an intepretation that the latitude/longitude or planar
%   x/y coordinates represent  E-W and N-S coordinates of the grid cell
%   edges, for example as is commonly the case for image-based data (e.g.
%   MODIS satellite imagery).

%   This function converts the Matlab spatial referencing object R from
%   type 'postings' to type 'cells'

%   Author: Matt Cooper, guycooper@ucla.edu, June 2019
%   Citation: Matthew Cooper (2019). matrasterlib
%   (https://www.github.com/mguycooper/matrasterlib), GitHub. Retrieved MMM
%   DD, YYYY.

refmat              =   worldFileMatrixToRefmat(Rpost.worldFileMatrix);
Zsize               =   Rpost.RasterSize;
Rcells              =   refmatToMapRasterReference(refmat,Zsize);
end

