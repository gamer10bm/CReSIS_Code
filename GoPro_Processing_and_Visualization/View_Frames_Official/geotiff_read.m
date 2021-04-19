function [data, Xdata, Ydata, Cdata] = geotiff_read(geotiff_fn)
% [data, Xdata, Ydata, Cdata] = geotiff_read(geotiff_fn)
%
% Function for reading geotiff files and outputs the X, Y, and C data
% needed for imagesc
%
%
% Outputs:
%  data: A structure with all X, Y, and C data
%
% Examples:
%   geotiff_fn = ct_filename_gis('greenland/Landsat-7/mzl7geo_90m_lzw.tif');
%   geotiff_read(geotiff_fn);
%
% Author: John Paden

% Get the projection information
proj = geotiffinfo(geotiff_fn);

% Read the image
switch lower(proj.ColorType)
  case {'indexed'}
    [X, CMAP, R] = geotiffread(geotiff_fn);
    error('Not supported.');
  otherwise
    [RGB, R, tmp] = geotiffread(geotiff_fn);
end
switch lower(proj.ModelType)
  case {'modeltypeprojected'}
    R = R/1e3; % Convert to km
  case {'modeltypegeographic'}
    R = R; % lat/lon in degrees
  otherwise
    error('Not supported.');
end
if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
  RGB = uint8(RGB);
end
if strcmpi(class(RGB),'int16') || strcmpi(class(RGB),'single')
  RGB = double(RGB);
end
if 0
  % DEBUG
  RGB(RGB<0) = NaN;
end

Xdata = R(3,1) + R(2,1)*(1:size(RGB,2)); 
Ydata = R(3,2) + R(1,2)*(1:size(RGB,1));
Cdata = RGB;

data.X = Xdata; data.Y = Ydata; data.C = Cdata;

return;
