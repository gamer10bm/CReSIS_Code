clc
clearvars -except gRadar 
% close all

if ~exist('gopro_data','var')
  %% Input checks
  params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

  % Example to run a specific segment and frame by overriding parameter spreadsheet values
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20200128_01');

  %% Load processed gopro data
  gopro_data = {};

  for pid = 1:length(params)
    if params(pid).cmd.generic
      param = params(pid);
      gopro_data{end+1} = load(ct_filename_support(param,'','gopro'));
    end
  end
end

%% Load the gopro images and make lat and lon map
gopro_ids = 345:2:355;%length(gopro_data{1}.path);

%Load paths elevs, lats, and lons
gopro_paths = gopro_data{1}.path(gopro_ids);
gopro_elevs = gopro_data{1}.elev(gopro_ids);
gopro_h_agls = gopro_data{1}.h_agl(gopro_ids);
gopro_headings = rad2deg(gopro_data{1}.heading(gopro_ids));
gopro_rolls = rad2deg(gopro_data{1}.roll(gopro_ids));
gopro_pitchs = rad2deg(gopro_data{1}.pitch(gopro_ids));
gopro_lats = gopro_data{1}.lat(gopro_ids);
gopro_lons = gopro_data{1}.lon(gopro_ids);

%make gopro structure
gopro_struct = struct('path',gopro_paths,'elev',gopro_elevs,'h_agl',gopro_h_agls,...
  'heading',gopro_headings,'roll',gopro_rolls,'pitch',gopro_pitchs,'lat',gopro_lats,...
  'lon',gopro_lons);

mosaic_plot(gopro_struct)

% %Load gopro_images and gopro_infos based on gopro_paths
% gopro_images = {}; gopro_infos = {};
% for gid = 1:length(gopro_paths)
%   %Load the images
%   gopro_images{end+1} = imread(gopro_paths{gid});
%   gopro_infos{end+1} = imfinfo(gopro_paths{gid});
% end
% 
% 
% lat_mosaic = []; lon_mosaic =[]; z_mosaic = []; cdat_mosaic = uint8([]);
% for mid = fliplr(1:length(gopro_images))
%   latgrid = []; longrid = []; cdata = [];
%   %debug
%   if 0
%     figure(mid)
%     imagesc(rot90(gopro_images{mid},2))
%     title('no changes')
%   end
%   %Get image state values
%   lat = gopro_lats(mid); lon = gopro_lons(mid); heading = gopro_headings(mid);
%   pitch_ang = gopro_pitchs(mid); roll_ang = gopro_rolls(mid);
%   %Get the x and y resolutions
%   ximg = convlength(gopro_infos{mid}.XResolution,'in','m'); %m
%   yimg = convlength(gopro_infos{mid}.YResolution,'in','m'); %m
% 
%   %Get the focal length
%   foc_leng = gopro_infos{mid}.DigitalCamera.FocalLength; %cm?
%   fin35 = gopro_infos{mid}.DigitalCamera.FocalLengthIn35mmFilm;
% 
%   %Get height above ground
%   H_agl_nadir = gopro_h_agls(mid); %m
% 
%   %Compensate for camangle
%   cam_ang = 5;
% %   H_agl_nadir = H_agl_nadir/cosd(cam_ang);
%   %Compensate for pitch
%   H_agl_y = H_agl_nadir/cosd(pitch_ang);
%   %Compensate for roll
%   H_agl_x = H_agl_nadir/cosd(roll_ang);
% 
%   %Convert X and Y to ground distances
%   xgrnd = ximg*H_agl_x/foc_leng; %m
%   ygrnd = yimg*H_agl_y/foc_leng; %m
% 
%   %FOV
%   x_fov = 122.6; %degrees
%   y_fov = 94.4; %degrees
% 
%   %Get the image origin indexes and the
%   img_rows = size(gopro_images{mid},1); img_cols = size(gopro_images{mid},2);
%   center_row = floor(img_rows./2);
%   center_col = floor(img_cols./2);
%   top_rows = 1:center_row; bot_rows = center_row:img_rows;
%   left_col = 1:center_col; rght_col = center_col:img_cols;
% 
%   %Make matrix for x-offsets and y-offsets (to the right and to the top are
%   %positive distances [will represent distance north and distance east])
%   topmaxdist = ygrnd*length(top_rows)./img_rows;
%   topdistvec = fliplr(linspace(0,topmaxdist,length(top_rows)));
%   botmaxdist = topmaxdist-ygrnd;
%   botdistvec = linspace(0,botmaxdist,length(bot_rows));
%   topfovvec = fliplr(linspace(0,y_fov/2,length(top_rows)));
%   botfovvec = linspace(0,-y_fov/2,length(bot_rows));
% 
%   ydistvec = [topdistvec botdistvec(2:end)];
%   yfovvec = [topfovvec botfovvec(2:end)];
% 
%   rghtmaxdist = xgrnd*length(rght_col)./img_cols;
%   leftmaxdist = rghtmaxdist-xgrnd;
%   leftdistvec = fliplr(linspace(0,leftmaxdist,length(left_col)));
%   rghtdistvec = linspace(0,rghtmaxdist,length(rght_col));
%   leftfovvec = fliplr(linspace(0,x_fov/2,length(left_col)));
%   rghtfovvec = linspace(0,-x_fov/2,length(rght_col));
% 
%   xdistvec = [leftdistvec rghtdistvec(2:end)];
%   xfovvec = [leftfovvec rghtfovvec(2:end)];
%   
%   %Adjust distvecs based on associated fov angle
%   [xdistgrid, ydistgrid] = meshgrid(xdistvec,ydistvec);
%   if 1
%     [xfovgrid, yfovgrid] = meshgrid(xfovvec,yfovvec);
%     xdistgrid = xdistgrid./cosd(xfovgrid);
%     ydistgrid = ydistgrid./cosd(yfovgrid);
%   else
%     xfovgrid = zeros(size(xdistgrid));
%     yfovgrid = zeros(size(ydistgrid));
%   end
% 
%   %Compensate for the heading
%   while abs(heading)>=360
%     if heading < 0
%       heading = heading + 360;
%     else
%       heading = heading - 360;
%     end
%   end
%   %Do switching for different quadrants
%   if (heading > 270 && heading < 360) || (heading < 0 && heading > -90)
%     %Quadrant -x +y
%     yn_cof = 1; xn_cof = 1;
%     ye_cof = -1; xe_cof = 1;
%   elseif (heading > 0 && heading <90) || (heading > -360 && heading < -270)
%     %Quadrant +x +y
%     yn_cof = 1; xn_cof = -1;
%     ye_cof = 1; xe_cof = 1;
%   elseif (heading > 90 && heading < 180) || (heading > -270 && heading < -180)
%     %Quadrant +x -y
%     yn_cof = 1; xn_cof = -1;
%     ye_cof = 1; xe_cof = 1;
%   else
%     %Quadrant -x -y
%     yn_cof = 1; xn_cof = 1;
%     ye_cof = -1; xe_cof = 1;
%   end
%   
%   northdistgrid = yn_cof*ydistgrid.*cosd(-heading) + xn_cof*xdistgrid.*sind(-heading);
%   eastdistgrid = ye_cof*ydistgrid.*sind(-heading) + xe_cof*xdistgrid.*cosd(-heading);
% 
%   %Estimate latitude grid values
%   Re_pol = 6357e3; Re_eq = 6378e3; %m Earth Radius at poles and equator
%   Ce_pol = Re_pol*2*pi; Ce_eq = Re_eq*2*pi; %m Earth Circumference
%   east_m2lat_deg = Ce_pol./360; %m/deg conversion from meters to degrees
%   latoffgrid = eastdistgrid./east_m2lat_deg; %degrees of latitude offset
%   latgrid = lat+latoffgrid;
% 
%   %Estimate longitude grid values
%   north_m2lon_deg = Ce_eq./360*cosd(lat);
%   lonoffgrid = northdistgrid./north_m2lon_deg;
%   longrid = lon+lonoffgrid;
%   
%   %Grab cdata
%   cdata = rot90(gopro_images{mid},2);
%   %Get inner image based on FOV limit
%   fovlim = 30;
%   if 1
%     fovlogmat = ones(size(latgrid));
%     xfovlog = xfovgrid<=fovlim;
%     yfovlog = yfovgrid<=fovlim;
%     fovlogmat = (xfovlog+yfovlog)==2;
%     latgrid(~fovlogmat) = nan;
%     longrid(~fovlogmat) = nan;
%     clogmat(:,:,1) = fovlogmat;
%     clogmat(:,:,2) = fovlogmat;
%     clogmat(:,:,3) = fovlogmat;
%     cdata(~clogmat) = nan;
%   end
%     
%   %Reduce image fidelity
%   if 1
%     subrows = floor(linspace(1,size(latgrid,1),500));
%     subcols = floor(linspace(1,size(latgrid,2),500));
%     latgrid = latgrid(subrows,subcols);
%     longrid = longrid(subrows,subcols);
%     cdata = cdata(subrows,subcols,:);
%   end
%   zgrid = zeros(size(latgrid));
%   
%   %Send to mosaic variables for surf plotting
%   lat_mosaic(end+1:end+size(latgrid,1),:) = latgrid;
%   lon_mosaic(end+1:end+size(latgrid,1),:) = longrid;
%   z_mosaic(end+1:end+size(latgrid,1),:) = zgrid;
%   cdat_mosaic(end+1:end+size(latgrid,1),:,:) = cdata;
%   %Add nan column for separating the different images
%   nanvec = nan(1,size(latgrid,2)); nanmat = nan([1 size(latgrid,2) size(cdat_mosaic,3)]);
%   lat_mosaic(end+1,:) = nanvec;
%   lon_mosaic(end+1,:) = nanvec;
%   z_mosaic(end+1,:) = nanvec;
%   cdat_mosaic(end+1,:,:) = nanmat;
%   
%   if 0
%     figure(100)
%     hold on
%     surf(longrid,latgrid,zgrid,cdata,'edgecolor','none')
%     hold off
%     set(gca,'YDir','normal')
%     view(2)
%     keyboard
%   end
% end
% 
% %Plot the completed mosaic
% if 1
%   figure(10+mid)
%   clf
%   geoshow(lat_mosaic,lon_mosaic,uint8(cdat_mosaic))
%   grid on
% %   title('Lat lon conversion')
%   xlabel(sprintf('Longitude ^o N'))
%   ylabel(sprintf('Latitude ^o E'))
% end

