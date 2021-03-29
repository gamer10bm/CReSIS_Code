function [output] = gopro_process(param, param_override)
% function [output] = process_gopro(param, param_override)

% The purpose of this function is to grab any gopro pictures and save 
% useful information by each day and segment turned on.

%% General Setup
c = physconst('Lightspeed'); %m/s
%Check if preprocessed
gopro_tmp_path = ct_filename_ct_tmp(param,fullfile('gopro',param.season_name,'data.mat'));
if ~exist(gopro_tmp_path,'file') || 1
  gopro_pics = gopro_preprocess(param);
else
  gopro_tmp = load(gopro_tmp_path);
  gopro_pics = struct('path',gopro_tmp.path,'date',gopro_tmp.date,'gpsinfo',gopro_tmp.gpsinfo);
end

% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isfield(param.records.file,'base_dir') || isempty(param.records.file.base_dir)
  error('Please specify the field param_override.records.file.base_dir')
end

%% Processing
%Get layer directory
laydir = ct_filename_out(param,'','layer'); layfnfmt = 'Data_%s_%03.0f.mat';
%Grab gps data and determine lower and upper bounds of time
seg_records = load(ct_filename_support(param,'','records'));
seg_gps_time = datetime(seg_records.gps_time,'convertfrom','posixtime')-hours(6); %Time is GMT %Hardcode fix
seg_records.frms = geodetic_to_along_track(seg_records.lat, seg_records.lon)./param.records.frames.length + 1;
low_time = seg_gps_time(1); 
upp_time = seg_gps_time(end);

%Get all gopro pictures for the given day segment
path = {}; pic_time = [];
recs = struct('id',[],'frm',[],'lat',[],'lon',[],'elev',[],'heading',[],...
  'pitch',[],'roll',[],'h_agl',[]);
gps_time =[]; lat =[]; lon=[]; elev=[]; frm = [];
heading= []; pitch = []; roll = []; h_agl = [];
id_rec = 1;
for id_pic = 1:length(gopro_pics)
  %Check if creation date falls on day_seg
  if isbetween(gopro_pics(id_pic).date,low_time,upp_time)
    %Load pic data
    path{id_rec} = gopro_pics(id_pic).path;
    %Load the pic time
    if isfield(gopro_pics(id_pic).gpsinfo,'GPSDateStamp')
      datestamp = gopro_pics(id_pic).gpsinfo.GPSDateStamp;
      timestamp = gopro_pics(id_pic).gpsinfo.GPSTimeStamp;
      timeraw = datetime(str2num(datestamp(1:4)), str2num(datestamp(6:7)), ...
        str2num(datestamp(9:end)), timestamp(1), timestamp(2), timestamp(3));
    else
      timeraw = gopro_pics(id_pic).date + hours(6);
    end
    timeraw = timeraw + seconds(15); %Hardcode
    pic_time(id_rec) = posixtime(timeraw); %In posixtime
    
    %Find closest index for gps_time
    [~, rec_id] = min(abs(timeraw-hours(6)-seg_gps_time));
    %Load state vectors
    gps_time(id_rec) = posixtime(seg_gps_time(rec_id)+hours(6)); %In posixtime
    lat(id_rec) = seg_records.lat(rec_id); lon(id_rec) = seg_records.lon(rec_id);
    elev(id_rec) = seg_records.elev(rec_id);
    heading(id_rec) = seg_records.heading(rec_id);
    pitch(id_rec) = seg_records.pitch(rec_id);
    roll(id_rec) = seg_records.roll(rec_id);
    frm(id_rec) = seg_records.frms(rec_id);
    
    
    %Load in gps info from gopro pic
    if isfield(gopro_pics(id_pic).gpsinfo,'GPSAltitude')
      %Adjust elevation lat and lon based on gopro information
      elev(id_rec) = gopro_pics(id_pic).gpsinfo.GPSAltitude;

      latval = gopro_pics(id_pic).gpsinfo.GPSLatitude;
      lat(id_rec) = latval(1)+latval(2)./60+latval(3)./3600;
      if strcmpi(gopro_pics(id_pic).gpsinfo.GPSLatitudeRef,'S')
        lat(id_rec) = lat(id_rec)*-1;
      end

      lonval = gopro_pics(id_pic).gpsinfo.GPSLongitude;
      lon(id_rec) = lonval(1)+lonval(2)./60+lonval(3)./3600;
      if strcmpi(gopro_pics(id_pic).gpsinfo.GPSLongitudeRef,'W')
        lon(id_rec) = lon(id_rec)*-1;
      end
    end
    
    %Grab layer data and determine height AGL
    layfn = fullfile(laydir,sprintf(layfnfmt,param.day_seg,frm(id_rec)));
    if exist(layfn,'file')
      layer = load(layfn);
      Rtosurface = layer.twtt(1,:)*c/2;
      id_lay = find(gps_time(id_rec)<=layer.gps_time,1);
      windwidth = 5;
      ind_wind = (id_lay-floor(windwidth/2)):(id_lay+ceil(windwidth/2));
      if ind_wind(1) < 1
        ind_wind = 1:windwidth;
      elseif ind_wind(end) > length(Rtosurface)
        ind_wind = (length(Rtosurface)-ind_width):length(Rtosurface);
      end
      h_agl(id_rec) = mean(Rtosurface(ind_wind));
    else
      layer = [];
      h_agl(id_rec) = convlength(1500,'ft','m');
    end
    
    if 1
      %Adjust lat and lon based on assumed cam angle to horizon (cam tilt + pitch) and elevation based on
      %records data
      camangle = deg2rad(0.5); R = 6371e3; %m Earth Radius
      horzangle = camangle+pitch(id_rec);
      dist_along = tan(horzangle)*h_agl(id_rec); %m
      lat1 = seg_records.lat(rec_id); lat2 = seg_records.lat(rec_id+1);
      lon1 = seg_records.lon(rec_id); lon2 = seg_records.lon(rec_id+1);
      bearing = atan2(sind(lon2-lon1)*cosd(lat2), cosd(lat1)*sind(lat2)-...
        sind(lat1)*cosd(lat2)*cosd(lon2-lon1)); %radians
      
      %Change heading to bearing value
      heading(id_rec) = bearing;
      
      %Make lat1 and lon1 the updated value
      latnew = asind(sind(lat1)*cos(dist_along/R)+cosd(lat1)*sin(dist_along/R)*cos(bearing));
      lonnew = lon1 + atan2d(sin(bearing)*sin(dist_along/R)*cosd(lat1), cos(dist_along/R)-sind(lat1)*sind(latnew));
      
      %Adjust lat and lon based on roll angle 
      if 1
        lon1 = lonnew; lat1 = latnew;
        dist_cross = tan(roll(id_rec))*h_agl(id_rec); %m
        bearing = bearing+pi/2; %rad changes bearing direction to be out of the right wing

        latnew = asind(sind(lat1)*cos(dist_cross/R)+cosd(lat1)*sin(dist_cross/R)*cos(bearing));
        lonnew = lon1 + atan2d(sin(bearing)*sin(dist_cross/R)*cosd(lat1), cos(dist_cross/R)-sind(lat1)*sind(latnew));
      end
      lat(id_rec) = latnew;
      lon(id_rec) = lonnew;
      
      %Adjust gps_time and pic_time based on instantaneous velocity and
        %distance change
      a = sind((lat1-lat2)/2)^2+cosd(lat1)*cosd(lat2)*sind((lon1-lon2)/2)^2;
      b = 2*atan2(sqrt(a),sqrt(1-a));
      dinst = R*b;%m
      vel_est = dinst/(seg_records.gps_time(rec_id+1)-seg_records.gps_time(rec_id)); %m/s
      time_off = dist_along./vel_est; %seconds
      
      timeraw = timeraw + seconds(time_off);
      pic_time(id_rec) = posixtime(timeraw); %In posixtime
    
      %Find closest index for gps_time
      [~, rec_id] = min(abs(timeraw-hours(6)-seg_gps_time));
      %Update gps time
      gps_time(id_rec) = posixtime(seg_gps_time(rec_id)+hours(6)); %In posixtime      
    end
    %Load state to records structure
    recs(id_rec) = struct('id',rec_id,'frm',frm(id_rec),'lat', lat(id_rec),...
      'lon', lon(id_rec), 'elev', elev(id_rec), 'heading',heading(id_rec),...
      'pitch',pitch(id_rec), 'roll',roll(id_rec),'h_agl',h_agl(id_rec));
    
    %Update id_rec
    id_rec = id_rec +1;
  end
end

%Sort the gopro fields based on pic_time
if ~isempty(pic_time)
  [pic_time, ind_sort] = sort(pic_time);
  path = path(ind_sort);
  gps_time = gps_time(ind_sort);
  lat = lat(ind_sort);
  lon = lon(ind_sort);
  elev = elev(ind_sort);
  heading = heading(ind_sort);
  pitch = pitch(ind_sort);
  roll = roll(ind_sort);
  h_agl = h_agl(ind_sort);
  recs = recs(ind_sort);
end

fprintf('%d files saved for %s\n',length(recs),param.day_seg)
%% Save results
out_path = ct_filename_support(param,'','gopro');
ct_save(out_path,'path','pic_time','recs','gps_time','lat','lon','elev',...
  'heading','pitch','roll','h_agl')
fprintf('Saved gopro info: %s (%s)\n',out_path,datestr(now))

%% Output useful info
output.param = param; output.paths = path; output.date = pic_time;
output.lat = lat; output.lon = lon; output.elev = elev; output.h_agl = h_agl;
output.heading = heading; output.pitch = pitch; output.roll = roll;
output.records = seg_records;


end