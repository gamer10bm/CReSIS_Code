function [] = gopro_postprocess(gopro_output)

% this function is used to save all processing outputs and for filling in
% missing segments

%% General Setup
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, gopro_output{1}(1).param.season_name, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks


%% Processing
%Initialize postprocessing save variables
path_all = {}; lat_all = []; lon_all = []; elev_all = []; id_out = 1;
date_all = []; head_all = []; pitch_all = []; roll_all = []; hagl_all = [];
postproc_log = zeros(1,length(gopro_output));
for id_gp = 1:length(gopro_output)
  proc_now = gopro_output{id_gp};
  %Check if output is empty
  if isempty(proc_now.paths)
    postproc_log(id_gp) = 1;
  else
    for id_load = 1:length(proc_now.paths)
      path_all{id_out} = proc_now.paths{id_load};
      date_all(id_out) = proc_now.date(id_load);
      lat_all(id_out) = proc_now.lat(id_load);
      lon_all(id_out) = proc_now.lon(id_load);
      elev_all(id_out) = proc_now.elev(id_load);
      head_all(id_out) = proc_now.heading(id_load);
      pitch_all(id_out) = proc_now.pitch(id_load);
      roll_all(id_out) = proc_now.roll(id_load);
      hagl_all(id_out) = proc_now.h_agl(id_load);
      id_out = id_out +1;
    end
  end
  proc_now = [];
end

% Post processing for empty segments
chk_num = 500; %Needs to be set in run_gopro
for id_post = 1:length(gopro_output)
  proc_now = gopro_output{id_post};
  %Check if empty
  if postproc_log(id_post)
    rec_now = proc_now.records;
    ids_short = floor(linspace(1,length(rec_now.gps_time),chk_num));
    lat_short = rec_now.lat(ids_short); lon_short = rec_now.lon(ids_short);
    %Search for closest gopro photo
    path = {}; pic_time = []; recs = struct('id',[],'frm',[],'lat',[],...
    'lon',[],'elev',[],'heading',[],'pitch',[],'roll',[],'h_agl',[]);
    gps_time =[]; lat =[]; lon=[]; elev=[]; frm = [];
    heading= []; pitch = []; roll = []; h_agl = [];
    id_rec = 1;
    for id_chk = 1:chk_num
      lat_diff = abs(lat_all-lat_short(id_chk));
      lon_diff = abs(lon_all-lon_short(id_chk));
      [min_diff, min_id] = min(sqrt(lat_diff.^2+lon_diff.^2));
      if min_diff < .5
        rec_id(id_rec) = ids_short(id_chk);
        path{id_rec} = path_all{min_id};
        pic_time(id_rec) = date_all(min_id);
        gps_time(id_rec) = rec_now.gps_time(ids_short(id_chk));
        lat(id_rec) = lat_all(min_id);
        lon(id_rec) = lon_all(min_id);
        elev(id_rec) = elev_all(min_id);
        frm(id_rec) = rec_now.frms(ids_short(id_chk));
        heading(id_rec) = head_all(min_id);
        pitch(id_rec) = pitch_all(min_id);
        roll(id_rec) = roll_all(min_id);
        h_agl(id_rec) = hagl_all(min_id);
        %Load state to records structure
        recs(id_rec) = struct('id',rec_id,'frm',frm(id_rec),'lat', lat(id_rec),...
          'lon', lon(id_rec), 'elev', elev(id_rec), 'heading',heading(id_rec),...
          'pitch',pitch(id_rec), 'roll',roll(id_rec),'h_agl',h_agl(id_rec));
        %Update id_rec
        id_rec = id_rec +1;
      end
    end
    if 0 %debug plot
      figure(1)
      plot(lat_short,lon_short)
      hold on
      plot(lat,lon,'ro','MarkerFaceColor','r')
      hold off
      legend({'Records','Gopro'})
      keyboard
    end
    fprintf('%d files saved for %s\n',length(rec_id),proc_now.param.day_seg)
    %Save the ouput
    out_path = ct_filename_support(proc_now.param,'','gopro');
    ct_save(out_path,'path','pic_time','recs','gps_time','lat','lon','elev',...
      'heading','pitch','roll','h_agl')
    fprintf('Saved gopro info: %s (%s)\n',out_path,datestr(now))
  end
end

end
    