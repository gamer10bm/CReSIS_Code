function [recdata] = init_view_frames(params,gisfn,tmpdays_offset,tmpfn_init_run_override)
%% Determine if tmp file exists or init_run
if ~exist('tmpdays_offset','var')
  tmpdays_offset = 2;
end

if ~exist('tmpfn_init_run_override')
  tmpfn_init_run_override = false;
end
fn_prefix = 'data';
tmpdir = ct_filename_ct_tmp(params(1),fullfile('selectframes',params(1).season_name));
tmpfn_srch=fullfile(tmpdir,sprintf('%s*.mat',fn_prefix));
tmpfn_tdy = fullfile(tmpdir,sprintf('%s_%s.mat',fn_prefix,strrep(date,'-','_')));
if exist(tmpfn_tdy,'file') && ~tmpfn_init_run_override
  tmpfn = tmpfn_tdy;
else
  %Check for existing data files
  currdiff = tmpdays_offset;
  tmpfns = dir(tmpfn_srch);
  for tmp_id = 1:length(tmpfns)
    tmp = tmpfns(tmp_id);
    tmpdate = tmp.name(length(fn_prefix)+2:end-length('.mat'));
    if ~isempty(tmpdate) && (datenum(tmpdate)-datenum(date()))<=currdiff
      %Grab valid file
      tmpfn = fullfile(tmpdir,tmp.name);
      currdiff = datenum(tmpdate)-datenum(date());
    else
      %Suggest deleting this file
      warning('We recommend deleting this file: copy and run the following lines')
      fprintf('\n\n%sRun this to delete the file:\n\tdelete(''%s'');\n\n,','%%',fullfile(tmpdir,tmp.name));
    end
  end
end

if ~exist('tmpfn','var') || tmpfn_init_run_override
  tmpfn = tmpfn_tdy;
  init_run = true;
else
  load(tmpfn);
  return
end

%Check for init_run
if exist(tmpfn,'file') && init_run
    delete(tmpfn)
end
if init_run
  %% Initialize data file for view_frames
  for param_idx = 1:length(params)
    params(param_idx).gisfn = gisfn;
    param = params(param_idx);
    rec_tmp = [];
    if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
      %Grab the records data
      fn = ct_filename_support(param,'','records');
      rec_tmp = load(fn);
      lontmp = []; lattmp = []; eltmp = []; gpstmp = [];
      try
        lontmp = rec_tmp.Longitude; lattmp = rec_tmp.Latitude; eltmp = rec_tmp.Elevation;
        gpstmp = rec_tmp.GPS_time;
      catch
        lontmp = rec_tmp.lon; lattmp = rec_tmp.lat; eltmp = rec_tmp.elev;
        gpstmp = rec_tmp.gps_time;
      end
      %Grab every nint^th element
      nint = 100;
      lontmp = lontmp(1:nint:end); lattmp = lattmp(1:nint:end); eltmp = eltmp(1:nint:end);
      gpstmp = gpstmp(1:nint:end);
      %Determine the frame number for each gpstime
      [~, frm_ids] = get_frame_id(param,gpstmp); 
      %Load to output structure
      if ~exist('recdata','var') || ~isstruct(recdata)
        recdata = struct('Longitude',lontmp,'Latitude',lattmp, ....
          'Elevation',eltmp,'GPS_time',gpstmp,'Frame_IDs',frm_ids,...
          'param',param,'filename',fn,'id',param.day_seg);
      else
        recdata(end+1) = struct('Longitude',lontmp,'Latitude',lattmp, ....
          'Elevation',eltmp,'GPS_time',gpstmp,'Frame_IDs',frm_ids,...
          'param',param,'filename',fn,'id',param.day_seg);
      end
      %Reinitialize
      rec_tmp=[];
    end
  end
  %Save the data each time
  ct_save(tmpfn,'recdata')
end
end