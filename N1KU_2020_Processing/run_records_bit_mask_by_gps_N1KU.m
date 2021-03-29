% Script for updating records.bit_mask based on segments, frms, and gps
% ranges input then appends to records file and saves it.
new_bit_val = 4;
segs = {}; gps_ints = {}; inptfrmt = 'yyyy-MM-dd HH:mm:ss.SS';
segs{end+1} ='20200128_01'; gps_ints{end+1} = {'2020-01-28 15:13:07.01', '2020-01-28 15:13:19.17'};
segs{end+1} ='20200128_01'; gps_ints{end+1} = {'2020-01-28 15:15:43.49', '2020-01-28 15:17:26.85'};
segs{end+1} ='20200128_03'; gps_ints{end+1} = {'2020-01-28 17:11:11.07', '2020-01-28 17:11:21.74'};
segs{end+1} ='20200128_04'; gps_ints{end+1} = {'2020-01-28 20:31:12.24', '2020-01-28 20:31:20.50'};
segs{end+1} ='20200128_04'; gps_ints{end+1} = {'2020-01-28 20:31:42.26', '2020-01-28 20:31:48.98'};
segs{end+1} ='20200128_05'; gps_ints{end+1} = {'2020-01-28 21:21:32.33', '2020-01-28 21:21:42.68'};
segs{end+1} ='20200129_01'; gps_ints{end+1} = {'2020-01-29 19:10:03.02', '2020-01-29 19:10:12.25'};
segs{end+1} ='20200129_02'; gps_ints{end+1} = {'2020-01-29 20:08:30.51', '2020-01-29 20:08:32.86'};
segs{end+1} ='20200201_01'; gps_ints{end+1} = {'2020-02-01 15:14:03.50', '2020-02-01 15:14:14.81'};
segs{end+1} ='20200209_01'; gps_ints{end+1} = {'2020-02-09 22:57:02.14', '2020-02-09 23:01:14.14'};
segs{end+1} ='20200209_01'; gps_ints{end+1} = {'2020-02-09 23:04:38.63', '2020-02-09 23:04:41.35'};
segs{end+1} ='20200201_02'; gps_ints{end+1} = {'2020-02-01 16:34:28.58', '2020-02-01 16:34:37.97'};
segs{end+1} ='20200201_02'; gps_ints{end+1} = {'2020-02-01 16:52:20.87', '2020-02-01 16:56:09.35'};
segs{end+1} ='20200209_02'; gps_ints{end+1} = {'2020-02-09 23:15:46.33', '2020-02-09 23:16:05.96'};



%Load params
params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

%Mark segs for processing
params = ct_set_params(params,'cmd.generic',0);
for pid = 1:length(params)
  for sid = 1:length(segs)
    if strcmp(segs{sid},params(pid).day_seg)
      %Set generic to one
      params = ct_set_params(params,'cmd.generic',1,'day_seg',segs{sid});
      %Update notes with gps_interval cells
      newnotes = {};
      if ~isfield(params(pid).cmd,'notes') || ~iscell(params(pid).cmd.notes)
        newnotes = {gps_ints{sid}};
      else
        newnotes = params(pid).cmd.notes;
        newnotes{end+1} = gps_ints{sid};
      end
      params = ct_set_params(params,'cmd.notes',newnotes,'day_seg',segs{sid});
    end
  end
end

% params = ct_set_params(params,'reset_override',true);
%% Automated section
% =====================================================================
%Formatting 
cmd_size = matlab.desktop.commandwindow.size;
line_width = cmd_size(1); lid = 1; line_str = '=';
while lid <line_width
  lid = lid +1; line_str = [line_str line_str(1)];
end
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  if isfield(param.cmd,'notes') && (~isempty(param.cmd.notes) || ~iscell(param.cmd.notes))
    %% Step 1: convert gps times ranges to records and make records index vector
    cmdnotes = param.cmd.notes;
    rec_vec = [];
    for cmd_id = 1:length(cmdnotes)
      %Convert gps_ints strings to gps_time
      gps_times = posixtime(datetime(cmdnotes{cmd_id},'InputFormat',inptfrmt));
      %Get the recs based on gps time
      [day_seg, frm_id, rec_ids] = get_frame_id(param,gps_times);
      %Make records id vec
      add_strt = floor(min(rec_ids)); add_end = ceil(max(rec_ids));
      rec_wind = 300; %floor(.2*length(rec_vec_add));
      if new_bit_val == 0 && exist('rec_wind','var')
        warning(sprintf('\n\tExtending record windows by %.0f',rec_wind*2))
        add_strt = min(rec_vec_add)-rec_wind; add_end = max(rec_vec_add)+rec_wind;
      end
      %Make append vector
%       if add_strt < 20; if new_bit_val==0; add_strt = 1; else add_strt = 2; end; end;
      if add_strt < 1; add_strt = 1; end
      rec_vec_add = add_strt:add_end;
      
      rec_vec = [rec_vec rec_vec_add];
    end
    %% Step 2: Load records file, set bit mask to two, then save
    %Load the records file
    recs_fn = ct_filename_support(param,'','records');
    fprintf('%s\nLoading records for %s...\n\n',line_str,param.day_seg)
    recs = load(recs_fn);
    %Set bit mask for rec range the current bit_mask
    bit_mask = recs.bit_mask;
    if isfield(param,'reset_override') && param.reset_override
      warning(sprintf('Resetting entire bit_mask of segment %s\n',param.day_seg));
      bit_mask(:) = 0;
    else
      rec_vec = rec_vec(rec_vec<length(bit_mask));
      bit_mask(rec_vec) = new_bit_val;
    end
    %Set the param_records radar.wfs.bad_value to nan
    param_records = ct_set_params(recs.param_records,'radar.wfs.bad_value',nan);
    %Save using append
    fprintf('Appending bitmask:\n\t%.0f/%.0f updated\n%s\n',length(rec_vec),length(bit_mask),line_str);
    ct_save(recs_fn,'-append','bit_mask','param_records')
  end
end