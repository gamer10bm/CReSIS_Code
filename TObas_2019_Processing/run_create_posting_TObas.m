% Script run_create_posting
%
% Loads the "post" worksheet from the parameter spreadsheet and then calls
% create_posting with this information.
%
% Authors: Theresa Stumpf, John Paden, Reece Mathews
%
% See also: create_posting.m

%% User Settings
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'));

%Set layer_dir to ops
params = ct_set_params(params,'post.layer_dir','ops'); %Empty for OPS by wiki
params = ct_set_params(params,'post.ops.layers',{'bottom','surface'});
params = ct_set_params(params,'post.ops.en',1);
layers = [struct('name', 'surface', 'source', 'ops')];
layers(end+1) = [struct('name', 'bottom', 'source', 'ops')];
params = ct_set_params(params,'post.layers',layers);

%Use qlook without coh_noise_removal
% params = ct_set_params(params,'post.data_dirs',{'qlook'});

%Use qlook with coh_noise_removal
% params = ct_set_params(params,'post.data_dirs',{'qlook_noise'});
% params = ct_set_params(params,'post.out_path','post_noise');

% % %Use qlook with deconv
% params = ct_set_params(params,'post.data_dirs',{'qlook_deconv'});
% params = ct_set_params(params,'post.out_path','post_deconv');

% Use standard after array processing
params = ct_set_params(params,'post.data_dirs',{'CSARP_post/standard'});
params = ct_set_params(params,'post.out_path','post_standard');

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
  '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
  '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
for did = 1:length(daysegs)
  params = ct_set_params(params,'cmd.generic',1,'day_seg',daysegs{did});
end
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20190129_01');
% params = ct_set_params(params,'cmd.frms',6);

%% Automated Section
% =====================================================================
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  create_posting(param,param_override);
end

%% Create season-wide concatenated and browse files (CSV and KML)
% =====================================================================
concatenate_csv_kml = false;
for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  if param.post.concat_en
    concatenate_csv_kml = true;
    param = params(param_idx);
    break;
  end
end
if concatenate_csv_kml
  fprintf('Creating season concatenated and browse files (%s)\n', datestr(now));
  if ~isempty(param.post.out_path)
    post_path = ct_filename_out(param,param.post.out_path,'',1);
  else
    post_path = ct_filename_out(param,'post','',1);
  end
  % Create concatenated and browse files for all data
  csv_base_dir = fullfile(post_path,'csv');
  kml_base_dir = fullfile(post_path,'kml');
  if ~exist(kml_base_dir,'dir')
    mkdir(kml_base_dir)
  end
  run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param);
  
  % Create concatenated and browse files for all data with thickness
  csv_base_dir = fullfile(post_path,'csv_good');
  kml_base_dir = fullfile(post_path,'kml_good');
  if ~exist(kml_base_dir,'dir')
    mkdir(kml_base_dir)
  end
  run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param);
end
