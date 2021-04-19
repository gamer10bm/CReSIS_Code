% Script run_post
%
% Loads the "post" worksheet from the parameter spreadsheet and then calls
% post with this information.
%
% Authors: Theresa Stumpf, John Paden, Reece Mathews
%
% See also: post.m

%% User Settings
param_override = [];

if 0
  params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
  
  % Example to run a specific segment and frame by overriding parameter spreadsheet values
  % params = ct_set_params(params,'cmd.generic',0);
  % params = ct_set_params(params,'cmd.generic',1,'day_seg','20190201_01');
  % params = ct_set_params(params,'cmd.frms',[]);
  % params = ct_set_params(params,'post.data_dirs',{'standard'});
  % params = ct_set_params(params,'post.out_path','post_standard');
  % params = ct_set_params(params,'post.data_dirs',{'qlook'});
  % params = ct_set_params(params,'post.out_path','post');
elseif 1
  %   params = read_param_xls(ct_filename_param('snow_param_2020_Arctic_Polar6.xls'));
  params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'));
  
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20210302');
  
  
%   params = ct_set_params(params,'post.img',0);
%   params = ct_set_params(params,'post.data_dirs',{'qlook'});
%   params = ct_set_params(params,'post.out_path','post_img1');
  params = ct_set_params(params,'post.img',0);
  params = ct_set_params(params,'post.data_dirs',{'qlook_noise_deconv'});
  params = ct_set_params(params,'post.out_path','post_img1');
%   params = ct_set_params(params,'post.img',2);
%   params = ct_set_params(params,'post.data_dirs',{'qlook_noise_deconv'});
%   params = ct_set_params(params,'post.out_path','post_img2');
%   params = ct_set_params(params,'post.img',3);
%   params = ct_set_params(params,'post.data_dirs',{'qlook_noise_deconv'});
%   params = ct_set_params(params,'post.out_path','post_img3');
%   params = ct_set_params(params,'post.img',4);
%   params = ct_set_params(params,'post.data_dirs',{'qlook_noise_deconv'});
%   params = ct_set_params(params,'post.out_path','post_img4');
  %   params = ct_set_params(params,'post.img',1);
  %   params = ct_set_params(params,'post.out_path','post_VV');
  %   params = ct_set_params(params,'post.img',2);
  %   params = ct_set_params(params,'post.out_path','post_HV');
  %   params = ct_set_params(params,'post.img',3);
  %   params = ct_set_params(params,'post.out_path','post_HH');
  %   params = ct_set_params(params,'post.img',4);
  %   params = ct_set_params(params,'post.out_path','post_VH');
  
  %John settings
  % Input images
  param_override.echogram_to_jpeg.data_type = 'qlook_noise_deconv';

  % Input top/bottom layers
  param_override.echogram_to_jpeg.layers = [struct('name', 'surface', 'source', 'layerdata', 'layerdata_source', 'layer', 'existence_check', false) ...
    struct('name', 'bottom', 'source', 'layerdata', 'layerdata_source', 'layer', 'existence_check', false)];

  % Output path
  param_override.echogram_to_jpeg.mat_out_path = 'CSARP_post/small_mat';
  param_override.echogram_to_jpeg.jpeg_out_path = 'CSARP_post/small_jpg';

  % Truncation parameters
  param_override.echogram_to_jpeg.N_before_surface = 500;
  param_override.echogram_to_jpeg.N_after_surface = 1500;
  param_override.echogram_to_jpeg.N_after_bottom = 100;

  % Decimation parameters
  param_override.echogram_to_jpeg.decimate_fasttime = 1;
  param_override.echogram_to_jpeg.decimate_slowtime = 1;

  
  params = ct_set_params(params,'post.maps_en',0);
  params = ct_set_params(params,'post.echo_en',1);
  params = ct_set_params(params,'post.layers_en',0);
  param_override.post.layers.existence_check = false;
  param_override.post.layers = [struct('name', 'surface', 'source', 'layerData')];
  params = ct_set_params(params,'post.echo_with_no_layer_en',false);
  params = ct_set_params(params,'post.data_en',0);
  params = ct_set_params(params,'post.csv_en',0);
  params = ct_set_params(params,'post.concat_en',0);
  params = ct_set_params(params,'post.pdf_en',1);
  
end

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
  post(param,param_override);
end

%% Create season-wide concatenated and browse files (CSV and KML)
% =====================================================================
concatenate_csv_kml = false;
for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if ct_generic_en(param)
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
