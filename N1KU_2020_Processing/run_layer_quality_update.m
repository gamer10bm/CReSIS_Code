% % script run_layer_quality_update
%
% Updates the quality value for the specified layer 
%

%% User Settings
% ----------------------------------------------------------------------
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

% Example to run specific segments and frames by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_01');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_03');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_04');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_05');

params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_02');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_03');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_04');

quality_now = 1;
quality_set = 2;
layer_set = 'bottom';

%% Automated Section

% Process each of the segments
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
   %Load the layer files
   laydir = ct_filename_out(param,'','layer');
   laydata_fns = dir(fullfile(laydir,'Data_**.mat'));
   layfrmt_fn = dir(fullfile(laydir,'layer**.mat'));
   layfrmt = load(fullfile(laydir,layfrmt_fn.name));
   %Get the desired layer_id
   layid = layfrmt.lyr_id(strcmp(layer_set,layfrmt.lyr_name));
   %Set the quality for each layer file and then resave the file
   for f_id = 1:length(laydata_fns)
     laypath = fullfile(laydir,laydata_fns(f_id).name);
     %Load the file
     laynow = load(laypath);
     %Set quality_now values to quality_set values
     laynow.quality(layid,laynow.quality(layid,:)==quality_now) = quality_set;
     quality = laynow.quality;
     %Save the file
     save(laypath,'quality','-append');
   end
  end
end
