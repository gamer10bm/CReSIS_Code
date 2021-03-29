% script run_gopro.m

% The purpose of this script is to run gopro preprocessing for all photos,
  % processing of photos by segment, and post processing for filling in
  % missing photos

%% User Settings
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20200128_01');


params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {}; proc_output = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    %Process gopro pictures by segments
    proc_output{end+1} = gopro_process(param,param_override);
  end
end

%Post process 
%Save all paths and gps locations
gopro_postprocess(proc_output)