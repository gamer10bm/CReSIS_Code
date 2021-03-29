% script run_records_create
%
% Script for running run_records_create (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_records_create.m, records_create.m,
%   run_records_create_sync.m, records_check.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));
% params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.records',1);
params = ct_set_params(params,'cmd.records',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.records',0,'day_seg','20200116');
params = ct_set_params(params,'records.file.base_dir','/cresis/snfs1/data/SnowRadar/2020_SouthDakota_CESSNA/');

dbstop if error;

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
for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.records
    records_create(param,param_override);
  end
end
