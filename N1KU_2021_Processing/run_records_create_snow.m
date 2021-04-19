% script run_create_records
%
% Script for running run_create_records (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records.m, create_records.m,
%   run_create_records_sync.m, check_records.m

%% User Setup
% =====================================================================
param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_GV.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2020_Arctic_Polar6.xls'));
params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.records',0);
% params = ct_set_params(params,'cmd.records',1,'day_seg','20170122_01');

params = ct_set_params(params,'cmd.records',0);
params = ct_set_params(params,'cmd.records',1,'day_seg','20210328');
% params = ct_set_params(params,'cmd.records',1,'day_seg','20210205');

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
