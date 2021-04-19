% script run_check_surface
%
% Runs check_surface.m
%
% cat /N/dcwan/projects/cresis/output/ct_tmp/check_surface/snow/2017_Greenland_P3/*.txt
%
% Author: John Paden

%% User Settings
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
out=1;
% params = ct_set_params(params,'radar.wfs(1).nz_valid',[0 1]);


% params = ct_set_params(params,'radar.nz_valid',[0 1 2 3]);

% params.cmd.generic=1;

params = ct_set_params(params,'cmd.generic',0);
daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
  '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
  '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
for did = 1:length(daysegs)
  params = ct_set_params(params,'cmd.generic',1,'day_seg',daysegs{did});
end
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20190129_02');

% out = 0;

% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','^sea.*');
% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','(?(?!^sea.*)^.*)');
% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

% param_override.check_surface.debug_plots = {'visible','twtt','gps','nz'};
param_override.check_surface.debug_plots = {'NA'};
param_override.check_surface.save_records_en = false;

%Set offsets
% param_override.check_surface.radar_twtt_offset = -169e-9;
% param_override.check_surface.radar_gps_time_offset = -18;

%Set tdc_adjust to 0 for output purposes
params = ct_set_params(params,'radar.wfs(1).Tadc_adjust',0);
params = ct_set_params(params,'radar.wfs(2).Tadc_adjust',0);

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
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    check_surface(param,param_override);
  end
end

%Output all results as table
if ~ispc && out
  checksur_dir = fileparts(ct_filename_ct_tmp(params(1),'','check_surface','*.txt'));
  system(sprintf('cat %s/*.txt',checksur_dir))
end