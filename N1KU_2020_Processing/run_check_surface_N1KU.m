% script run_check_surface
%
% Runs check_surface.m
%
% cat /N/dcwan/projects/cresis/output/ct_tmp/check_surface/snow/2017_Greenland_P3/*.txt
%
% Author: John Paden

%% User Settings
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

% Example to run specific segments and frames by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_03'); 
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200209_02'); 
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200208_02'); 
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_06'); 
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200131_03'); 
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200129_02'); params = ct_set_params(params,'cmd.frms',[11],'day_seg','20200129_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200209_02'); params = ct_set_params(params,'cmd.frms',[35, 36],'day_seg','20200209_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200208_02'); params = ct_set_params(params,'cmd.frms',[27, 28],'day_seg','20200208_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_06'); params = ct_set_params(params,'cmd.frms',[52, 53],'day_seg','20200202_06');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200131_03'); params = ct_set_params(params,'cmd.frms',[12],'day_seg','20200131_03');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_02'); params = ct_set_params(params,'cmd.frms',[49, 50],'day_seg','20200128_02');


params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200128_01');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_02');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_03');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_04');


% %Check Surface settings
%Radar layer
param_override.check_surface.radar_layer_params.name = 'surface';
param_override.check_surface.radar_layer_params.source = 'layerdata';

%Ref Layer
params = ct_set_params(params,'check_surface.ref_layer_params.name','USGS_surface');
params = ct_set_params(params,'check_surface.ref_layer_params.source','layerdata');
params = ct_set_params(params,'check_surface.ref_layer_params.layerdata_source','layer');
params = ct_set_params(params,'check_surface.ref_layer_params.existence_check',false);

%Override for now
% params = ct_set_params(params,'post.ops.location','usa');

param_override.check_surface.debug_plots = {'twtt','gps','nz'};
% param_override.check_surface.debug_plots = {'visible','twtt','gps','nz'};
% param_override.check_surface.debug_plots = {'NA'};
% param_override.check_surface.save_records_en = true;

%Set offsets
% param_override.check_surface.radar_twtt_offset = -103e-9;
% param_override.check_surface.radar_gps_time_offset = -18;

%Set tdc_adjust to 0 for output purposes
% params = ct_set_params(params,'radar.wfs(1).t_ref',103e-9);

%Set nz_valid
params = ct_set_params(params,'radar.nz_valid',[0 1 2 3]);
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
out = 1;
if ~ispc && out
  checksur_dir = fileparts(ct_filename_ct_tmp(params(1),'','check_surface','*.txt'));
  system(sprintf('cat %s/*.txt',checksur_dir))
end