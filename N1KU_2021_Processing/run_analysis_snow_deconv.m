% script run_analysis_snow_deconv
%
% Script for running analysis
% https://ops.cresis.ku.edu/wiki/index.php/Analysis
%
% Authors: John Paden
%
% See also: master.m, run_analysis.m, analysis.m, analysis_task.m

%% User Setup
% =====================================================================
param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_GV.xls'),'',{'analysis_spec' 'analysis'});
params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'),'',{'analysis_spec' 'analysis'});

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20210219_01');

params = ct_set_params(params,'analysis.imgs',{[1 1; 1 2; 1 3; 1 4]});
params = ct_set_params(params,'analysis.imgs',{[1 2; 1 3; 1 4]});

% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow');
% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow_threshold'); % TEMP DEBUG
% params = ct_set_params(params,'radar.wfs(1).BW_window',[2800000000 7480000000]);
% params = ct_set_params(params,'analysis.out_path','paden/analysis_narrow');

% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow2');
% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow2_threshold'); % TEMP DEBUG
% params = ct_set_params(params,'radar.wfs(1).BW_window',[2400000000 8480000000]);
% params = ct_set_params(params,'analysis.out_path','paden/analysis_narrow2');

% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis');
% params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_threshold'); % TEMP DEBUG
% params = ct_set_params(params,'radar.wfs(1).BW_window',[2800000000 17480000000]);
% params = ct_set_params(params,'analysis.out_path','paden/analysis');
  
  params = ct_set_params(params,'qlook.out_path','qlook_noise_threshold');
  params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
  params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis');

% param.radar.wfs.coh_noise_method = '';

% ENABLE ON SECOND PASS OF COHERENT NOISE
% params = ct_set_params(params,'analysis.cmd{1}.threshold','');
params = ct_set_params(params,'analysis.out_path','analysis_adjust');

% ENABLE AFTER TWO COHERENT NOISE PASSES ARE COMPLETED
% params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
% params = ct_set_params(params,'radar.wfs(2).coh_noise_method','analysis');
% params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis_threshold');
% params = ct_set_params(params,'radar.wfs(2).coh_noise_arg.fn','analysis_threshold');

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
param_override.cluster.cpu_time_mult  = 4;
% param_override.cluster.mem_mult  = 2;

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
    ctrl_chain{end+1} = analysis(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

