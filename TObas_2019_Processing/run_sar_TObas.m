% script run_sar
%
% Script for running sar.m (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'');

% Example to run specific segments and frames by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.sar',0);
daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
  '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
  '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
for did = [1:3,5:length(daysegs)]
  params = ct_set_params(params,'cmd.sar',1,'day_seg',daysegs{did});
end
% params = ct_set_params(params,'cmd.sar',1,'day_seg','20190131_01');
% params = ct_set_params(params,'cmd.frms',[]);

%Set after deconv
params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
params = ct_set_params(params,'radar.wfs(2).deconv.en',true);

%Also remember to add for coherent noise removal:
params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
params = ct_set_params(params,'radar.wfs(2).coh_noise_method','analysis');
params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis_threshold');
params = ct_set_params(params,'radar.wfs(2).coh_noise_arg.fn','analysis_threshold');

% dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
param_override.cluster.cpu_time_mult  = 10;
% param_override.cluster.mem_mult  = 3;

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
  if param.cmd.sar
    ctrl_chain{end+1} = sar(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
