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

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

% Example to run specific segments and frames by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.sar',0);
daysegs = {'20200202_05'}; frms = {[8:12]};
for did = 1:length(daysegs)
  params = ct_set_params(params,'cmd.sar',1,'day_seg',daysegs{did});
  params = ct_set_params(params,'cmd.frms',frms{did},'day_seg',daysegs{did});
end
%Adjusted SAR settings
params = ct_set_params(params,'sar.sigma_x',0.5);
params = ct_set_params(params,'sar.out_path','sar2');
%Adjusted array settings
params = ct_set_params(params,'array.line_rng',-5:5);
params = ct_set_params(params,'array.dline',6);
params = ct_set_params(params,'array.in_path','sar2');
params = ct_set_params(params,'array.out_path','standard2');


%Set after deconv
% params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
% params = ct_set_params(params,'radar.wfs(2).deconv.en',true);

%Also remember to add for coherent noise removal:
params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis');
params = ct_set_params(params,'radar.wfs(1).deconv.en',false);

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
