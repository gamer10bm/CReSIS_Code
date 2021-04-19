% script run_analysis
%
% Script for running analysis.m
% https://ops.cresis.ku.edu/wiki/index.php/Analysis
%
% Authors: John Paden
%
% See also: master.m, run_analysis.m, analysis.m, analysis_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'',{'analysis_noise','analysis'});
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'),'','analysis');

% Example to run specific segments and frames by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',0);
daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
  '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
  '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
for did = [2,3,6:length(daysegs)]%3:length(daysegs)
  params = ct_set_params(params,'cmd.generic',1,'day_seg',daysegs{did});
end
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20190129_01');
% params = ct_set_params(params,'cmd.frms',[]);

% % ENABLE ON SECOND PASS
% params = ct_set_params(params,'analysis.cmd{1}.threshold','');
% params = ct_set_params(params,'analysis.out_path','analysis_threshold');

dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
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

