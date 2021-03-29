% script run_qlook
%
% Script for running qlook.m (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_qlook.m, qlook.m,
%   qlook_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Antarctica_Ground.xls'));

% Example to run specific segments and frames by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.qlook',0);
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200202_02');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200202_03');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200202_04');
params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200128_05');
params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200131_01');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200128_01');	params = ct_set_params(params,'cmd.frms',[11],'day_seg','20200128_01');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200128_06');	params = ct_set_params(params,'cmd.frms',[045],'day_seg','20200128_06');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200131_02');	params = ct_set_params(params,'cmd.frms',[016],'day_seg','20200131_02');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200201_01');	params = ct_set_params(params,'cmd.frms',[034],'day_seg','20200201_01');
%  params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200202_02'); params = ct_set_params(params,'cmd.frms',[079],'day_seg','20200202_01');
%  params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200208_03');	params = ct_set_params(params,'cmd.frms',[006],'day_seg','20200208_03');
% params=ct_set_params(params,'cmd.qlook',1,'day_seg','20200209_01');	params = ct_set_params(params,'cmd.frms',[050],'day_seg','20200209_01');

%Run specific segments and frames
% segs = {};frms={};
% segs{end+1} ='20200128_01'; frms{end+1} = 11;
% segs{end+1} ='20200128_01'; frms{end+1} = 13:15;
% segs{end+1} ='20200128_03'; frms{end+1} = 1;
% segs{end+1} ='20200128_04'; frms{end+1} = 1;
% segs{end+1} ='20200128_04'; frms{end+1} = 1;
% segs{end+1} ='20200128_05'; frms{end+1} = 1;
% segs{end+1} ='20200129_01'; frms{end+1} = 1;
% segs{end+1} ='20200129_02'; frms{end+1} = 1;
% segs{end+1} ='20200201_01'; frms{end+1} = 1;
% segs{end+1} ='20200209_01'; frms{end+1} = 43:47;
% segs{end+1} ='20200209_01'; frms{end+1} = 50;
% segs{end+1} ='20200201_02'; frms{end+1} = 1;
% segs{end+1} ='20200201_02'; frms{end+1} = 17:25;
% segs{end+1} ='20200209_02'; frms{end+1} = 1;
% 
% 
% %Mark segs for processing
% params = ct_set_params(params,'cmd.qlook',0);
% for pid = 1:length(params)
%   for sid = 1:length(segs)
%     if strcmp(segs{sid},params(pid).day_seg)
%       %Set generic to one
%       params = ct_set_params(params,'cmd.qlook',1,'day_seg',segs{sid});
%       %Update frms
%       old_frms = params(pid).cmd.frms;
%       params(pid).cmd.frms = [old_frms frms{sid}];
%     end
%   end
% end

% % Run all segments 
% params = ct_set_params(params,'cmd.qlook',1);
% params = ct_set_params(params,'cmd.qlook',0,'day_seg','20191211');
% params = ct_set_params(params,'cmd.qlook',0,'day_seg','20200116');
% % params = ct_set_params(params,'cmd.qlook',0,'day_seg','20200128_01');
% % params = ct_set_params(params,'cmd.qlook',0,'day_seg','20200202_02');
% % params = ct_set_params(params,'cmd.qlook',0,'day_seg','20200202_03');
% % params = ct_set_params(params,'cmd.qlook',0,'day_seg','20200202_04');


% QLOOK initial settings
params = ct_set_params(params,'radar.wfs(1).coh_noise_method','');
params = ct_set_params(params,'radar.wfs(1).deconv.en',false);
params = ct_set_params(params,'qlook.out_path','qlook');

%Turn on motion comp and resampling
params = ct_set_params(params,'qlook.resample', [2 1; 1 1]);
params = ct_set_params(params,'qlook.motion_comp',true);
% params = ct_set_params(params,'qlook.out_path','qlook_comp');

%Turn on bad_value
params = ct_set_params(params,'radar.wfs.bad_value',nan);

% QLOOK_NOISE
% params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
% params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis');
% params = ct_set_params(params,'radar.wfs(1).deconv.en',false);
% params = ct_set_params(params,'qlook.out_path','qlook_noise');
 
% QLOOK_DECONV
params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis');
params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
params = ct_set_params(params,'qlook.out_path','qlook_deconv');


% dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
param_override.cluster.cpu_time_mult  = 3;
param_override.cluster.mem_mult  = 2;
param_override.cluster.max_mem_mode = 'auto';

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
  if param.cmd.qlook
    ctrl_chain{end+1} = qlook(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
