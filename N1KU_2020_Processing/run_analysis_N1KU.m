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

if 1
  %Load analysis_spec
  params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'),'',{'analysis_spec','analysis'});
  %Adjust peakiness value
  params = ct_set_params(params,'analysis.cmd{1}.threshold',32);
  %Adjust rlines bin size
  rlines_set = 32; %Must be whole number when divided by 4 and greater than 51
  params = ct_set_params(params,'analysis.cmd{1}.rlines',rlines_set); %Default is 128
%   %Adjust signal_doppler_bins based on rline
%   sigdif = ceil(3/128*rlines_set); %default is 3
%   sigdop_bins_set = [1:1+sigdif rlines_set-sigdif:rlines_set]; %Default is [1:4 125:128]
%   params = ct_set_params(params,'analysis.cmd{1}.signal_doppler_bins',sigdop_bins_set);
%   %Adjust noise_doppler_bins based on rline
%   noise_strt = ceil(12/128*rlines_set); noise_end = floor(117/128*rlines_set);
%   noisedop_bins_set = [noise_strt:noise_end]; %default [12:117];
%   params = ct_set_params(params,'analysis.cmd{1}.noise_doppler_bins',noisedop_bins_set);
  %Enable for analysis_spec
  params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
  params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis');
elseif 0
  %Load analysis_noise
  params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'),'',{'analysis_noise','analysis'});
end

% Example to run specific segments and frames by overriding parameter spreadsheet values
% segs = {}; 
% segs{end+1} ='20200128_01'; 
% segs{end+1} ='20200128_01'; 
% segs{end+1} ='20200128_03';
% segs{end+1} ='20200128_04';
% segs{end+1} ='20200128_04';
% segs{end+1} ='20200128_05';
% segs{end+1} ='20200129_01';
% segs{end+1} ='20200129_02';
% segs{end+1} ='20200201_01'; 
% segs{end+1} ='20200209_01';
% segs{end+1} ='20200209_01';
% segs{end+1} ='20200201_02';
% segs{end+1} ='20200201_02';
% segs{end+1} ='20200209_02';
% 
% 
% %Mark segs for processing
% params = ct_set_params(params,'cmd.generic',0);
% for pid = 1:length(params)
%   for sid = 1:length(segs)
%     if strcmp(segs{sid},params(pid).day_seg)
%       %Set generic to one
%       params = ct_set_params(params,'cmd.generic',1,'day_seg',segs{sid});
%     end
%   end
% end

params = ct_set_params(params,'cmd.generic',0);
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_01');
params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_05');
% % params=ct_set_params(params,'cmd.generic',1,'day_seg','20200129_01');
params=ct_set_params(params,'cmd.generic',1,'day_seg','20200131_01');
% % params=ct_set_params(params,'cmd.generic',1,'day_seg','20200208_03');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200201_02');
% % params=ct_set_params(params,'cmd.generic',1,'day_seg','20200209_01');
% % params=ct_set_params(params,'cmd.generic',1,'day_seg','20200202_01');


% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_02');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_03');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_04');
% % params = ct_set_params(params,'cmd.generic',0,'day_seg','20200128_01');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200128_03');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200128_06');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200131_02');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200201_01');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200202_01');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200208_03');
% params=ct_set_params(params,'cmd.generic',0,'day_seg','20200209_01');

% dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;
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
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    ctrl_chain{end+1} = analysis(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

