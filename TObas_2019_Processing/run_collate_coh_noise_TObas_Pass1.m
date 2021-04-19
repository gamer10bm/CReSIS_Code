% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];


params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'',{'analysis_noise','analysis'});
% 
params = ct_set_params(params,'cmd.generic',0);
daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
  '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
  '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
for did = [2,3,6:length(daysegs)]%3:length(daysegs)
  params = ct_set_params(params,'cmd.generic',1,'day_seg',daysegs{did});
end
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20190129_01');

if 1
  param_override.collate_coh_noise.method = 'firdec';
  param_override.collate_coh_noise.firdec_fs = 1/7.5;
  param_override.collate_coh_noise.firdec_fcutoff = @(t) 1/30;
else
  param_override.collate_coh_noise.method = 'dft';
  param_override.collate_coh_noise.dft_corr_time = inf;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  for wf = 1:length(params(param_idx).radar.wfs)
    params(param_idx).collate_coh_noise.threshold_eval{wf} = 'threshold = max(min(-100,threshold + 20),10*log10(abs(noise.dft(:,1)).^2)+6);';
  end
end

param_override.collate_coh_noise.debug_plots = {'cn_plot','threshold_plot'};
% param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot'};


% ENABLE NEXT TWO LINES ON SECOND PASS
% param_override.collate_coh_noise.in_path = 'analysis_threshold';
% param_override.collate_coh_noise.out_path = 'analysis_threshold';
%% Debug
% param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot'}; % Debugging
param_override.collate_coh_noise.debug_plots = {'cn_plot','threshold_plot'}; % Typical setting when not debugging
% param_override.collate_coh_noise.debug_plots = {}; % Necessary if plots are too large for memory

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
%   collate_coh_noise(param,param_override);
  collate_coh_noise
  
end
