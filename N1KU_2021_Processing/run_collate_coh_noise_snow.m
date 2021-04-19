% script run_collate_coh_noise_snow
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_GV.xls'),'',{'analysis_noise' 'analysis'});
params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'),'',{'analysis_noise' 'analysis'});

params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20210205_01');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20210219_01');

params = ct_set_params(params,'analysis.imgs',{[1 1; 1 2; 1 3; 1 4]});
params = ct_set_params(params,'collate_coh_noise.imgs',1);

% param_override.collate_coh_noise.in_path = 'analysis';
% param_override.collate_coh_noise.out_path = 'analysis';
% param_override.collate_coh_noise.out_path = 'paden/analysis_threshold';

% param_override.collate_coh_noise.in_path = 'paden/analysis_narrow';
% param_override.collate_coh_noise.out_path = 'paden/analysis_narrow';
% param_override.collate_coh_noise.out_path = 'paden/analysis_narrow_threshold';

% param_override.collate_coh_noise.in_path = 'paden/analysis_narrow2';
% param_override.collate_coh_noise.out_path = 'paden/analysis_narrow2';
% param_override.collate_coh_noise.out_path = 'paden/analysis_narrow2_threshold';

% ENABLE ON SECOND PASS
param_override.collate_coh_noise.in_path = 'analysis_adjust';
param_override.collate_coh_noise.out_path = 'analysis_adjust';

% param_override.collate_coh_noise.in_path = 'analysis';
% param_override.collate_coh_noise.out_path = 'analysis';

param_override.collate_coh_noise.min_samples = 0.98;

if 0
  param_override.collate_coh_noise.method = {'firdec'};
  param_override.collate_coh_noise.firdec_fs = {1/640};
  param_override.collate_coh_noise.firdec_fcutoff = {@(t) 1/160};
elseif 1
  param_override.collate_coh_noise.method = {'firdec','firdec','firdec','firdec'};
  param_override.collate_coh_noise.firdec_fs = {1/160,1/160,1/160,1/160};
  param_override.collate_coh_noise.firdec_fcutoff = {@(t) 1/40,@(t) 1/40,@(t) 1/40,@(t) 1/40};
else
  param_override.collate_coh_noise.method = {'dft'};
  param_override.collate_coh_noise.dft_corr_time = inf;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
%   for wf = 1:length(params(param_idx).radar.wfs)
%     params(param_idx).collate_coh_noise.threshold_eval{wf} = 'threshold = max(min(-100,threshold + 20),10*log10(abs(dft_noise(:,1)).^2)+6);';
%   end
  for wf = 1:length(params(param_idx).radar.wfs)
    params(param_idx).collate_coh_noise.threshold_eval{wf} = 'threshold = max(-95,max_filt1(10*log10(abs(dft_noise(:,1)).^2)+6,15));';
  end
end

param_override.collate_coh_noise.debug_plots = {'cn_plot','threshold_plot'};
% param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','threshold_plot'};
% param_override.collate_coh_noise.debug_plots = {'cn_plot'};
% param_override.collate_coh_noise.debug_plots = {'visible','cn_plot'};
% param_override.collate_coh_noise.debug_plots = {'visible','cn_plot','reuse'};

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
%   collate_coh_noise(param,param_override);
  collate_coh_noise
  
end
