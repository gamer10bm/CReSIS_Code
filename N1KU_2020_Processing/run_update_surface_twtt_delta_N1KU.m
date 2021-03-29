% script run_update_surface_twtt_delta
%
% Runs update_surface_twtt_delta.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200128_02');
% params=ct_set_params(params,'cmd.generic',1,'day_seg','20200129_02'); 

params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'qlook_noise'});
% params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'qlook','deconv','qlook_uwb','qlook_kuband'}); % Snow Radar
params = ct_set_params(params,'update_surface_twtt_delta.imgs',[0]); % 0: combined img, 1+: img_II files
% params = ct_set_params(params,'update_surface_twtt_delta.update_adc_gains_dB',0);

% Example to override parameters
% for param_idx = 1:length(params)
%   %param = params(param_idx);
%   for wf = 1:length(params(param_idx).radar.wfs)
%     params(param_idx).radar.wfs(wf).Tsys = [0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9;
%     params(param_idx).radar.wfs(wf).chan_equal_dB = [6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2];
%     params(param_idx).radar.wfs(wf).chan_equal_deg = [-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6];
%   end
% end

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
  update_surface_twtt_delta(param,param_override);
end
