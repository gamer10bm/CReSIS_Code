% Script run_collate_deconv.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'),'',{'analysis_spec','analysis'});

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20200128_05');
params = ct_set_params(params,'cmd.generic',1,'day_seg','20200131_01');
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');
% 
% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_02');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_03');
% params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_04');


params = ct_set_params(params,'collate_deconv.f0',3e9);
params = ct_set_params(params,'collate_deconv.f1',5e9);
params = ct_set_params(params,'collate_deconv.rbins',{[-160 40]});
params = ct_set_params(params,'collate_deconv.abs_metric',[58 9.8 -25 -35 inf inf]);
params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);

% STEP 1: Check peakiness to ensure that enough waveforms qualify. If
% peakiness threshold has to be adjusted to let more waveforms in, then
% analysis spec must be run again.
param_override.collate_deconv.debug_plots = {'peakiness','metric'}; param_override.collate_deconv.stage_two_en = false;

% STEP 2: Use the "metric" table output to choose debug_rlines
%param_override.collate_deconv.debug_plots = {'metric','visible'}; param_override.collate_deconv.stage_two_en = false;

% STEP 3: To evaluate individual waveforms, set debug_rlines to these
% waveforms and enable rbins (evaluate SNR for the rbins setting you have
% chosen) and/or deconv (evaluate the SL_guard_bins, abs_metric, and
% sidelobe suppression achieved):
%param_override.collate_deconv.debug_rlines = [227];
%param_override.collate_deconv.debug_plots = {'deconv','metric','visible'}; param_override.collate_deconv.stage_two_en = false;
%param_override.collate_deconv.debug_plots = {'rbins','deconv','metric','visible'}; param_override.collate_deconv.stage_two_en = false;

% STEP 4: Once rbins are set and waveforms appear to be deconvolving well,
% run stage one and stage two (recommend disabling "visible" if many
% segments or wf_adc pairs).
%param_override.collate_deconv.debug_plots = {'metric','final','visible'};
%param_override.collate_deconv.debug_plots = {'metric','final'};

if 0
  % For debugging, use this to select specific images and wf_adc pairs to
  % collate instead of doing them all
  param_override.collate_deconv.wf_adcs = {[1],[1],[1]};
  param_override.collate_deconv.imgs = [1];
end
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
  collate_deconv(param,param_override);
%   collate_deconv
  
end
