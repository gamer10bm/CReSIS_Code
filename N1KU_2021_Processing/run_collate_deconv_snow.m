% Script run_collate_deconv_snow.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_GV.xls'),'',{'analysis_spec' 'analysis'});
params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'),'',{'analysis_spec' 'analysis'});

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20210219_01');
% params = ct_set_params(params,'cmd.generic',1,'cmd.notes','deconv');
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');

% param_override.collate_deconv.debug_plots = {'peakiness','rbins','deconv','metric','final','visible'};
param_override.collate_deconv.debug_plots = {'metric','final','visible','rbins_best','deconv_best'};
% param_override.collate_deconv.debug_plots = {'final','visible'};

% param_override.collate_deconv.debug_rlines = [2 44];

param_override.collate_deconv.imgs = 1;
param_override.collate_deconv.rbins{1} = [-400 375];

param_override.collate_deconv.surf_layer.layerdata_source = 'layer';
% param_override.collate_deconv.surf_layer.layerdata_source = 'paden/layer';

if 1
  params = ct_set_params(params,'analysis.imgs',{[1 1; 1 2; 1 3; 1 4]});
  params = ct_set_params(params,'analysis.imgs',{[1 2; 1 3; 1 4]});
  params = ct_set_params(params,'collate_deconv.imgs',1);
  
  params = ct_set_params(params,'collate_deconv.f0',2560000000+145e6);
  params = ct_set_params(params,'collate_deconv.f1',6208000000);
%   param_override.collate_deconv.in_path = 'analysis';
%   param_override.collate_deconv.out_path = 'analysis';
  param_override.collate_deconv.in_path = 'analysis_adjust';
  param_override.collate_deconv.out_path = 'analysis_adjust';
  params = ct_set_params(params,'collate_deconv.SL_guard_bins',9);
  params = ct_set_params(params,'collate_deconv.abs_metric',[90 9.5 -34 -34 inf inf]);
  
%   param_override.collate_deconv.decimate_table_sec = 1;
  %params = ct_set_params(params,'collate_deconv.gps_times',[1.568037672174584e+09 -inf 1.568037826300003e+09; 1.568038282413137e+09 1.568037826300003e+09 inf],'day_seg','20190909_02');
  %params = ct_set_params(params,'collate_deconv.gps_times',[1568037672.05 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
  %params = ct_set_params(params,'collate_deconv.gps_times',[1568037670.24 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
%   params = ct_set_params(params,'collate_deconv.gps_times',[1568037719.61 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
  
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',2800000000);
  params = ct_set_params(params,'collate_deconv.f1',17000000000);
  param_override.collate_deconv.in_path = 'paden/analysis';
  param_override.collate_deconv.out_path = 'paden/analysis_uwb';
  %   params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);
  %   params = ct_set_params(params,'collate_deconv.abs_metric',[90 3.8 -34 -34 inf inf]);
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',2800000000);
  params = ct_set_params(params,'collate_deconv.f1',7480000000-500e6);
  param_override.collate_deconv.in_path = 'paden/analysis';
  param_override.collate_deconv.out_path = 'paden/analysis_deconv';
    params = ct_set_params(params,'collate_deconv.SL_guard_bins',11);
    params = ct_set_params(params,'collate_deconv.abs_metric',[90 9.5 -34 -34 inf inf]);
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',2800000000);
  params = ct_set_params(params,'collate_deconv.f1',8400000000);
  param_override.collate_deconv.in_path = 'paden/analysis';
  param_override.collate_deconv.out_path = 'paden/analysis_deconv3';
    params = ct_set_params(params,'collate_deconv.SL_guard_bins',9);
    params = ct_set_params(params,'collate_deconv.abs_metric',[90 9.5 -34 -34 inf inf]);
    
  param_override.collate_deconv.decimate_table_sec = 1;
  %params = ct_set_params(params,'collate_deconv.gps_times',[1.568037672174584e+09 -inf 1.568037826300003e+09; 1.568038282413137e+09 1.568037826300003e+09 inf],'day_seg','20190909_02');
  %params = ct_set_params(params,'collate_deconv.gps_times',[1568037672.05 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
  %params = ct_set_params(params,'collate_deconv.gps_times',[1568037670.24 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
  params = ct_set_params(params,'collate_deconv.gps_times',[1568037719.61 -inf 1.568037826300003e+09; 1568038282.54 1.568037826300003e+09 inf],'day_seg','20190909_02');
  
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',12000000000);
  params = ct_set_params(params,'collate_deconv.f1',17000000000);
  param_override.collate_deconv.in_path = 'paden/analysis';
  param_override.collate_deconv.out_path = 'paden/analysis_kuband';
  param_override.collate_deconv.out_path = 'paden/analysis_kuband_frm12';
  param_override.collate_deconv.out_path = 'paden/analysis_kuband_frm29';
  params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);
  %param_override.collate_deconv.rbins{1} = [-30 30];
  % 
  param_override.collate_deconv.decimate_table_sec = 1;
  param_override.collate_deconv.gps_times = [1.568037672497153e+09 -inf inf]; % Use frame 12
  param_override.collate_deconv.gps_times = [1.568038279510004e+09 -inf inf]; % Use frame 29
  params = ct_set_params(params,'collate_deconv.gps_times',[1.568037672497153e+09 -inf 1.568037826300003e+09; 1.568038279510004e+09 1.568037826300003e+09 inf],'day_seg','20190909_02');
%   param_override.collate_deconv.gps_times = [1.568037672497153e+09 -inf 1.568037826300003e+09; 1.568038279510004e+09 1.568037826300003e+09 inf];
%   param_override.collate_deconv.gps_times = [1.568037672497153e+09 -inf 1.568037865300003e+09; 1.568038121708597e+09 1.568037865300003e+09 inf];
%   param_override.collate_deconv.gps_times = [1.568037672497153e+09 -inf 1.568037826310001e+09; 1.568038121708597e+09 1.568037826310001e+09 inf];
  %param_override.collate_deconv.gps_times = [1.568037671658473e9 -inf 1.568037826310001e+09; 1.568038279445490e+09 1.568037826310001e+09 inf];
  params = ct_set_params(params,'collate_deconv.abs_metric',[65 8.9 -30 -30 inf inf]);
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
  %param_override.collate_deconv.surf_layer.layerdata_source = 'CSARP_post/layer';
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',2400000000);
  params = ct_set_params(params,'collate_deconv.f1',8000000000);
  params = ct_set_params(params,'collate_deconv.gps_times',[1.568037672174584e+09 -inf 1.568037826300003e+09; 1.568038279574518e+09 1.568037826300003e+09 inf],'day_seg','20190909_02');
  param_override.collate_deconv.in_path = 'paden/analysis_narrow2';
  param_override.collate_deconv.out_path = 'paden/analysis_narrow2';
  param_override.collate_deconv.rbins{1} = [-287 80];
  % params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);
  params = ct_set_params(params,'collate_deconv.abs_metric',[61 3.44 -34 -34 inf inf]);
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
  
elseif 0
  params = ct_set_params(params,'collate_deconv.f0',2800000000);
  params = ct_set_params(params,'collate_deconv.f1',7480000000-500e6);
  param_override.collate_deconv.in_path = 'paden/analysis_narrow';
  param_override.collate_deconv.out_path = 'paden/analysis_narrow';
  param_override.collate_deconv.rbins{1} = [-287 80];
  % params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);
  params = ct_set_params(params,'collate_deconv.abs_metric',[61 3.44 -34 -34 inf inf]);
  % 759	3.80	515608	78.8	3.4	-36.9	-35.5	-24.5	-25.8		27.8	0
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
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  collate_deconv(param,param_override);
  %collate_deconv
  
end
