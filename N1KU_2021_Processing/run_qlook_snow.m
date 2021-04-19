% script run_qlook
%
% Script for running qlook (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_qlook.m, qlook.m,
%   qlook_task.m

%% User Setup
% =====================================================================
param_override = [];

if 0
  params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
  % params = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'));
  
  % Example to run a specific segment and frame by overriding parameter spreadsheet values
  params = ct_set_params(params,'cmd.qlook',0);
  params = ct_set_params(params,'cmd.qlook',1,'day_seg','20120330_04');
  params = ct_set_params(params,'cmd.frms',[194 195]);
  
  % params = ct_set_params(params,'qlook.out_path','qlook');
  % params = ct_set_params(params,'qlook.out_path','qlook_noise');
  % params = ct_set_params(params,'qlook.out_path','qlook_deconv');
  
  % params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
  
  params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
  
elseif 0
  params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_GV.xls'),'',{'analysis_noise' 'analysis'});
  
  params = ct_set_params(params,'cmd.qlook',0);
%   params = ct_set_params(params,'cmd.qlook',1,'day_seg','20190914_01');
%   params = ct_set_params(params,'cmd.frms',[18]);
%   params = ct_set_params(params,'cmd.frms',[364]);
%   params = ct_set_params(params,'cmd.frms',[200]);
%   params = ct_set_params(params,'cmd.frms',[18 19 277 364]);
%   params = ct_set_params(params,'cmd.frms',[19]);
  params = ct_set_params(params,'cmd.qlook',1,'day_seg','20190909_02');
  params = ct_set_params(params,'cmd.frms',[14]);
%   params = ct_set_params(params,'cmd.frms',[16 17 18]);
%  params = ct_set_params(params,'cmd.frms',[12 14 29]);
 params = ct_set_params(params,'cmd.frms',[12 14 29 39]);
  
  %params = ct_set_params(params,'radar.wfs.BW_window',[2.4e9 8e9]);
%   params = ct_set_params(params,'radar.wfs.BW_window',[2.4e9 17e9]);

    params = ct_set_params(params,'qlook.surf.en',false);
    params = ct_set_params(params,'qlook.surf_layer.layerdata_source','paden/layer');

  params = ct_set_params(params,'radar.wfs.coh_noise_method','analysis');
  
  params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis');
  params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_threshold'); % TEMP DEBUG
  params = ct_set_params(params,'radar.wfs(1).BW_window',[2800000000 17480000000]);
  
  if 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',false);
    params = ct_set_params(params,'qlook.out_path','paden/qlook');
  elseif 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_deconv');
    params = ct_set_params(params,'qlook.resample',[209 367]);
    params = ct_set_params(params,'qlook.out_path','paden/qlook_deconv');
  elseif 1
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_deconv3');
    params = ct_set_params(params,'qlook.resample',[280 367]);
    params = ct_set_params(params,'qlook.out_path','paden/qlook_deconv3');
  elseif 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'qlook.resample',[710 367]);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_uwb');
    params = ct_set_params(params,'qlook.out_path','paden/qlook_uwb');
  elseif 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'qlook.resample',[250 367]);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_kuband');
    params = ct_set_params(params,'qlook.out_path','paden/qlook_kuband');
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_kuband_frm12');
    params = ct_set_params(params,'qlook.out_path','paden/qlook_kuband_frm12');
%     params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_kuband_frm29');
%     params = ct_set_params(params,'qlook.out_path','paden/qlook_kuband_frm29');
  elseif 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',false);
    params = ct_set_params(params,'qlook.out_path','paden/qlook_narrow');
    params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow');
    params = ct_set_params(params,'radar.wfs(1).BW_window',[2800000000 7480000000]);
  elseif 1
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_narrow2');
    params = ct_set_params(params,'qlook.resample',[35 19]);
    params = ct_set_params(params,'qlook.out_path','paden/qlook_narrow2_deconv');
    params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow2_threshold');
    params = ct_set_params(params,'radar.wfs(1).BW_window',[2400000000 8480000000]);
  elseif 0
    params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
    params = ct_set_params(params,'radar.wfs(1).deconv.fn','paden/analysis_narrow');
    params = ct_set_params(params,'qlook.resample',[209 117]);
    params = ct_set_params(params,'qlook.out_path','paden/qlook_narrow_deconv');
    params = ct_set_params(params,'radar.wfs.coh_noise_arg.fn','paden/analysis_narrow');
    params = ct_set_params(params,'radar.wfs(1).BW_window',[2800000000 7480000000]);
  end
  
elseif 0
  params = read_param_xls(ct_filename_param('snow_param_2020_Arctic_Polar6.xls'));
  
  params = ct_set_params(params,'cmd.qlook',0);
  params = ct_set_params(params,'cmd.qlook',1,'day_seg','20200206');
  
elseif 1
  params = read_param_xls(ct_filename_param('snow_param_2020_SouthDakota_N1KU.xls'));
  params = ct_set_params(params,'cmd.qlook',0);
  %params = ct_set_params(params,'cmd.qlook',1,'day_seg','20210204');
%   params = ct_set_params(params,'cmd.qlook',1,'day_seg','20210219');
%   params = ct_set_params(params,'cmd.qlook',1,'day_seg','20210225_03');
%   params = ct_set_params(params,'cmd.qlook',1,'day_seg','20210328');
  params = ct_set_params(params,'cmd.qlook',1,'day_seg','20210219_01');
  params = ct_set_params(params,'cmd.frms',[41:68],'day_seg','20210219_01');
%   params = ct_set_params(params,'cmd.frms',[32 33]);

%   params = ct_set_params(params,'qlook.out_path','qlook_noise');
  params = ct_set_params(params,'radar.wfs(1).coh_noise_method','analysis');
  params = ct_set_params(params,'radar.wfs(1).coh_noise_arg.fn','analysis_adjust');
  
  params = ct_set_params(params,'qlook.out_path','qlook_noise_adjust_deconv');
  params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
  params = ct_set_params(params,'radar.wfs(1).deconv.fn','analysis');
  params = ct_set_params(params,'qlook.resample',[2 1]);
  
  params = ct_set_params(params,'qlook.motion_comp',true);
  
%   params = ct_set_params(params,'qlook.imgs',{[1 1],[1 2],[1 3],[1 4]});
%   params = ct_set_params(params,'qlook.img_combine',[]);
end

% dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
param_override.cluster.rerun_only = false;
% param_override.cluster.desired_time_per_job  = 240*60;
param_override.cluster.cpu_time_mult  = 2;
param_override.cluster.mem_mult  = 2;
param_override.cluster.max_mem_mode = 'truncate';

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

ctrl_chain = cluster_run(ctrl_chain);