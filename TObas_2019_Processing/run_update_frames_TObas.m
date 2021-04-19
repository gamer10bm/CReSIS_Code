% script run_update_frames
%
% Script for running "update_frames.m"

clear param;
if 0 %Singular qlook product for 2018 Antartica TObas
  param.radar_name = 'accum';
  
  param.day_seg = '20190129_01';
  
  param.season_name = '2018_Antarctica_TObas';
  param.image_out_dir = {
    {ct_filename_out(rmfield(param,'day_seg'),'','post_noise/images',1)}};
  param.mat_out_dir = {{fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_post','CSARP_standard')}};
  param.mat_out_img = {0};
  
  param.img_type = 'echo.jpg';
  
  % RDS image debugging mode (for snow, kuband radars)
  fmcw_img_debug_mode = false;
  noise_time_buffer = 250e-9;
  noise_time_duration = 45e-9;
  img_sidelobe = -40;
  noise_threshold_offset_dB = 5.2;
  
elseif 1 %Comparison between qlook w/ coh_noise_removal and w/o removal
  param.radar_name = 'accum';
  
  param.day_seg = '20190129_02';
  
  param.season_name = '2018_Antarctica_TObas';
  param.image_out_dir = {
    {ct_filename_out(rmfield(param,'day_seg'),'','post_deconv/images',1)},... %w/ deconv
    {ct_filename_out(rmfield(param,'day_seg'),'','post_noise/images',1)}... %w/ coherent noise removal
    {ct_filename_out(rmfield(param,'day_seg'),'','post/images',1)}... %w/o coherent noise removal
    };
  param.mat_out_dir = {
    {fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_qlook_deconv')},... %w/ deconv
    {fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_qlook_noise')}... %w/ coherent noise removal
    {fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_qlook')}... %w/o coherent noise removal
    };
  param.mat_out_img = {0,0,0};
  
  param.img_type = 'echo.jpg';
  
  % RDS image debugging mode (for snow, kuband radars)
  fmcw_img_debug_mode = false;
  noise_time_buffer = 250e-9;
  noise_time_duration = 45e-9;
  img_sidelobe = -40;
  noise_threshold_offset_dB = 5.2;  
  
elseif 0 %Comparison between qlook w/ coh_noise_removal and w/o removal
  param.radar_name = 'accum';
  
  param.day_seg = '20190131_03';
  
  param.season_name = '2018_Antarctica_TObas';
  param.image_out_dir = {
    {ct_filename_out(rmfield(param,'day_seg'),'','post_noise/images',1)},... %w/ coherent noise removal
    {ct_filename_out(rmfield(param,'day_seg'),'','post/images',1)}... %w/o noise removal
    };
  param.mat_out_dir = {
    {fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_qlook_noise')},... %w/ coherent noise removal
    {fullfile(gRadar.out_path,param.radar_name,param.season_name,'CSARP_qlook')}... %w/o coherent noise removal
    };
  param.mat_out_img = {0,0};
  
  param.img_type = 'echo.jpg';
  
  % RDS image debugging mode (for snow, kuband radars)
  fmcw_img_debug_mode = false;
  noise_time_buffer = 250e-9;
  noise_time_duration = 45e-9;
  img_sidelobe = -40;
  noise_threshold_offset_dB = 5.2;  
elseif 0
  param.radar_name = 'snow';
  
  param.day_seg = '20170330_01';
  
  param.season_name = '2017_Greenland_P3';
  param.image_out_dir = {{ct_filename_out(rmfield(param,'day_seg'),'','post/images',1)},
    {ct_filename_out(rmfield(param,'day_seg'),'','post/images_uwb',1)},
    {ct_filename_out(rmfield(param,'day_seg'),'','post/images_kuband',1)},
    {ct_filename_out(rmfield(param,'day_seg'),'','post/images_deconv',1)}};
  param.mat_out_dir = {{'dummy_string'},
    {'/cresis/snfs1/dataproducts/ct_data/snow/2017_Greenland_P3/CSARP_post/CSARP_qlook_uwb'},
    {'/cresis/snfs1/dataproducts/ct_data/snow/2017_Greenland_P3/CSARP_post/CSARP_qlook_kuband'},
    {'/cresis/snfs1/dataproducts/ct_data/snow/2017_Greenland_P3/CSARP_post/CSARP_deconv'}};
  param.mat_out_img = {0,0,0,0};
  
  param.img_type = 'echo.jpg';
  
  % FMCW image debugging mode (for snow, kuband radars)
  fmcw_img_debug_mode = true;
  noise_time_buffer = 250e-9;
  noise_time_duration = 45e-9;
  img_sidelobe = -40;
  noise_threshold_offset_dB = 5.2;
end

% update_field: string containing the field in frames file to update.
% Currently there are two options:
%               ('proc_mode', 'quality')
if 1
  update_field = 'proc_mode'; update_field_type = 'double';
else
  update_field = 'quality'; update_field_type = 'mask';
  update_field_mask = {'turn all masks off','coherent noise','deconvolution artifact', ...
    'raised noise floor/vertical stripes','missing data','no good data','low SNR','unclassified','land or iceberg','sea ice'};
end

audio_tone_for_nonzero_nonisnan = true;
audio_tone_check_code = '~isnan(frames.(update_field)(frm)) && frames.(update_field)(frm) ~= 0';

%update_field_match = [1 2 3];
update_field_match = [];

update_frames;
