clc
close all

%Load params
params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20200128_01');


qlook_suffs = {'noise','noise_threshold','noise_60'}; %


for id_p = 1:length(params)
  param = params(id_p);  
  if param.cmd.generic
    %% Plot the qlooks
    dat.fh = figure('numbertitle','off','name',fn,'keypressfcn',@(h,e)butt_press(h,e,@(id)qlook_check(param, qlook_suffs, id)));
    dat.fn_id = 1;
    dat.fn_id = qlook_check(param, qlook_suffs, 1);
    guidata(dat.fh, dat)
  end
end