clc
close all

%Load params
params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20190203_01');

for id_p = 1:length(params)
  qlook ={};
  param = params(id_p);
  
  if param.cmd.generic
    %% Load qlook
    qlook_dir = ct_filename_out(param,'qlook');
    qlook_fns = dir(qlook_dir);
    fn = qlook_fns(round(length(qlook_fns)/2)).name;
    
    qlook{end+1} = load(fullfile(qlook_dir,fn));
    
    %% Load qlook_noise_50
    qlook_dir = ct_filename_out(param,'qlook_noise');
    
    qlook{end+1} = load(fullfile(qlook_dir,fn));
    
%     %% Load qlook_noise_60
%     qlook_dir = ct_filename_out(param,'qlook_noise_60');
%     
%     qlook{end+1} = load(fullfile(qlook_dir,fn));
    %% Plot the two
    figure(22)
    ql = length(qlook);
    for qid = 1:ql
      subplot(1,ql,qid)
      imagesc(lp(qlook{qid}(1).Data))
    end
  end
end
