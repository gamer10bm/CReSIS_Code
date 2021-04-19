% Script run_collate_deconv.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'',{'analysis_spec' 'analysis'});

% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

params = ct_set_params(params,'cmd.generic',0);
% daysegs = {'20190129_01','20190129_02','20190130_01','20190131_01','20190131_02',...
%   '20190131_03','20190201_01','20190203_01','20190204_01','20190204_02','20190204_03',...
%   '20190205_01','20190206_01','20190206_02','20190207_01','20190207_02'};
% for did = [2,3,6:length(daysegs)]%3:length(daysegs)
%   params = ct_set_params(params,'cmd.generic',1,'day_seg',daysegs{did});
% end
params = ct_set_params(params,'cmd.generic',1,'day_seg','20190131_01');

params = ct_set_params(params,'collate_deconv.f0',615e6);
params = ct_set_params(params,'collate_deconv.f1',875e6);
params = ct_set_params(params,'collate_deconv.abs_metric',[58 5 -25 -35 inf inf]);
params = ct_set_params(params,'collate_deconv.SL_guard_bins',6);
params = ct_set_params(params,'collate_deconv.threshold',82);
params = ct_set_params(params,'collate_deconv.rbins',{[-200 650],[-200 650]});
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
