% script run_make_layer_files
%
% Calls the make_layer_files function
%
% Author: John Paden
%
% See also: make_layer_files

%% User Settings
% =========================================================================
param_override = [];

% Parameters spreadsheet to use for updating
params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'));
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20190130_01');
% params = ct_set_params(params,'cmd.frms',[]);

% .echogram_input = ct_filename_out path argument for which 
%   radar echograms to use for grabbing the initial surface values and the
%   time axis from. Typical values are shown here.
layer.echogram_input = 'qlook';
% layer.echogram_input = 'standard';
% layer.echogram_input = 'CSARP_post/qlook';
% layer.echogram_input = 'CSARP_post/standard';

% 0: Data_YYYYMMDD_SS_FFF.mat file, I: Data_II_YYYYMMDD_SS_FFF.mat file
layer.echogram_img = 1;

layer.layer_output = 'layerData';
% layer.layer_output = 'CSARP_post/layerData';

layer.do_not_overwrite_layer_files = true;

%% Automated section
% =========================================================================

param_override.make_layer_files = layer;

global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  make_layer_files(param,param_override);
end
