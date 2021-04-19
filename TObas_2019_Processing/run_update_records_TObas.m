% script run_update_records
%
% This script setsup the parameters and calls update_records.  Make
% a local copy of the file in your personal folder.
%
% Author: John Paden
%
% See also update_records.m

%% User Settings
% =========================================================================
param_override = [];

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));

%% Automated section
% =========================================================================

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
  
  update_records(param,param_override);
end
