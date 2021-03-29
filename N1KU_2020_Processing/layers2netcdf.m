% % % script layers2netcdf
% clearvars -except gRadar
% close all
% clc

%% User Settings
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20200202_05');


params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_02');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_03');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200202_04');


%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Get valid layer directories
laydirs = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    laydirs{end+1} = ct_filename_out(param,'','layer');
  end
end

%Extract useful information from layer files
lay_fields = {'surface','USGS_surface','bottom'};
lay_names = {'Ground Range', 'USGS Surface Range','Snow Range'};
pull_fields = {'lat','lon','elev','gps_time'};
pull_names = {'Latitude (deg)', 'Longitude (deg)', 'Flight Elevation (MSL [m])',''};
layouts = {};
for lay_idx = 1:length(laydirs)
  %Find all .mat files in laydirs
  laydatafns = dir(fullfile(laydirs{lay_idx},'Data**.mat'));
  layfrmfn = dir(fullfile(laydirs{lay_idx},'layer**.mat'));
  %Load all laydatafns
  layeroutput = [];
  if ~isempty(layfrmfn) && ~isempty(laydatafns)
    %Load the formatting file (layer_**.mat)
    layfrmt = load(fullfile(laydirs{lay_idx},layfrmfn.name));
    for fn_id = 1:length(laydatafns)
      laynow = load(fullfile(laydirs{lay_idx},laydatafns(fn_id).name));
      %Grab the basic fields
      for fid = 1:length(pull_fields)
        outvec = []; laydat = laynow.(pull_fields{fid});      
        if isfield(layeroutput,pull_fields{fid})
          outvec = [layeroutput.(pull_fields{fid}) laydat];
        else
          outvec = laydat;
        end
        layeroutput.(pull_fields{fid}) = outvec;
      end
      %Grab each layerpull layer
      layout_fields = {}; layout_names = {};
      for tp_id = 1:length(lay_fields)
        %Determine row of the desired layer
        pulllogvec = layfrmt.lyr_id(strcmp(layfrmt.lyr_name,lay_fields{tp_id}));
        %Grab layer data
        if any(pulllogvec)
          %Grab twtt vector
          twttvec = laynow.twtt(pulllogvec,:);
          %Grab quality vectory
          qualvec = laynow.quality(pulllogvec,:);
        else
          %Make twtt nan and quality nan
          twttvec = nan(size(laynow.(pull_fields{1})));
          qualvec = nan(size(laynow.(pull_fields{1})));          
        end
        %Add twtt vec to layeroutput
        outvec = [];
        layf_twtt = sprintf('%s_twtt',lay_fields{tp_id});
        twttunits = '(s)';
        layn_twtt = sprintf('%s %s',lay_names{tp_id},twttunits);
        if isfield(layeroutput,layf_twtt)
          outvec = [layeroutput.(layf_twtt) twttvec];
        else
          outvec = twttvec;
        end
        layeroutput.(layf_twtt) = outvec;
        %Add quality vec to layeroutput
        outvec = [];
        layf_qual = sprintf('%s_qual',lay_fields{tp_id});
        layn_qual = sprintf('%s Quality',lay_names{tp_id});
        if isfield(layeroutput,layf_qual)
          outvec = [layeroutput.(layf_qual) qualvec];
        else
          outvec = qualvec;
        end
        layeroutput.(layf_qual) = outvec;
        %Add new fields and names for layout
        layout_fields = [layout_fields{:} {layf_twtt} {layf_qual}];
        layout_names = [layout_names{:} {layn_twtt} {layn_qual}];
      end
    end
  end
  %Send to cell array of outputs
  layouts{end+1} = layeroutput;
end
%Send all the layerdata to a single large structure
pull_fields = [pull_fields layout_fields];
pull_names = [pull_names layout_names];
layerdata_full = cell2struct(cell(size(pull_fields)),pull_fields,2);
for lid = 1:length(layouts)
  if ~isempty(layouts{lid})
    for fid = 1:length(pull_fields)
      layerdata_full.(pull_fields{fid}) = [layerdata_full.(pull_fields{fid}) layouts{lid}.(pull_fields{fid})];
    end
  end
end

%Write raw layer data to netcdf file
outfn = sprintf('%s_layers.nc',param.season_name);
if exist(outfn,'file')
  delete(outfn)
end
for fid = 1:length(pull_fields)
  varname = pull_fields{fid};
  if ~isempty(pull_names{fid})
    varname = pull_names{fid};
  end
  tmpvec = layerdata_full.(pull_fields{fid});
  nccreate(outfn,varname,'Dimensions',{'r',size(tmpvec,1),'c',size(tmpvec,2)});
  ncwrite(outfn,varname,tmpvec);
end

%% Math Section
physical_constants;
field_add = {}; name_add = {}; data_add = {};
% Add surface in meters
field_add{end+1} = 'surface_m'; name_add{end+1} = 'Flight Elevation (AGL [m])';
surface_range_m = layerdata_full.surface_twtt*c/2;
data_add{end+1} = surface_range_m;
% Add AGL Elevation
field_add{end+1} = 'h_agl'; name_add{end+1} = 'Ground Elevation (MSL [m])';
elev_msl_m = layerdata_full.elev;
data_add{end+1} = elev_msl_m-surface_range_m;
% Add Depth
field_add{end+1} = 'depth'; name_add{end+1} = 'Snow Depth (m)';
bottom_range_m = layerdata_full.bottom_twtt*c/2;
data_add{end+1} = surface_range_m- bottom_range_m;
% Add Depth Quality
field_add{end+1} = 'depth_qual'; name_add{end+1} = 'Snow Depth Quality';
  %Make all nan
  qualvec = nan(size(data_add{end}));
  valdepth = 2; %m %Based on statistical distribution
  %Make all valid values (x>0 and x<valdepth) the bottom quality
  logvec_valdepth = and(data_add{end}>=0,data_add{end}<=valdepth);
  qualvec(logvec_valdepth) = layerdata_full.bottom_qual(logvec_valdepth);
  %Make all negative values 3
  qualvec(or(data_add{end}<0,data_add{end}>valdepth)) = 3;
  data_add{end+1} = qualvec;

%Create variables in netcdf file
for did = 1:length(data_add)
  varname = name_add{did};
  tmpvec = data_add{did};
  nccreate(outfn,varname,'Dimensions',{'r',size(tmpvec,1),'c',size(tmpvec,2)});
  ncwrite(outfn,varname,tmpvec);
end