% script run_view_frames.m

% The purpose of this script is to run gopro plot functions
clc
close all
clearvars -except gRadar params
%% User Settings
param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));
if ~exist('params','var')
  params = read_param_xls(ct_filename_param('snow_param_2019_SouthDakota_N1KU.xls'));
end

% Example to run a specific segment and frame by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20200131_01');


params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'day_seg','20191211');
params = ct_set_params(params,'cmd.generic',0,'day_seg','20200116');

params = ct_set_params(params,'gopro.qlooks',{'qlook','qlook_noise','qlook_noise_threshold'});

%Define window locations
% start.lat = 43.93; start.lon = -103.98;
% stop.lat = 44.21; stop.lon =  -103.57;

%Load reference locations
refs = struct('lon',[],'lat',[],'name',[],'mark',[]);
snotelmark = struct('type','s','color','g','size',10,'tag','marker');
snocoursemark = struct('type','^','color','m','size',10,'tag','marker');
if 1
%     2 Snotel (North Rapid Creek and Blind Park) and 4 Snow Course (Mallo, Upper Spearfish, Ditch Creek, Little Bear Run). 
    refs(1) = struct('lon',-68.703056,'lat',76.531111,'name','Thule','mark',[]);
    
    refs(end+1) = struct('lat', [44.206188, 44.1077501], 'lon', [-103.787557, -103.9768665], 'name', 'Snotel','mark',snotelmark);
    refs(end+1) = struct('lat', [44.1088951 44.2000038 43.8666741 43.9744436], ...
      'lon', [-104.066101 -103.9999853 -103.783331 -104.055273], 'name', 'Snow Course','mark',snocoursemark);
end

%GIS information
gisfn= fullfile(gRadar.gis_path,'SouthDakota','Model_Landsat.tif');
%Qlook settings
params = ct_set_params(params,'gopro.qlooks',{'qlook','qlook_noise','qlook_noise_threshold'});


%% Automated Section
%Load the records data for the selection window
recdata = init_view_frames(params,gisfn);
%Open the settings window to select the different frames
view_setup(recdata,[],[],refs);