function [gopro_pics] = gopro_preprocess(param)
% function [gopro_pics] = gopro_preprocess(param, param_override)

% The purpose of this function is to grab any gopro pictures and save 
% useful information by each day and segment turned on.

%% General Setup
% =====================================================================

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.season_name, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

if ~isfield(param.records.file,'base_dir') || isempty(param.records.file.base_dir)
  error('Please specify the field param_override.records.file.base_dir')
end

%Turn off tiff warning
tiffwarnid = 'MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero';
warning('off',tiffwarnid)
%% Preprocess all gopro pics
goprodir= fullfile(param.records.file.base_dir,'gopro');

%Get all gopro pictures
filefmt = '.JPG';
cmd_call = sprintf('find %s -name "*%s"',goprodir,filefmt);
[~, sysoutpics] = system(cmd_call);
gopropics = regexp(sysoutpics,'\n','split');

%Get useful information from pics
path={}; date={}; gpsinfo={};
for id_p = 1:length(gopropics)
  if ~isempty(gopropics{id_p}) 
    path{id_p} = gopropics{id_p};
    infonow = imfinfo(path{id_p});
    if ~isfield(infonow,'GPSInfo');
      gpsinfo{id_p} = [];
    else
      gpsinfo{id_p} = infonow.GPSInfo;
    end
    date{id_p} = datetime(infonow.DateTime, 'InputFormat', 'yyyy:MM:dd HH:mm:ss');
    clearvars infonow
  end
end

%% Save data
fprintf('=====================================================================\n');
fprintf('Saving %s GoPro Preprocess: (%s)\n',param.season_name, datestr(now));
fprintf('=====================================================================\n');
gopro_tmp_path = ct_filename_ct_tmp(param,fullfile('gopro',param.season_name,'data.mat'));
gopro_tmp_dir = fileparts(gopro_tmp_path);
if ~exist(gopro_tmp_dir,'dir')
  mkdir(gopro_tmp_dir)
end
save(gopro_tmp_path, 'path', 'date', 'gpsinfo')

%Format output
gopro_pics = struct('path',path,'date',date,'gpsinfo',gpsinfo);
end