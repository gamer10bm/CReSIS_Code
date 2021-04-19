% script run_preprocess_BAS.m
%
% Support script for run_preprocess.m

param.config.default = [];
singleday = false;
%% ACCUM3 TObas
if singleday
    date_strs = {'20190129','20190131'};
else
    %Make the base directory string
    base_dir = fullfile(gRadar.data_path,'Accum_Data','2018_Antarctica_TObas');
    %Search the base directory for subfolders
    base_files = dir(base_dir);
    base_subdirs = base_files([base_files.isdir]);
    %Grab the subfolder names, truncate to YYYYMMDD, and find unique days
    iddates = 1;
    for idsub = 1:length(base_subdirs)
        if ~isnan(str2double(base_subdirs(idsub).name(1)))
            base_dates{iddates} = base_subdirs(idsub).name;
            iddates = iddates +1;
        end
    end
    date_strs = base_dates;
end

for idx = 1:length(date_strs)
    %Get current date string
    date_str = date_strs{idx};
    %Load the param structure
    cur_idx = length(param.config.default)+1;
    param.config.default{cur_idx} = default_radar_params_2018_Antarctica_TObas_accum();
    param.config.base_dir{cur_idx} = base_dir;
    param.config.config_folder_names{cur_idx} = date_str;
    param.config.board_folder_names{cur_idx} = sprintf('%s%s',date_str,'/%b');
    param.config.date_strs{cur_idx} = date_str;
end
