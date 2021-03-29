clc
clearvars -except gRadar
close all

%% Check N1KU GPS
gps_dir = fullfile(gRadar.support_path,'gps','2019_SouthDakota_N1KU');
gps_fns = dir(gps_dir);

for id_gps = 1:length(gps_fns)
  gps_path = fullfile(gps_dir,gps_fns(id_gps).name);
  if ~any(strcmp(gps_fns(id_gps).name,{'.','..'}))
    %Load it
%     gps_now = load();
    %Plot it 
    gps = gps_plot(gps_path);
    keyboard
  end
end
