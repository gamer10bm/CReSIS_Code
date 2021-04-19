%% Make images for March Update report
close all
clearvars -except gRadar
savedir = '/cresis/snfs1/scratch/bmiller/N1KU_Processing/Deployment2/MarchUpdate_Images';
if ~exist(savedir,'dir')
  fprintf('\nMaking directory: %s\n',savedir)
  mkdir(savedir)
end
%Wyoming echogram
run_load_data_by_gps_time_N1KU_Wyoming
saveas(gcf,fullfile(savedir,'Echo_Wyoming_0219.png'))
saveas(gcf,fullfile(savedir,'Echo_Wyoming_0219.fig'))
%Lake echogram
run_load_data_by_gps_time_N1KU_Lake
saveas(gcf,fullfile(savedir,'Echo_Lake_0219.png'))
saveas(gcf,fullfile(savedir,'Echo_Lake_0219.fig'))
%Snotel echogram
run_load_data_by_gps_time_N1KU_Snotel
saveas(gcf,fullfile(savedir,'Echo_Snotel_0225.png'))
saveas(gcf,fullfile(savedir,'Echo_Snotel_0225.fig'))
run_load_data_by_gps_time_N1KU_SnotelZoom
saveas(gcf,fullfile(savedir,'Echo_SnotelZoom_0225.png'))
saveas(gcf,fullfile(savedir,'Echo_SnotelZoom_0225.fig'))

%% Make coverage map
coverage_maps_from_gps_or_records_N1KU
saveas(gcf,fullfile(savedir,'CoverageMap.png'))
saveas(gcf,fullfile(savedir,'CoverageMap.fig'))

%% Zip .fig
zip(fullfile(savedir,'Fig_Files'),'*.fig',savedir)
