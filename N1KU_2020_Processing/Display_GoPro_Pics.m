clc
clearvars -except gRadar pics
close all

%% GoPro Pics
goprodir= fullfile(gRadar.data_path,'SnowRadar','2020_SouthDakota_CESSNA','gopro');
picdirs = dir(goprodir);

%Get all gopro pictures
if ~exist('pics','var') || 1
  pics = struct('path',[],'info',[]); id_pic = 1;
  for pd = 1:3%length(picdirs)
    dirname = picdirs(pd).name;
    if ~any(strcmp(dirname,{'.','..'}))
      picnow=dir(fullfile(goprodir,dirname));

      for pic_id = 1:length(picnow)
        picname = picnow(pic_id).name;
        if ~any(strcmp(picname,{'.','..'}))
          if any(strcmp(picname(end-3:end),{'.JPG','JPEG'}))
            pics(id_pic).path=fullfile(goprodir,dirname,picname);
            pics(id_pic).info=imfinfo(pics(id_pic).path);
            id_pic=id_pic+1;
          end
        end
      end
    end
  end
end

%Get all of the datetimes for pics
datetimevec = [];
for id_pic = 1:length(pics)
  datetimevec(id_pic) = datetime(pics(id_pic).info.DateTime, ...
    'InputFormat','yyyy:MM:dd HH:mm:ss');
end
%debug
if 1
  figure(10)
  plot(datetimevec)
end