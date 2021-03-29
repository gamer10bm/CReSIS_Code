function snow_thickness_estimate(seg_ovrrde)
global gRadar
% snotelthick = 30; %in
% snotelthick = convlength(snotelthick,'in','m');
physical_constants

% twttlayer = 3.58e-6;
% twttestdiff_air =(2*snotelthick)/c;
% twttestdiff_ice =(2*snotelthick*sqrt(er_ice))/c;
%% Settings 
season_name = '2019_SouthDakota_N1KU'; radar_name = 'snow'; out_pdt = 'qlook_noise'; 
layerData = 'layer'; post_dir = [];

lay_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2019_SouthDakota_N1KU/CSARP_layer';
ql_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2019_SouthDakota_N1KU/CSARP_qlook_noise';
if 0
  area_str = 'Blind Park'; seg_ids = {'20200201_01', '20200128_03', '20200131_02', '20200128_06','20200208_03','20200209_01'};frm_ids = {[34 35], [27 28], [16 17], [45 46],[6 7],[50 51]};
elseif 0
  area_str = 'Duck Lake SW to NE'; seg_ids = {'20200209_02','20200202_06','20200128_05','20200128_02','20200208_02'}; frm_ids = {[43:50] [77:81] [39:43] [43:47] [21:25]};
  seg_dirs = [1 1 -1 1 1];
elseif 0
  area_str = 'Duck Lake SE to NW'; seg_ids = {'20200202_06','20200202_05','20200202_05','20200202_05'}; frm_ids = {[78:79] [62:66] [8:12] [44:48]}; seg_dirs = [1 1 -1 -1];
elseif 0
  area_str = 'Duck Lake Crossovers'; seg_ids = {'20200202_05','20200202_05','20200202_05','20200202_06','20200202_06','20200202_06'}; 
  frm_ids = {[44] [10] [63] [78] [79] [79]}; seg_dirs = [1 1 1 1 1 1]; seg_rlines = {[200:225] [115:140] [35:60] [250:270] [105:130] [217:242]};
elseif 0
  area_str = 'Duck Lake Crossovers'; seg_ids = {'20200202_05','20200202_05','20200202_05','20200202_06'}; 
  frm_ids = {[44] [10] [63] [78:79]}; seg_dirs = [1 1 1 1]; seg_rlines = {[200:225] [115:140] [35:60] []};
elseif 1
  area_str = 'Duck Lake Grid'; seg_ids = {'20200202_06','20200202_05','20200202_05','20200202_05','20200202_05','20200202_07'}; frm_ids = { [75:81] [62:66] [8:12] [44:48],[26:30] [6:12]};
  
end
%% Input Check and Sorting
if ~exist('seg_rlines','var')
  seg_rlines = cell(1,length(seg_ids));
end
if ~exist('frm_ids','var')
  frm_ids = cell(1,length(seg_ids));
end
if ~exist('seg_dirs','var')
  seg_dirs = ones(1,length(seg_ids));
end
%Sort the seg_ids
segs = cellfun(@(c)str2num(c(1:8)),seg_ids);
[~,sortind] = sort(segs,'ascend');
seg_ids = seg_ids(sortind); frm_ids = frm_ids(sortind); seg_dirs = seg_dirs(sortind);

layer_dats = []; lyd_id = 0;
%% Load the layer data and qlook data
if ~exist('seg_ovrrde','var')
  seg_vec = 1:length(seg_ids);
else
  seg_vec = seg_ovrrde;
end
for sid = seg_vec
  lyd_id = lyd_id +1;
  seg_id = seg_ids{sid};
  seg_laydir = fullfile(lay_dir,seg_id);
  seg_qldir = fullfile(ql_dir,seg_id);
  %Load the layer format  file
  layfrmt = load(fullfile(seg_laydir,sprintf('layer_%s.mat',seg_id)));
  frms = frm_ids{sid};
  for frmnow = frms
    %Make fns for this frame
    lay_fn = fullfile(seg_laydir,sprintf('Data_%s_%03.0f.mat',seg_id,frmnow));
    ql_fn = fullfile(seg_qldir,sprintf('Data_%s_%03.0f.mat',seg_id,frmnow));
    %% Layer load
    laydat = load(lay_fn);
    if frmnow == frms(1) 
      if lyd_id == 1
        layer_dats = struct();
      end
    end
    %Load layers to layer_dats
    for lay_id = 1:length(layfrmt.lyr_name);
      layname = layfrmt.lyr_name{lay_id};
      if frmnow == frms(1)
        layer_dats(lyd_id).(layname)=[];
      end
      layer_dats(lyd_id).(layname)=[layer_dats(lyd_id).(layname) laydat.twtt(lay_id,:)];
    end
    %Load the gps_time
    gps_tmp = laydat.gps_time; 
    lat_tmp = laydat.lat; lon_tmp = laydat.lon;
    frm_tmp = frmnow*ones(1,length(gps_tmp));
    if frmnow == frms(1)
      layer_dats(lyd_id).gps_time = gps_tmp;
      layer_dats(lyd_id).lat = lat_tmp;
      layer_dats(lyd_id).lon = lon_tmp;
      layer_dats(lyd_id).frms = frm_tmp;
    else
      layer_dats(lyd_id).gps_time = [layer_dats(lyd_id).gps_time gps_tmp];
      layer_dats(lyd_id).lat = [layer_dats(lyd_id).lat lat_tmp];
      layer_dats(lyd_id).lon = [layer_dats(lyd_id).lon lon_tmp];
      layer_dats(lyd_id).frms = [layer_dats(lyd_id).frms  frm_tmp];
    end
    %% qlook load
    qldat = load(ql_fn);
    if frmnow == frms(1)
      if lyd_id == 1
        ql_dats = struct('data',qldat.Data,'gps_time',qldat.GPS_time,'time',qldat.Time,'lat',qldat.Latitude,'lon',qldat.Longitude);
      else
        ql_dats(lyd_id) = struct('data',qldat.Data,'gps_time',qldat.GPS_time,'time',qldat.Time,'lat',qldat.Latitude,'lon',qldat.Longitude);
      end
    else
      if size(ql_dats(lyd_id).data,1) < size(qldat.Data,1)
        dat_tmp = zeros(size(qldat.Data,1),size(ql_dats(lyd_id).data,2)+size(qldat.Data,2));
        dat_tmp(1:size(ql_dats(lyd_id).data,1),1:size(ql_dats(lyd_id).data,2)) = ql_dats(lyd_id).data;
        dat_tmp(:,size(ql_dats(lyd_id).data,2)+1:end)=qldat.Data;
      else
        dat_tmp = [ql_dats(lyd_id).data(:,:) qldat.Data(:,:)];
      end
      ql_dats(lyd_id) = struct('data',dat_tmp,'gps_time',[ql_dats(lyd_id).gps_time qldat.GPS_time],...
        'time',qldat.Time,'lat',[ql_dats(lyd_id).lat(:,:) qldat.Latitude],'lon',[ql_dats(lyd_id).lon(:,:) qldat.Longitude]);
    end
  end
  %Load special stuff
  layer_dats(lyd_id).seg_id = seg_id;
  layer_dats(lyd_id).dir = seg_dirs(sid);
  layer_dats(lyd_id).rlines = seg_rlines{sid};
  layer_dats(lyd_id).frm_ids = frms;
end


%% Take the difference between surface and bottom
botvecs = cell(1,length(layer_dats)); survecs = cell(1,length(layer_dats));
alongvecs = cell(1,length(layer_dats)); depthvecs = cell(1,length(layer_dats));
gpstimevecs = cell(1,length(layer_dats)); diffvecs = cell(1,length(layer_dats));
latvecs = cell(1,length(layer_dats)); lonvecs = cell(1,length(layer_dats));
depthvecs_filt = cell(1,length(layer_dats));
seg_dats = cell(1,length(layer_dats)); 
frm_ids = cell(1,length(layer_dats));
for lid =  1:length(layer_dats)
  nonnan_logvec = [];
  %Get rid of all nan values
  if ~isempty(layer_dats(lid).rlines)
    ind_strt = find(layer_dats(lid).gps_time>=ql_dats(lid).gps_time(layer_dats(lid).rlines(1)),1,'first');
    ind_end = find(layer_dats(lid).gps_time<=ql_dats(lid).gps_time(layer_dats(lid).rlines(end)),1,'last');
  else
    ind_strt = find(~isnan(layer_dats(lid).bottom),1,'first');
    ind_end = find(~isnan(layer_dats(lid).bottom),1,'last');
  end
  nonnan_logvec = ind_strt:ind_end;
  botvecs{lid} = layer_dats(lid).bottom(nonnan_logvec);
  survecs{lid} = layer_dats(lid).surface(nonnan_logvec);
  latvecs{lid} = layer_dats(lid).lat(nonnan_logvec);
  lonvecs{lid} = layer_dats(lid).lon(nonnan_logvec);
  alongvec_tmp = geodetic_to_along_track(latvecs{lid},lonvecs{lid});
  alongvecs{lid} = alongvec_tmp-min(alongvec_tmp);
  gpstimevecs{lid} = layer_dats(lid).gps_time(nonnan_logvec);
  frm_ids{lid} = layer_dats(lid).frms(nonnan_logvec);
  %Determine snow depth
  diffvecs{lid} = abs(survecs{lid}-botvecs{lid});
  %Load seg_dats with qldata
  inds_strt = find(ql_dats(lid).gps_time>=min(gpstimevecs{lid}),1,'first');
  inds_end = find(ql_dats(lid).gps_time<=max(gpstimevecs{lid}),1,'last');
  inds_seg = inds_strt:inds_end;
  %Modulate for start and finish
  %Make lat_lon_str
  llid = 0; lat_lon_strs = cell(1,length(inds_seg));
  for ind_id = inds_seg
    llid = llid +1;
    lat_lon_frmt = '%.3fN%s%.3fE';
    lat_lon_strs{llid} = sprintf(lat_lon_frmt,ql_dats(lid).lat(ind_id),'\n', ql_dats(lid).lon(ind_id));
  end
    
  seg_dats{lid} = struct('data',ql_dats(lid).data(:,inds_seg),...
    'gps_time',ql_dats(lid).gps_time(inds_seg),'time',ql_dats(lid).time,...
    'lat_lon_strs',{lat_lon_strs},'frm_ids',frm_ids{lid},'seg_id',layer_dats(lid).seg_id,...
    'lat',ql_dats(lid).lat(inds_seg),'lon',ql_dats(lid).lat(inds_seg));
%   seg_dats{lid}.lat_lon_strs=lat_lon_strs;
  %Estimate snowdepth for this day
  depthvecs{lid} = convlength(diffvecs{lid}*c/(2*sqrt(er_ice)),'m','in'); %in
  %Do some filtering
  fillength = 40;
  while fillength > length(depthvecs{lid})
    fillength = fillength-4;
  end
  for did = 1:length(depthvecs{lid})
    windstrt = did-floor(fillength/2);
    windend = did+floor(fillength/2);
    if windstrt < 1
      windstrt = 1;
    end
    if windend>length(depthvecs{lid})
      windend = length(depthvecs{lid});
    end
    wind_ind = windstrt:windend;
    %Take average of window and log it
    depthvecs_filt{lid}(did) = mean(depthvecs{lid}(wind_ind));
  end
  %Post result
  %make frm string
  unq_frms = unique(seg_dats{lid}.frm_ids);
  frm_str = '';
  for f_id = 1:length(unq_frms)
    frm_str = sprintf('%s%.0f,',frm_str,unq_frms(f_id));
  end
  if lid == 1
    fprintf('\nSegment\tFrames\tAverage Depth (in)\tAlong Track (km)\t\tArea Description\n')
  end
  fprintf('\n%s\t[%s]\t%.2f\t%.2f\t\t%s',layer_dats(lid).seg_id,frm_str(1:end-1),nanmean(depthvecs{lid}),alongvecs{lid}(end)*1e-3,area_str)
end
fprintf('\n\n')

%% Determine any crossovers
for sid1 = 1:length(seg_dats)
  latref = seg_dats{sid1}.lat;
  lonref = seg_dats{sid1}.lon;
  crossovers = cell(1,length(seg_dats));
  for sid2 = 1:length(seg_dats)
    %Determine any crossovers for this segment
    latchk = seg_dats{sid2}.lat; lonchk = seg_dats{sid2}.lon;
    minnow = .1; minid = 0;
    xfrm = 0; xgps_time = 0; xseg_id =[];
    for llid = 1:length(latref)
      distval = sqrt((latref(llid)-latchk).^2+(lonref(llid)-lonchk).^2);
      if any(distval<.000005) && min(distval)<minnow
        [minnow, minind] = min(distval);
        xfrm = seg_dats{sid2}.frm_ids(minind);
        xgps_time = seg_dats{sid1}.gps_time(llid);
        xseg_id = seg_dats{sid2}.seg_id;
      end
    end
    crossovers{sid2} = struct('seg_id',xseg_id,'gps_time',xgps_time,'frm_id',xfrm);
  end
  seg_dats{sid1}.crossovers = crossovers;
end

%% Make 2D line plot for color of depths
if 1
  figure(1000);clf;
  %plot gis map
  %GIS information
  gisfn= fullfile(gRadar.gis_path,'SouthDakota','Model_Landsat.tif');
  gisdata = geotiff_read(gisfn);
  mapim = imagesc(gisdata.X, gisdata.Y, gisdata.C);
  hold on;
  %Flip y dir
  set(gca,'YDir','normal')
  minlat = nan; minlon = nan; maxlat =nan; maxlon = nan;
  for sid = 1:length(seg_dats)
    lon = lonvecs{sid}' ; %// extract "X" column
    lat = latvecs{sid}' ; %// same for "Y"
    depth = depthvecs_filt{sid}' ; %// extract color index for the custom colormap
    %Extract bounds
    if isnan(minlat) || min(lat)<minlat
      minlat = min(lat);
    end
    if isnan(minlon) || min(lon)<minlon
      minlon = min(lon);
    end
    if isnan(maxlat) || max(lat)>maxlat
      maxlat = max(lat);
    end
    if isnan(maxlon) || max(lon)>maxlon
      maxlon = max(lon);
    end
    %% // Prepare matrix data
    lonmat=[lon lon];           %// create a 2D matrix based on "X" column
    latmat=[lat lat];           %// same for Y
    zz=zeros(size(lonmat)); %// everything in the Z=0 plane
    depthvec =[depth depth] ;         %// matrix for "CData"

    %// draw the surface (actually a line)
    hs=surf(lonmat,latmat,zz,depthvec,'EdgeColor','interp','FaceColor','none') ;

    colormap('spring') ;     %// assign the colormap
    shading flat                    %// so each line segment has a plain color
    view(2) %// view(0,90)          %// set view in X-Y plane
  end
  %format figue
  xlabel('Longitude (^oE)')
  ylabel('Latitude (^oN)')
  title(sprintf('Depth Map Near %s',area_str))
  set(findall(gca,'Type','Surface'),'Linewidth',5)
  
  c = colorbar;
  c.Label.String = 'Snow Depth (in.)';
  c.FontSize = 12;
  caxis([10 60])
  lonoff = .02; latoff = .02;
  ylim([minlat-latoff maxlat+latoff]);
  xlim([minlon-lonoff maxlon+lonoff]);
  set(findall(gcf,'Type','Text'),'Fontsize',12)
end
%% Make layer comparison figure
if 0
  figure(100); clf; leg = {};
  for lid =  1:length(layer_dats)
    %Post to figure
    plot(alongvecs{lid},depthvecs_filt{lid})
    leg{lid} = strrep(layer_dats(lid).seg_id,'_',' ');
    hold on
  end
  hold off
  %Format figure
  title(sprintf('Snow Depths near %s',area_str))
  xlabel('Along Track (m)')
  ylabel('Snow Depth (in)')
  legend(leg)
  set(findall(gcf,'Type','text'),'Fontsize',12)
  grid on
end
%% Make figures for the qlook data and the layers
if 1
  for sid = 1:length(seg_dats)
    %make frm string
    unq_frms = unique(seg_dats{sid}.frm_ids);
    frm_str = '';
    for f_id = 1:length(unq_frms)
      frm_str = sprintf('%s%.0f,',frm_str,unq_frms(f_id));
    end
    %Make new figure
    figure(100+sid); clf; subm = 2; subn = 1;
    yoff = .2e-7;
    %Plot the qlook with layers
    subplot(subm,subn,1)
    imagesc(seg_dats{sid}.gps_time, seg_dats{sid}.time, lp(seg_dats{sid}.data))
    colormap(flipud(gray))
    hold on
    leg = {};
    plot(gpstimevecs{sid},survecs{sid},'--')
    leg{end+1} = 'Surface';
    plot(gpstimevecs{sid},botvecs{sid},'--')
    leg{end+1} = 'Bottom';
    %Plot crossovers
    for sid2 = 1:length(seg_dats{sid}.crossovers)
      cross_now = seg_dats{sid}.crossovers{sid2};
      if ~isempty(cross_now.seg_id) && ~strcmp(cross_now.seg_id,seg_dats{sid}.seg_id)
        gpsvec = cross_now.gps_time;
        xvec = []; yvec = [];
        for xid = 1:length(gpsvec)
          xvec = [xvec gpsvec(xid) gpsvec(xid)];
          if rem(xid,2) == 1
            yvec = [yvec 0 1];
          else
            yvec = [yvec 1 0];
          end
        end
        %Plot the crossover
        plot(xvec,yvec)
        leg{end+1} = sprintf('Cross Seg: %s Frm: %.0f',strrep(cross_now.seg_id(5:end),'_',' '),cross_now.frm_id);
      end
    end
    hold off
    ylim([min(survecs{sid})-yoff max(botvecs{sid})+yoff])
    legend(leg,'Location','Best')
    title(sprintf('Segment: %s\nFrame: %s',strrep(seg_dats{sid}.seg_id,'_',' '),frm_str(1:end-1)))
    xlabel('Lat (^oN) Lon (^oE)')
    ylabel('Fast Time (sec)')

    %Plot the qlook with layers
    subplot(subm,subn,2)
    imagesc(seg_dats{sid}.gps_time, seg_dats{sid}.time, lp(seg_dats{sid}.data))
    colormap(flipud(gray))
    ylim([min(survecs{sid})-yoff max(botvecs{sid})+yoff])
    xlabel('Lat (^oN) Lon (^oE)')
    ylabel('Fast Time (sec)')

    set(findall(gcf,'Type','Line'),'Linewidth',2)
    %Change the tick labels
    axs = findall(gcf,'type','axes');
    for a_id = 1:length(axs)
      xtickval = axs(a_id).XTick;
      xtickstr = {}; xtickvalnew = []; nid = 0;
      for xt_id = 1:length(xtickval)
        nid = nid +1;
        if any(xt_id==1:2:length(xtickval))
          str_ind = find(seg_dats{sid}.gps_time>=xtickval(xt_id),1,'first');
          if isempty(str_ind)
            str_ind = length(seg_dats{sid}.gps_time);
          end
          str_now = seg_dats{sid}.lat_lon_strs{str_ind};
          brk_ind = strfind(str_now,'\n');
          if ~isempty(brk_ind)
            xtickstr{nid} = sprintf(str_now(1:brk_ind-1));
            xtickvalnew(nid) = xtickval(xt_id);
            nid = nid +1;
            xtickstr{nid} = sprintf(str_now(brk_ind+2:end));
            xtickvalnew(nid) = xtickvalnew(nid-1)+.15*median(diff(xtickval));
          else
            xtickstr{nid} = sprintf(str_now);xtickvalnew(nid) = xtickval(xt_id);
          end
        else
          xtickstr{nid} = ' ' ;
          xtickvalnew(nid) = xtickval(xt_id);
        end
      end
      %Make pretty label
      set(axs(a_id),'XTickLabelMode','manual')
      set(axs(a_id),'XTick',xtickvalnew)
      set(axs(a_id),'XTickLabel',xtickstr);
      set(axs(a_id),'XTickLabelRotation',45);
      if layer_dats(sid).dir == -1
        set(axs(a_id),'XDir','reverse')
      end
      %Make expressive tick labels
    end
    %Make figure fullscreen
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %Save figure as .png
    basedir = '~/scripts/matlab/N1KU_Processing';
    areadir = fullfile(basedir,strrep(area_str,' ','_'));
    if ~exist(areadir,'dir')
      mkdir(areadir)
    end
    seg_fn = fullfile(areadir,sprintf('Echo_%s__%s',seg_dats{sid}.seg_id,strrep(frm_str(1:end-1),',','_')));
    saveas(gcf,[seg_fn '.png'])
    saveas(gcf,[seg_fn '.fig'])
  end
end