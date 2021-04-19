function recdata = view_frames(recdata,start,stop,refs)
set_fig = findall(0,'type','figure','tag','Setup');
%Change the pointer
oldpoint = waitpoint(set_fig);
%Get view_frames figure
vwfrm_tag = 'View Frames';
vwfrm_fig = findall(0,'type','figure','tag',vwfrm_tag);
if isempty(vwfrm_fig) || 1
  %Figure is not open
  vwfrm_fig = init_view(recdata,start,stop,refs,vwfrm_tag);
end

% Check initial settings figure
chksettings(true)

waitpoint(set_fig,oldpoint);
end

%% Data access functions
function [recdata, start, stop, refs, selections] = savefigdata(recdata,start,stop,refs,selections)
%Send recdata to userdata
set(gcf,'UserData',recdata);
%Check for selections
if ~exist('selections','var') || isempty(selections) || ~isfield(selections,'seg_id')
  %Initialize selection
  selections.seg_id = {recdata(1).id};
  selections.param = {recdata(1).param};
  selections.frames = {floor(recdata(1).Frame_IDs(1))};
  selections.lat = recdata(1).Latitude(1);
  selections.lon = recdata(1).Longitude(1);
  selections.gpstime = recdata(1).GPS_time(1);
end
gdat = struct('start',start,'stop',stop,'refs',refs,'selections',selections);
guidata(gcf,gdat);
end

function [recdata, start, stop, refs, selections] = getfigdata()
%Get recdata from userdata
recdata = get(gcf,'UserData');
%Get the rest from guidata
gdat = guidata(gcf);
start = gdat.start; stop = gdat.stop; refs = gdat.refs; selections = gdat.selections;
end

function [botlayer, surlayer] = loadlayers(selections)
if (~exist('selections','var') || isempty(selections))
  [~, ~, ~, ~, selections] = getfigdata;
end
%Get layer directory
laydir = ct_filename_out(selections.param{1},'','layer'); layfnfmt = 'Data_%s_%03.0f.mat';
%Grab layer data and determine height AGL
layfn = fullfile(laydir,sprintf(layfnfmt,selections.seg_id{1},selections.frames{1}));
if exist(layfn,'file')
  layer = load(layfn);
else
  layer = [];
end
layalong = geodetic_to_along_track(layer.lat,layer.lon);
%Plot the surface
botid = 1; surid = 2;
botlayer = struct('alongtrack',layalong,'lat',layer.lat,'lon',layer.lon,'twtt',layer.twtt(botid,:));
surlayer = struct('alongtrack',layalong,'lat',layer.lat,'lon',layer.lon,'twtt',layer.twtt(surid,:));
end

%% Initialization Functions

function vwfrm_fig = init_view(recdata,start,stop,refs,vwfrm_tag)
% Make and format the figure
vwfrm_fig = figure('numbertitle','off','name',vwfrm_tag,...
    'Units','normalized','outerposition',[0 0 1 1],'tag',vwfrm_tag);
set(vwfrm_fig,'CloseRequestFcn',@(s,e)close_view(s,e))
%Change fig to watch pointer
oldpoint = waitpoint(vwfrm_fig);
% Save the initial guidata
[recdata, start, stop, refs, selections] = savefigdata(recdata,start,stop,refs);

% Make gopro uipanel
gopan_pos = [0 0 .5 .5];
gopro_pan = uipanel(vwfrm_fig,'Units','normalized','Position',gopan_pos,'Tag','gopro','Title','GoPro');
gopro_ref = init_gopro(gopro_pan,recdata,selections);
%Add gopro to refs
refs(end+1) = gopro_ref;

% Make qlook uipanel
qlpan_pos = [.5 0 .5 .5];
ql_pan = uipanel(vwfrm_fig,'Units','normalized','Position',qlpan_pos,'Tag','qlook','Title','Qlook');
init_qlook(ql_pan,recdata,refs,selections);

% Make map and settings panel (must be last)
setpan_pos = [0 .5 .25 .5];
% set_pan = uipanel(vwfrm_fig,'Units','normalized','Position',setpan_pos,'Tag','settings','Title','Settings');
% init_set(set_pan)

mappan_pos = [0 .5 1 .5];
if exist('setpan_pos','var') && exist('set_pan','var') && ~isempty(setpan_pos) && isnum(setpan_pos(1))
  mappan_pos = [setpan_pos(1)+setpan_pos(3) mappan_pos(2) mappan_pos(3)-setpan_pos(3) mappan_pos(4)];
end
map_pan = uipanel(vwfrm_fig,'Units','normalized','Position',mappan_pos,'Tag','map','Title','Map');
init_map(map_pan,recdata,start,stop,refs,selections);


%Change point back
waitpoint(vwfrm_fig,oldpoint);
end

function gopro_ref = init_gopro(gopro_pan,recdata,selections)
%Make gopro axes
ax_pos = [0 0 1 1];
go_tag = 'gopro';
go_ax = axes('Parent',gopro_pan,'Tag',go_tag);
%Plot the selections
gopro_ref = plotgopro(recdata,selections);
%Format the gopro axis
set(go_ax,'Tag',go_tag,'outerposition',ax_pos);
end

function init_qlook(ql_pan,recdata,refs,selections)
%Make qlook axes
ax_pos = [0 0 1 1];
ql_tag = 'qlook';
ql_ax = axes('Parent',ql_pan,'Tag',ql_tag);
%Plot the selections
plotqlook(refs,selections)
%Format the qlook axis
set(ql_ax,'Tag',ql_tag,'outerposition',ax_pos)
end


function init_map(map_pan,recdata,start,stop,refs,selections)
% Map the geotiff (if available)
ax_pos = [0 0 .75 1];
if isfield(recdata(1).param,'gisfn')
  gisdata = geotiff_read(recdata(1).param.gisfn);
  mapax = axes('Parent',map_pan,'Tag','map');
  mapim = imagesc(gisdata.X, gisdata.Y, gisdata.C);
  set(mapim,'PickableParts','none','Tag','map')
  hold('on'); grid('on');
  %Adjust axes settings
  set(mapax,'YDir','normal','Units','Normalized','outerposition',ax_pos)
end
% Make selector panel
sel_tag = 'map_selector';
sel_pos = [ax_pos(1)+ax_pos(3) ax_pos(2) 1-ax_pos(3) ax_pos(4)];
make_selector(map_pan,recdata,sel_pos,sel_tag);
%Plot the frame
[axlim] = plotframes(recdata,start,stop,refs,selections);
%Set mapax tag for full effect
set(mapax,'Tag','map')
end

function init_set(set_pan)
qlook_pass = {'but'};
button_qlook(set_pan,qlook_pass);
button_gopro(set_pan);
end
%% Gopro button group function
function go_out = button_gopro(gopro_pan)
%Make a gopro selection panel
but_pos = [0 .65 1 .2];

go_out.bg_gopro = uipanel('Parent',gopro_pan, 'Position', but_pos,...
  'Title','GoPro Options','FontSize', 12);

pos_pos = [0.05 .6 .75 .3];
go_out.pos_but = uicontrol(go_out.bg_gopro, 'Style', 'checkbox', ...
  'String', 'Position Display','Units','normalized','Position',pos_pos);

flt_pos = [0.05 .2 .75 .3];
go_out.flt_but = uicontrol(go_out.bg_gopro, 'Style', 'checkbox', ...
  'String', 'Flightline Display','Units','normalized','Position',flt_pos);

go_butlist = {go_out.pos_but, go_out.flt_but};
tog_pos = [.8 .4 .2 .5];
go_out.tog_gopro = uicontrol(go_out.bg_gopro,'Style','togglebutton',...
  'FontSize', 8,'String','On','Value',1,'Callback',@(h,e)toggle_group(h,e,go_butlist),...
  'Units','normalized','Position', tog_pos);
end

%% Qlook button group function
function qlook_out = button_qlook(set_pan, qlook_pass)
%Make a qlook selection panel
but_pos = [0 .35 1 .3];

qlook_out.bg_qlook = uipanel('Parent',set_pan, 'Position', but_pos,...
  'Title','qLook Options','FontSize', 12);

drop_pos = [0.05 .15 .8 .2];
stringer ={'Select...', qlook_pass{:}};
qlook_out.qlook_drop = uicontrol(qlook_out.bg_qlook, 'Style', 'popupmenu',...
  'String', stringer,'Units','normalized','Position',drop_pos);

surf_pos = [0.05 .7 .9 .2];
qlook_out.surface_but = uicontrol(qlook_out.bg_qlook, 'Style', 'checkbox',...
  'String', 'Surface Display','Units','normalized','Position',surf_pos);

camv_pos = [0.05 .5 .9 .2];
qlook_out.camview_but = uicontrol(qlook_out.bg_qlook, 'Style', 'checkbox',...
  'String', 'CamView Display','Units','normalized','Position',camv_pos);

tog_pos = [.8 .55 .2 .35];
ql_butlist = {qlook_out.camview_but, qlook_out.surface_but, qlook_out.qlook_drop};
if isempty(qlook_pass)
  tog_fcn = @(h,e)warning('No qlooks available');
else
  tog_fcn = @(h,e)toggle_group(h,e,ql_butlist);
end
qlook_out.tog_qlook = uicontrol(qlook_out.bg_qlook,'Style','togglebutton',...
  'FontSize', 8,'String','On','Value',1,'Callback',tog_fcn,...
  'Units','normalized','Position', tog_pos);

end

%% Plotting and Quit Button Group Function
function start_out = button_start(set)

end

%% Button group toggle
function toggle_group(h,e,but_list)
%Determine if valued changed
if strcmp(h.String,'On')
  h.String = 'Off';
else
  h.String = 'On';
end
%Turn all buttons to handle value
for bid = 1:length(but_list)
  but_list{bid}.Visible = lower(h.String);
end
end

%% Map plotting functions
function [axlim, lonvec, latvec] = plotframes(recdata,start,stop,refs,selections)
%Input check
if ~exist('recdata','var') || ~exist('start','var') || ~exist('stop','var') || ~exist('selections','var')
  [recdata, start, stop, refs, selections] = getfigdata();
end
%Get the map axes
vwfrm_fig = findall(0,'Type','figure','Tag','View Frames');
map_ax = findall(vwfrm_fig,'type','axes','tag','map');
%Set to current axes
set(vwfrm_fig,'CurrentAxes',map_ax);
%Modulate start and stop
if start.lat >stop.lat
  maxlat = start.lat; minlat = stop.lat;
else
  minlat = start.lat; maxlat = stop.lat;
end
if start.lon >stop.lon
  maxlon = start.lon; minlon = stop.lon;
else
  minlon = start.lon; maxlon = stop.lon;
end
%Defined axis offsets
offperc = .5;
latoff = (maxlat-minlat)*offperc;
lonoff = (maxlon-minlon)*offperc;
%Make the window
latvec = [minlat, maxlat, maxlat, minlat, minlat];
lonvec = [minlon, minlon, maxlon, maxlon, minlon];
%Generate the axis limits
axlim = [minlon-lonoff, maxlon+lonoff, minlat-latoff, maxlat+latoff];
%Set hold to on
hold('on')
%Plot the segments
pltinit = struct('Longitude',[],'Latitude',[],'Elevation',[],'GPS_time',[],...
  'Frame_IDs',[],'param',[],'filename',[],'id',[]);
param_season = cell(length(recdata),1); frm_leg = {};
colors2choose = {[0 0 255],  [0 255 0], [255 255 0], [7 252 255]}; chid =0;
for recid = 1:length(recdata)
  %Get season
  param_season{recid} = recdata(recid).param.season_name;
  %Reinitialize
  segtmp = []; plttmp = pltinit;
  %Load the data
  segtmp = recdata(recid);
  %Get logical vector of valid points
  window_log_vec = (segtmp.Latitude>minlat & segtmp.Latitude<maxlat & segtmp.Longitude>minlon & segtmp.Longitude<maxlon);
  %Add the full frames that may extend past the geowin
  valid_frms = unique(floor(segtmp.Frame_IDs(window_log_vec)));
  %Iterate through frames and make plttmp variable for segment
  frm_log_vec = zeros(1,length(segtmp.Frame_IDs));
  for frm_idx = valid_frms
    frm_log_tmp = (segtmp.Frame_IDs>frm_idx & segtmp.Frame_IDs<frm_idx+1);
    frm_log_vec = frm_log_vec|frm_log_tmp;
  end
  %Modulate plotting based on frames in window
  seg_chkbx = findall(vwfrm_fig,'Style','checkbox','Tag',recdata(recid).id);
  if ~isempty(seg_chkbx) && any(frm_log_vec) && ~seg_chkbx.Value
    %Reset chkbx text
    set(seg_chkbx,'String',recdata(recid).id)
    %Modulate segtmp based on logical vector
    plttmp = struct('Longitude',segtmp.Longitude(frm_log_vec),'Latitude',segtmp.Latitude(frm_log_vec),'Elevation',segtmp.Elevation(frm_log_vec),...
      'GPS_time',segtmp.GPS_time(frm_log_vec),'Frame_IDs',floor(segtmp.Frame_IDs(frm_log_vec)),'param',segtmp.param,'filename',segtmp.filename,'id',segtmp.id);

    %Plot each segment
    gn = plot(gca,plttmp.Longitude,plttmp.Latitude,'.');
    set(gn,'ButtonDownFcn',@(s,e)checkclick(s,e),'UserData',plttmp,...
      'Tag',plttmp.id,'LineWidth',2,'HandleVisibility','off')
    hold on
    
    %Check for adding legend entry and season color
    if isempty(frm_leg) || ~any(cellfun(@(c)strcmp(param_season{recid},c),frm_leg))
      set(gn,'HandleVisibility','on');
      frm_leg{end+1} = param_season{recid};
      chid = chid + 1;
    end
    %Set color
    set(gn,'Color',colors2choose{chid}./255);
        
    %Highlight selected frames and make unpickable
    chk_select = cellfun(@(c)strcmp(c,plttmp.id),selections.seg_id);
    if any(chk_select)
      %Plot each selected frame
      selec_frms = selections.frames{chk_select};
      for frm_id = selec_frms
        sel_log_vec = (plttmp.Frame_IDs == frm_id);
        %Find frame within plot tmp
        plot(gca,plttmp.Longitude(sel_log_vec),plttmp.Latitude(sel_log_vec),...
          'LineWidth',1.5,'Pickableparts','none','Color','r','Tag','Select','HandleVisibility','off')
      end
    end
  else
    %Change box to hidden
    hide_str = '(hidden)';
    if ~strcmp(seg_chkbx.String(end-length(hide_str)+1:end),hide_str)
      set(seg_chkbx,'String',sprintf('%s %s',seg_chkbx.String,hide_str))
    end
    set(seg_chkbx,'Value',1)
  end
end
%Turn Hold off
hold('off')
%Make the legend entries
for lid = 1:length(frm_leg)
  frm_leg{lid} = strrep(frm_leg{lid},'_',' ');
end
if isempty(frm_leg)
  wid = 'MATLAB:legend:PlotEmpty';
  warning('off',wid)
end
%Plot Reference location
hold on
ref_leg = plotrefs(refs,axlim);
%Format the plot
make_legend(ref_leg,frm_leg);
xlabel('Longitude (^o E or W)')
ylabel('Latitude (^o N or S)')
hold off
%Move map down
mapim = findall(gcf,'Type','image','Tag','map');
if isempty(mapim)
  uistack(mapim,'bottom')
end
drawnow
end

function [ref_leg] = plotrefs(refs,axlim)
%Delete any current refernce markers
delete(findall(gca,'Type','Line','Tag','Ref'));
%Plot refs based on update
if exist('refs','var') && ~isempty(refs(1).lat)
  ref_leg = {};
  for id_r = 1:length(refs)
    ref = refs(id_r);
    %Check if within viewing window
    if any(axlim(1)<ref.lon & ref.lon<axlim(2) & axlim(3)<ref.lat & ref.lat<axlim(4))
      %Plot based on settings
      if ~isempty(ref.mark)
        rp = plot(ref.lon,ref.lat,ref.mark.type,'Color',ref.mark.color);
        switch ref.mark.tag
          case 'marker'
            set(rp,'MarkerSize',ref.mark.size,'LineWidth', 2.5,...
              'Color','k','MarkerFaceColor',ref.mark.color);
          case 'line'
            set(rp,'LineWidth',ref.mark.size);
        end
      else
        rp = plot(ref.lon,ref.lat,'*');
      end
      %Add ref tag
      set(rp,'Tag','ref');
      ref_leg{end+1} = ref.name;
    end
  end
end
end

function make_legend(ref_leg, frm_leg)
%Make cell array of legend strings
legnow_str = frm_leg; rid = 0;
while rid < length(ref_leg)
  rid = rid +1;
  legnow_str{end+1} = ref_leg{rid};
end
%Check for old legend
legtag = 'map';
legs = findall(gcf,'Type','legend','Tag',legtag);
if ~isempty(legs)
  leg_str = legs.String;
  make_bool = all(cellfun(@(c)any(strcmp(c,leg_str)),legnow_str));
else
  make_bool = true;
end

if make_bool
  [legnow,legobjs] = legend(legnow_str,'Fontsize',12,'Location','Best','Visible','on','Tag',legtag);

  %Get frm_leg entries
  ch = [];
  for l_id = 1:length(legnow.String(1:end-length(ref_leg)))
    chnow = findobj(legobjs,'Type','line','Tag',legnow.String{l_id});
    ch = [ch chnow];
  end
  %Adjust frm_leg entry marker size
  set(ch,'Markersize',35);
end
end

function [] = make_selector(map_pan,recdata,sel_pos,sel_tag)
%Make selection panel for displaying lines
main_pan = uipanel(map_pan,'Units','Normalized','Position',sel_pos,...
  'Title','Select frames to hide...','Tag',sel_tag);
%Make tabgroup
tabgrp = uitabgroup(main_pan);
%Make a mini panel for each unique season
param_season = cell(length(recdata),1);
for r_id = 1:length(recdata)
  param_season{r_id} = recdata(r_id).param.season_name;
end
unq_seasons = unique(param_season);
unq_pans = []; season_strs = cell(1,length(unq_seasons));
for uid = 1:length(unq_seasons)
  u = uitab(tabgrp,'Units','Normalized','Title',strrep(unq_seasons{uid},'_',' '));
  unq_pans = [unq_pans u];
  %Find segments for each season
  season_strs{uid} = {};
  for r_id = 1:length(recdata)
    
    if strcmp(recdata(r_id).param.season_name,unq_seasons{uid})
      season_strs{uid}{r_id} = recdata(r_id).id;
    end
  end
end

%Make checkboxes for each season
cbfunc = @(s,e)updateframes(s,e);
for uid = 1:length(season_strs)
  numcols = 2;
  [seg_strs, seg_logvec] = sort(season_strs{uid});
  if numcols > length(seg_strs)
    numcols = length(seg_strs);
  end
  for colid = 1:numcols
    seg_str_now = seg_strs(colid:numcols:end);
    for cid = 1:length(seg_str_now)
      %Make unique position
      posnow = [(.025+1/numcols*(colid-1)) (.95-1/length(seg_str_now)*(cid-1))  (1/numcols) .05];
      %Make checkbox
      cbnow = uicontrol(unq_pans(uid),'Style','checkbox','String',seg_str_now{cid},'Units','Normalized',...
        'Position',posnow,'Value',0,'Callback',cbfunc,'Tag',seg_str_now{cid});
    end
  end
end
end

%% Qlook plotting functions
function plotqlook(refs,selections)
if (~exist('selections','var') || isempty(selections)) || (~exist('refs','var') || isempty(refs))
  [~, ~, ~, refs, selections] = getfigdata;
end
%Get the qlook axes
ql_tag = 'qlook';
ql_ax = findall(0,'Type','axes','Tag',ql_tag);
%Set to current axes
set(findall(0,'Type','figure','tag','View Frames'),'CurrentAxes',ql_ax);
%Get the qlook_type from the dropdown (Make soon)
qlook_type = 'qlook_noise';
%Get the qlook files for each frame in selections
sid = 1; qlook_data = [];
%Get the proper directory
ql_dir = ct_filename_out(selections.param{sid},qlook_type);
for fid = 1:length(selections.frames{sid})
  %Generate the qlook file name
  ql_fn = fullfile(ql_dir,sprintf('Data_%s_%03.0f.mat',selections.seg_id{sid},...
    selections.frames{sid}(fid)));
  %Load the qlook file
  ql_now = load(ql_fn);
  %Pass to qlook variable
  if fid == 1
    qlook_data = ql_now;
  else
    keyboard
  end
end
%Get the x-axis values
xvals = geodetic_to_along_track(qlook_data(1).Latitude,qlook_data(1).Longitude);
%Get the y-axis values
yvals = qlook_data(1).Time;
%Plot the qlook data now
qlim = imagesc(xvals, yvals, lp(qlook_data(1).Data));
colormap(gca,flipud(gray))
hold on
%Set the button down function, tag, and userdata
qldata = struct('Latitude',qlook_data(1).Latitude,'Longitude',qlook_data(1).Longitude,...
  'GPS_time',qlook_data(1).GPS_time);
set(qlim,'ButtonDownFcn',@(s,e)checkclick(s,e),'Tag','qlook','UserData',qldata);
%Plot lines for the center of the gopro pictures
legs = {};
for rid = 1:length(refs)
  if strcmp(refs(rid).name, 'GoPro Coord.')
    for lid = 1:length(refs(rid).lat)
      %Find nearest record
      gonow = struct('lat',refs(rid).lat(lid),'lon',refs(rid).lon(lid));
      rec_id = near_record(gonow,qlook_data(1));
      %Plot the line showing the nearest click
      rplt = plot([xvals(rec_id) xvals(rec_id)], [0 1e-5],'--k','linewidth',2,'pickableparts','none','HandleVisibility','off');
      if lid == 1
        set(rplt,'HandleVisibility','on');
        legs{end+1} = refs(rid).name;
      end
    end
  end
end

%Plot the layers
[botlayer, surlayer] = loadlayers(selections);
plot(botlayer.alongtrack,botlayer.twtt,'--b','Linewidth',1.5)
legs{end+1} = 'Bottom Layer';
plot(surlayer.alongtrack,surlayer.twtt,'--m','Linewidth',2)
legs{end+1} = 'Surface Layer';
hold off
%Format plot
if fid > 1
  frmstr = sprintf('[%03.0f - %03.0f]', selections.frames{sid}(1), selections.frames{sid}(end));
else
  frmstr = sprintf('%03.0f',selections.frames{sid}(1));
end
legend(legs,'Fontsize',12,'Location','Best','Visible','on','Tag','qlook');
title(strrep(sprintf('%s: %s_%s',qlook_type,selections.seg_id{sid},frmstr),'_','\_'));
xlabel('Along Track (m)')
ylabel('Time (\mu s)')
ylim([min(qlook_data(1).Time) max(qlook_data(1).Time)])
%Reset axes tag
set(gca,'Tag',ql_tag);
end

%% Gopro plotting functions
function gopro_ref = plotgopro(recdata,selections)
%Get the gopro axes
if (~exist('selections','var') || isempty(selections)) || (~exist('recdata','var') || isempty(recdata))
  [recdata, ~, ~, ~, selections] = getfigdata;
end

%Get the qlook axes
go_tag = 'gopro';
go_ax = findall(0,'Type','axes','Tag',go_tag);
%Set to current axes
set(findall(0,'Type','figure','tag','View Frames'),'CurrentAxes',go_ax);
%Load the gopro master file
sid = 1;
gopro_data = load(ct_filename_support(selections.param{sid},'','gopro'));
%Load relevant gopro pics
if ~isfield(selections,'lat') || isempty(selections.lat)
  gopro_ids = 1;
else
  gopro_ids = near_record(selections,gopro_data);
end
%Make vector of gopro_ids
mos_width = 3; %number of images plotted
gopro_ids = gopro_ids-ceil(mos_width./2):2:gopro_ids+ceil(mos_width./2);

if min(gopro_ids) <=0
  gopro_ids = 1:2:2*mos_width;
elseif  max(gopro_ids) > length(gopro_data.lat)
  gopro_ids = fliplr(length(gopro_data.lat):-2:(length(gopro_data.lat)-mos_width*2+1));
end
%Make gopro mosaic
[lat_mos, lon_mos] = mosaic_plot(gopro_data,gopro_ids);
hold on
%Make records flightline on gopro plot
segtmp = recdata(strcmp(selections.seg_id,{recdata(:).id}));
sel_log_vec = (floor(segtmp.Frame_IDs) == selections.frames{:});
lon = segtmp.Longitude(sel_log_vec); lat = segtmp.Latitude(sel_log_vec);
if 1
  physical_constants
  %Load the layers
  [botlayer, surlayer] = loadlayers(selections);
  lon = botlayer.lon; lat = botlayer.lat;
  %Determine snow depth
  diffvec = abs(surlayer.twtt-botlayer.twtt);
  %Estimate snowdepth for this day
  depthvec = convlength(diffvec*c/(2*sqrt(er_ice)),'m','in'); %in
  %Do some filtering
  fillength = 40;
  while fillength > length(depthvec)
    fillength = fillength-4;
  end
  for did = 1:length(depthvec)
    windstrt = did-floor(fillength/2);
    windend = did+floor(fillength/2);
    if windstrt < 1
      windstrt = 1;
    end
    if windend>length(depthvec)
      windend = length(depthvec);
    end
    wind_ind = windstrt:windend;
    %Take average of window and log it
    depthvecs_filt(did) = mean(depthvec(wind_ind));
  end
  %Plot currently selected frame with depth colors
  
  depth = depthvecs_filt ; %// extract color index for the custom colormap
  lonmat=[lon' lon'];           %// create a 2D matrix based on "X" column
  latmat=[lat' lat'];           %// same for Y
  zz=zeros(size(lonmat)); %// everything in the Z=0 plane
  depthvec =[depth' depth'] ;         %// matrix for "CData"

  %// draw the surface (actually a line)
  hs=surf(lonmat,latmat,zz,depthvec,'EdgeColor','interp','FaceColor','none','Tag','gopro') ;

  colormap('spring') ;     %// assign the colormap
  shading flat                    %// so each line segment has a plain color
  view(2) %// view(0,90)          %// set view in X-Y plane
else
  %Plot just the records line
  plot(gca,lon,lat,'LineWidth',1.5,'Pickableparts','none','Color','r','Tag','gopro')
end

%Make gopro reference structure
gopro_ref = struct('lon',gopro_data.lon(gopro_ids),'lat',gopro_data.lat(gopro_ids),...
  'name','GoPro Coord.','mark',struct('type','s','color','k','size',10,'tag','marker'));
%Plot the gopro_ref locations
plot(gopro_ref.lon,gopro_ref.lat,'Marker',gopro_ref.mark.type,...
  'Color',gopro_ref.mark.color,'MarkerFaceColor',gopro_ref.mark.color,...
  'MarkerSize',gopro_ref.mark.size,'LineStyle','none')

hold off
xlim([min(min(lon_mos)) max(max(lon_mos))])
ylim([min(min(lat_mos)) max(max(lat_mos))])
%plot the picture
% gopro_fn = gopro_data.path{gopro_ids};
% if ispc
%   %use Z: drive
%   gopro_fn = strrep(gopro_fn,'/cresis/snfs1/data','Z:');
% end

% pic_data = imread(gopro_fn);
% goim = image(pic_data);

%format figue
xlabel('Longitude (^oE)')
ylabel('Latitude (^oN)')
title('Depth Map With Gopro')
set(findall(gca,'Type','Surface','Tag','gopro'),'Linewidth',3)

c = colorbar;
c.Label.String = 'Snow Depth (in.)';
c.FontSize = 12;
caxis([10 100])

%Reset axes tag
set(gca,'Tag',go_tag)
end

function [lat_mosaic, lon_mosaic, cdat_mosaic] = mosaic_plot(gopro_data,gopro_ids)
  %Get the inputs fields
  if iscell(gopro_data)
    gopro_data = gopro_data{1};
  end
  if ~exist('gopro_ids','var')
    gopro_ids = 1:length(gopro_data.lat);
  end
  gofields = fieldnames(gopro_data);
  gopro_elevs = gopro_data.elev(gopro_ids);
  gopro_h_agls = gopro_data.h_agl(gopro_ids);
  gopro_headings = gopro_data.heading(gopro_ids);
  gopro_rolls = gopro_data.roll(gopro_ids);
  gopro_pitchs = gopro_data.pitch(gopro_ids);
  gopro_lats = gopro_data.lat(gopro_ids);
  gopro_lons = gopro_data.lon(gopro_ids);
  
  pathfn = 'path';
  gopro_paths = {gopro_data(:).(pathfn){:}};
  gopro_paths = gopro_paths(gopro_ids);
  if ischar(gopro_paths) || length(gopro_paths)~=length(gopro_lats);
    %Grab the paths in a sneaky way
    gopro_cells = struct2cell(gopro_data);
    gopro_paths = reshape(gopro_cells(strcmp(pathfn,gofields),:,:),size(gopro_lats));
  end
  
  %Modulate radians and degrees for heading, rolls, and pitchs
  if max(abs(gopro_headings))<2*pi
    gopro_headings = rad2deg(gopro_headings);
  end
  if max(abs(gopro_rolls))<pi
    gopro_rolls = rad2deg(gopro_rolls);
  end
  if max(abs(gopro_pitchs))<pi
    gopro_pitchs = rad2deg(gopro_pitchs);
  end

  %Load gopro_images and gopro_infos based on gopro_paths
  gopro_images = {}; gopro_infos = {};
  for gid = 1:length(gopro_paths)
    %Load the images
    gopro_images{end+1} = imread(gopro_paths{gid});
    gopro_infos{end+1} = imfinfo(gopro_paths{gid});
  end

  lat_mosaic = []; lon_mosaic =[]; z_mosaic = []; cdat_mosaic = uint8([]);
  for mid = fliplr(1:length(gopro_images))
    latgrid = []; longrid = []; cdata = [];
    %debug
    if 0
      figure(mid)
      imagesc(rot90(gopro_images{mid},2))
      title('no changes')
    end
    %Get image state values
    lat = gopro_lats(mid); lon = gopro_lons(mid); heading = gopro_headings(mid);
    pitch_ang = gopro_pitchs(mid); roll_ang = gopro_rolls(mid);
    %Get the x and y resolutions
    ximg = convlength(gopro_infos{mid}.XResolution,'in','m'); %m
    yimg = convlength(gopro_infos{mid}.YResolution,'in','m'); %m

    %Get the focal length
    foc_leng = gopro_infos{mid}.DigitalCamera.FocalLength; %cm?
    fin35 = gopro_infos{mid}.DigitalCamera.FocalLengthIn35mmFilm;

    %Get height above ground
    H_agl_nadir = gopro_h_agls(mid); %m

    %Compensate for camangle
    cam_ang = 5;
  %   H_agl_nadir = H_agl_nadir/cosd(cam_ang);
    %Compensate for pitch
    H_agl_y = H_agl_nadir/cosd(pitch_ang);
    %Compensate for roll
    H_agl_x = H_agl_nadir/cosd(roll_ang);

    %Convert X and Y to ground distances
    xgrnd = ximg*H_agl_x/foc_leng; %m
    ygrnd = yimg*H_agl_y/foc_leng; %m

    %FOV
    x_fov = 122.6; %degrees
    y_fov = 94.4; %degrees

    %Get the image origin indexes and the
    img_rows = size(gopro_images{mid},1); img_cols = size(gopro_images{mid},2);
    center_row = floor(img_rows./2);
    center_col = floor(img_cols./2);
    top_rows = 1:center_row; bot_rows = center_row:img_rows;
    left_col = 1:center_col; rght_col = center_col:img_cols;

    %Make matrix for x-offsets and y-offsets (to the right and to the top are
    %positive distances [will represent distance north and distance east])
    topmaxdist = ygrnd*length(top_rows)./img_rows;
    topdistvec = fliplr(linspace(0,topmaxdist,length(top_rows)));
    botmaxdist = topmaxdist-ygrnd;
    botdistvec = linspace(0,botmaxdist,length(bot_rows));
    topfovvec = fliplr(linspace(0,y_fov/2,length(top_rows)));
    botfovvec = linspace(0,-y_fov/2,length(bot_rows));

    ydistvec = [topdistvec botdistvec(2:end)];
    yfovvec = [topfovvec botfovvec(2:end)];

    rghtmaxdist = xgrnd*length(rght_col)./img_cols;
    leftmaxdist = rghtmaxdist-xgrnd;
    leftdistvec = fliplr(linspace(0,leftmaxdist,length(left_col)));
    rghtdistvec = linspace(0,rghtmaxdist,length(rght_col));
    leftfovvec = fliplr(linspace(0,x_fov/2,length(left_col)));
    rghtfovvec = linspace(0,-x_fov/2,length(rght_col));

    xdistvec = [leftdistvec rghtdistvec(2:end)];
    xfovvec = [leftfovvec rghtfovvec(2:end)];

    %Adjust distvecs based on associated fov angle
    [xdistgrid, ydistgrid] = meshgrid(xdistvec,ydistvec);
    if 1
      [xfovgrid, yfovgrid] = meshgrid(xfovvec,yfovvec);
      xdistgrid = xdistgrid./cosd(xfovgrid);
      ydistgrid = ydistgrid./cosd(yfovgrid);
    else
      xfovgrid = zeros(size(xdistgrid));
      yfovgrid = zeros(size(ydistgrid));
    end

    %Compensate for the heading
    while abs(heading)>=360
      if heading < 0
        heading = heading + 360;
      else
        heading = heading - 360;
      end
    end
    %Do switching for different quadrants
    if (heading > 270 && heading < 360) || (heading < 0 && heading > -90)
      %Quadrant -x +y
      yn_cof = 1; xn_cof = 1;
      ye_cof = -1; xe_cof = 1;
    elseif (heading > 0 && heading <90) || (heading > -360 && heading < -270)
      %Quadrant +x +y
      yn_cof = 1; xn_cof = -1;
      ye_cof = 1; xe_cof = 1;
    elseif (heading > 90 && heading < 180) || (heading > -270 && heading < -180)
      %Quadrant +x -y
      yn_cof = 1; xn_cof = -1;
      ye_cof = 1; xe_cof = 1;
    else
      %Quadrant -x -y
      yn_cof = 1; xn_cof = 1;
      ye_cof = -1; xe_cof = 1;
    end

    northdistgrid = yn_cof*ydistgrid.*cosd(-heading) + xn_cof*xdistgrid.*sind(-heading);
    eastdistgrid = ye_cof*ydistgrid.*sind(-heading) + xe_cof*xdistgrid.*cosd(-heading);

    %Estimate latitude grid values
    Re_pol = 6357e3; Re_eq = 6378e3; %m Earth Radius at poles and equator
    Ce_pol = Re_pol*2*pi; Ce_eq = Re_eq*2*pi; %m Earth Circumference
    east_m2lat_deg = Ce_pol./360; %m/deg conversion from meters to degrees
    latoffgrid = eastdistgrid./east_m2lat_deg; %degrees of latitude offset
    latgrid = lat+latoffgrid;

    %Estimate longitude grid values
    north_m2lon_deg = Ce_eq./360*cosd(lat);
    lonoffgrid = northdistgrid./north_m2lon_deg;
    longrid = lon+lonoffgrid;

    %Grab cdata
    cdata = rot90(gopro_images{mid},2);
    %Get inner image based on FOV limit
    fovlim = 30;
    if 1
      fovlogmat = ones(size(latgrid));
      xfovlog = xfovgrid<=fovlim;
      yfovlog = yfovgrid<=fovlim;
      fovlogmat = (xfovlog+yfovlog)==2;
      latgrid(~fovlogmat) = nan;
      longrid(~fovlogmat) = nan;
      clogmat(:,:,1) = fovlogmat;
      clogmat(:,:,2) = fovlogmat;
      clogmat(:,:,3) = fovlogmat;
      cdata(~clogmat) = nan;
    end

    %Reduce image fidelity
    if 1
      subrows = floor(linspace(1,size(latgrid,1),500));
      subcols = floor(linspace(1,size(latgrid,2),500));
      latgrid = latgrid(subrows,subcols);
      longrid = longrid(subrows,subcols);
      cdata = cdata(subrows,subcols,:);
    end
    zgrid = zeros(size(latgrid));

    %Send to mosaic variables for surf plotting
    lat_mosaic(end+1:end+size(latgrid,1),:) = latgrid;
    lon_mosaic(end+1:end+size(latgrid,1),:) = longrid;
    z_mosaic(end+1:end+size(latgrid,1),:) = zgrid;
    cdat_mosaic(end+1:end+size(latgrid,1),:,:) = cdata;
    %Add nan column for separating the different images
    nanvec = nan(1,size(latgrid,2)); nanmat = nan([1 size(latgrid,2) size(cdat_mosaic,3)]);
    lat_mosaic(end+1,:) = nanvec;
    lon_mosaic(end+1,:) = nanvec;
    z_mosaic(end+1,:) = nanvec;
    cdat_mosaic(end+1,:,:) = nanmat;

    if 0
      figure(100)
      hold on
      surf(longrid,latgrid,zgrid,cdata,'edgecolor','none')
      hold off
      set(gca,'YDir','normal')
      view(2)
      keyboard
    end
  end

  %Plot the completed mosaic
  geoshow(lat_mosaic,lon_mosaic,uint8(cdat_mosaic))
  grid on
  xlabel(sprintf('Longitude ^o N'))
  ylabel(sprintf('Latitude ^o E'))
end

%% Clicking effects
function [] = checkclick(src, evt)
%Get source type
srctype = get(src,'type');
%Get the lat lon of the click
% Switch between selection functions
if strcmp(srctype, 'line')
  if ~isempty(get(src,'UserData'))
    clicksegment(src,evt)
  end
elseif strcmp(srctype, 'image') && strcmp(src.Tag, 'qlook')
  clickqlook(src,evt)
elseif strcmp(srctype, 'image')
  %Warn user to click the line
  warning('Nothing happened... Please select either the red window or the records lines')
else
  keyboard
end
end

function [] = clickqlook(src,evt)
[recdata, start, stop, refs, selections] = getfigdata;
%Grab current selection data
frames_select = selections.frames{1};
seg_select = selections.seg_id{1};
param_select = selections.param{1};
%Get the x index of the clicked value
xclk = evt.IntersectionPoint(1);
xdat = src.XData;
rec_id = find(xdat>=xclk,1);
%Get values from src.UserData
sel_lat = src.UserData.Latitude(rec_id); sel_lon = src.UserData.Longitude(rec_id);
sel_gpstime = src.UserData.GPS_time(rec_id);
%Load selections structure
selections.seg_id = {seg_select};
selections.param = {param_select};
selections.frames = {frames_select};
selections.lat = sel_lat;
selections.lon = sel_lon;
selections.gpstime = sel_gpstime;
%Save selections
savefigdata(recdata, start, stop, refs, selections);
%Update the frames
updateframes(src,evt)
end

function [] = clicksegment(src,evt)
[recdata, start, stop, refs, selections] = getfigdata;
%Check the button click
ebut = evt.Button;
%Determine closest frame and all displayed frames
chk_dat = struct('Longitude',src.XData,'Latitude',src.YData);
chker = struct('Longitude', evt.IntersectionPoint(1), 'Latitude', evt.IntersectionPoint(2));
rec_id = near_record(chker,chk_dat);
%Get values from src.UserData
seg_frm_ids = unique(src.UserData.Frame_IDs);
pickfrm_id = src.UserData.Frame_IDs(rec_id);
seg_select = src.UserData.id; param_select = src.UserData.param;
sel_lat = src.UserData.Latitude(rec_id); sel_lon = src.UserData.Longitude(rec_id);
sel_gpstime = src.UserData.GPS_time(rec_id);
%Check if select entire segment or one frame
perframe = true;
idchk = {};
if perframe
  frames_select = pickfrm_id;
else
  frames_select = seg_frm_ids;
end
%Make full id for warnings
frm_str = '';
for f = frames_select
  frm_str = sprintf('%s%s, ',frm_str,num2str(f));
end
full_id = sprintf('%s_[%s]',seg_select,frm_str(1:end-2));
%Load selections structure
selections.seg_id = {seg_select};
selections.param = {param_select};
selections.frames = {frames_select};
selections.lat = sel_lat;
selections.lon = sel_lon;
selections.gpstime = sel_gpstime;
%Save selections
savefigdata(recdata, start, stop, refs, selections);
%Update the frames
updateframes(src,evt)
end

%% Update functions
function updateframes(s, e)
%Change pointer to waiting
vwfrm_fig = findall(0,'type','figure','tag','View Frames');
oldpoint = waitpoint(vwfrm_fig);
% Load the current data
[recdata, start, stop, refs, selections] = getfigdata;
%Plot the gopro (should be first)
gopro_ref = plotgopro(recdata,selections);
if ~isempty(gopro_ref) && strcmp(refs(end).name, gopro_ref.name)
  refs(end) = gopro_ref;
else
  refs(end+1) = gopro_ref;
end
%Plot the qlook (after gopro)
plotqlook(refs,selections);
%Plot the frames (must be last)
  %Delete any lines on the map
  delete(findall(findall(gcf,'Type','Axes','Tag','map'),'Type','Line'))
plotframes(recdata,start,stop,refs, selections);
%Save data
savefigdata(recdata,start,stop,refs,selections);
%Change pointer back
waitpoint(vwfrm_fig,oldpoint);
end

%% Figure functions (close)
function chksettings(minimize_bool,open_bool)
%This function determines if the settings window is open (opens if not)
jwarnid = 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame';
warning('off',jwarnid)
set_fig = findall(0,'type','figure','tag','Setup');
if strcmp(get(gcf,'Tag'),'View Frames')
  vwfrm_fig = findall(0,'type','figure','tag','View Frames');
else
  vwfrm_fig = gcf;
end
if isempty(set_fig) && (~exist('open_bool','var') || open_bool)
  recdata = getfigdata();
  view_setup([],recdata);
  set_fig = findall(0,'type','figure','tag','Setup');
end
%Maximize or minimize
if ~isempty(set_fig)
  if nargin > 0
    if minimize_bool
      %Minimize it
      if ~set_fig.JavaFrame.isMinimized
        set_fig.JavaFrame.setMinimized(1);
      end
    else
      %Unminimize it
      if set_fig.JavaFrame.isMinimized
        set_fig.JavaFrame.setMinimized(0);
      end
    end
  else
    %switch it
    if ~set_fig.JavaFrame.isMinimized
      set_fig.JavaFrame.setMinimized(1);
    else
      set_fig.JavaFrame.setMinimized(0);
    end
  end
end

      
end

function close_view(s,e)
%Determine if initial settings window is open
chksettings(false,true)
%Close this figure
delete(findall(0,'type','figure','tag','View Frames'))
end

%% Parsing functions
function [closeid, Xdata, Ydata] = near_id(src,evt)
cornx = evt.IntersectionPoint(1); corny = evt.IntersectionPoint(2);
%Get the x and y data
Xdata = get(src,'Xdata'); Ydata = get(src,'Ydata');
%Get the closest corner
distdat = sqrt((cornx-Xdata).^2+(corny-Ydata).^2);
[~, closeid] = min(distdat);
end

function [nearid] = near_record(chk_struct,dat_struct)
%Grab the check values
if isfield(chk_struct,'Latitude') && isfield(chk_struct,'Longitude')
  %for qlook data
  latchk = chk_struct.Latitude; lonchk = chk_struct.Longitude;
elseif isfield(chk_struct,'lat') && isfield(chk_struct,'lon')
  %for gopro data
  latchk = chk_struct.lat; lonchk = chk_struct.lon;
else
  nearid = 1;
  return
end
%Grab the data values
if isfield(dat_struct,'Latitude') && isfield(dat_struct,'Longitude')
  %for qlook data
  lattmp = dat_struct.Latitude; lontmp = dat_struct.Longitude;
elseif isfield(dat_struct,'lat') && isfield(dat_struct,'lon')
  %for gopro data
  lattmp = dat_struct.lat; lontmp = dat_struct.lon;
else
  nearid = 1;
  return
end
%Get the nearest record
[~, nearid] = min(sqrt((latchk-lattmp).^2+(lonchk-lontmp).^2));
end