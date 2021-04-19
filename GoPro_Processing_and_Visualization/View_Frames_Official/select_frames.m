% script plot flightlines
%
% This script uses data loaded through run_gopro_plot to display the
% flightlines that intercept the geographic window defined in check_region.
% 
% Plot controls include:
%   Select Flightline:      Left Click
%   De-Select Flightline:   Right Click
%   Move Window:            Left Click Window Corner Markers 
%                               Then Click Again in New Position
%   Print Annotation Data:  Left Click Annotation
%   Clear Annotation:       Right Click Annotation
%
% Author: Bailey Miller
%
% See Also: run_gopro_plot.m

function [selections] = select_frames(src,evt,recdata,start,stop,refs,selec_fcn)
fignum = 100;
%Load guidata
if ~isempty(findall(0,'Type','Figure','Number',fignum))
  set(0,'currentfigure',fignum)
  data = guidata(gcf);
else
   %% Format Figure
  %Settings
  fig_pos = [0 0 1 1];
  figure(fignum)
  set(gcf,'Units','Normalized','outerposition',fig_pos,'Name','Select_Frames','Numbertitle','off')
  data = [];
end

%% Input Check
if nargin == 1
  axchng = src;
else
  axchng = false;
end
if ~exist('recdata','var') || isempty(recdata)
  udata = get(gcf,'UserData');
  if ~isfield(udata,'param')
    try
      load('chktmp.mat')
    catch
      error('Temp file does not exist. Please run check_region before proceeding.')
    end
  else
    recdata = get(gcf,'UserData');
  end
else
  set(gcf,'UserData',recdata);
end

if ~exist('start','var') 
  if ~isfield(data,'start')
    start.lat = []; start.lon = [];
  else
    start = data.start;
  end
end

if ~exist('stop','var') 
  if ~isfield(data,'stop')
    stop.lat = []; stop.lon = [];
  else
    stop = data.stop;
  end
end

if ~exist('refs','var') 
  if ~isfield(data,'refs')
    refs.lat = []; refs.lon = [];
  else
    refs = data.refs;
  end
end


%% Plot the different windows
  ax_pos = [0.0643    0.1100    0.6477    0.8150];
  sel_strt = ax_pos(1)+ax_pos(3)+.01;
  sel_pos = [sel_strt .2 .99-sel_strt .75];
  butt_pos = [sel_strt .05 .99-sel_strt .15];
%Load guidata
if isempty(data)
  %Initialize the framescheck variable
  selections.seg_id = cell(1); selections.param = cell(1); selections.frames = cell(1);
  data = struct('selections',selections,'refs',refs, 'start', start, 'stop', stop);
else
  selections = data.selections;
end
guidata(gcf,data);

%Plot the map make axes settings
mapim = findall(gca,'Type','Image','Tag','map');
if isfield(recdata(1).param,'gisfn') && isempty(mapim)
  gisdata = geotiff_read(recdata(1).param.gisfn);
  mapim = imagesc(gisdata.X, gisdata.Y, gisdata.C);
  set(mapim,'PickableParts','none','Tag','map')
  hold('on'); grid('on');
  %Adjust axes settings
  set(gca,'YDir','normal','Units','Normalized','Position',ax_pos)
end

%Post segment selector
selector_tag = 'Selector_Panel';
if isempty(findall(gcf,'Tag',selector_tag))
  make_selector(recdata,sel_pos,selector_tag);
end
%Post select all button
butts_tag = 'Buttons_Panel';
if isempty(findall(gcf,'Tag',butts_tag))
  make_buttons(butt_pos,butts_tag,selec_fcn)
end
%% Plot the segments
%Delete any lines
delete(findall(gca,'Type','Line'))
%Plot the frame
[axlim, lonvec, latvec] = plotframes(recdata,start,stop,selections);
hold on
%Make the window
plotwind(lonvec,latvec)
%Plot Reference location
plotrefs(refs,axlim)
%Post annotations
win_annot(selections,start,stop)
%Format
if axchng
  axis(axlim)
end
xlabel('Longitude (^o E or W)')
ylabel('Latitude (^o N or S)')
% Make the title
updatetitle
hold off
%Move map down
uistack(mapim,'bottom')
drawnow
end
%% Panel build functions
function [] = make_selector(recdata,sel_pos,sel_tag)
%Make selection panel for displaying lines
main_pan = uipanel(gcf,'Units','Normalized','Position',sel_pos,...
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
  pos_now = [.01 (1/length(unq_seasons)*(uid-1)+.025) .98 1/length(unq_seasons)-.05];
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
for uid = 1:length(season_strs)
  seg_strs = sort(season_strs{uid});
  for cid = 1:length(seg_strs)
    posnow = [.025 (.96-1/length(seg_strs)*(cid-1))  1 .05];
    cbnow = uicontrol(unq_pans(uid),'Style','checkbox','String',seg_strs{cid},'Units','Normalized',...
      'Position',posnow,'Value',0,'Callback',@(src,evt)select_frames(src,evt),'Tag',seg_strs{cid});
  end
  %Make hide all button
  uicontrol(unq_pans(uid),'Units','Normalized','Position',[.7 .8 .2 .1],...
    'String','Hide All','Callback',@(src,evt)hide_all_segs(src,evt))
end
end


function [] = make_buttons(pos,tag,selec_fcn)
but_pan = uipanel(gcf,'Units','Normalized','Position',pos,'Tag',tag);
%Make select all button
uicontrol(but_pan,'Units','Normalized','String','Select All Shown',...
  'Callback',@(src,evt)but_select_all(src,evt),'Position',[.05 .55 .9 .4])
%Make gopro plot button
uicontrol(but_pan,'Units','Normalized','String','Plot Selections...',...
  'Callback',selec_fcn,'Position',[.05 .05 .9 .4])
end

function [] = hide_all_segs(src,evt)
%Get all the checkboxes
chkbxs = findall(get(src,'Parent'),'Style','checkbox');
if strcmp('Hide All',src.String)
  makeval = 1; makestr = 'Show All';
else
  makeval = 0; makestr = 'Hide All';
end
set(src,'String',makestr);
set(chkbxs,'Value',makeval);
select_frames;
end

function [] = but_select_all(src,evt)
%Get guidata
data = guidata(gcf);
selections = data.selections;
%Get the shown segments and add to select now
seg_lines = findall(gcf,'Type','Line','Tag','Segment');

for s_id = 1:length(seg_lines)
  src = seg_lines(s_id);
  frames_select = unique(src.UserData.Frame_IDs);
  seg_select = src.UserData.id; param_select = src.UserData.param.season_name;
  %Make full id for warnings
  frm_str = '';
  for f = frames_select
    frm_str = sprintf('%s%s, ',frm_str,num2str(f));
  end
  %Determine if segment already selected
  if isfield(selections,'param') && ~isempty(selections.param{1})
    seg_log_vec = cellfun(@(c)strcmp(c,seg_select),selections.seg_id);
    if any(seg_log_vec) %== Already selected segment
      %Determine if frame already selected
      frms_vec = selections.frames{seg_log_vec};
      frm_log_vec = zeros(1,length(frms_vec));
      for frm_now = frames_select
        frm_log_vec = frm_log_vec | frms_vec == frm_now;
      end
      % Add all newly selected frms
      frms_vec = sort(unique([frms_vec frames_select]));
      selections.frames{seg_log_vec} = frms_vec;

    else %== Newly selected segment
      %Add new selected segment and frames
      selections.seg_id{end+1} = seg_select;
      selections.frames{end+1} = frames_select;
      selections.param{end+1} = param_select;
    end
  else %No lines shown
    selections.seg_id = {seg_select};
    selections.frames = {frames_select};
    selections.param= {param_select};
  end
end
% Output framescheck
data.selections = selections;
guidata(gcf,data);
% Update Figure
select_frames;
end

%% Plotting functions
function [] = plotwind(lonvec,latvec)
%Delete previous window
delete(findall(gca,'Type','Line','Tag','geowin'))
%Plot new window
gwin = plot(gca,lonvec,latvec,'--sr','MarkerFaceColor','r','MarkerSize',20,'Tag','geowin');
set(gwin,'ButtonDownFcn',@(s,e)checkclick(s,e))
set(gwin,'HandleVisibility','off')
end

function [axlim, lonvec, latvec] = plotframes(recdata,start,stop,selections)
%Input check
if ~exist('recdata','var')
  recdata = get(gcf,'UserData');
end
if ~exist('start','var') || ~exist('stop','var') || ~exist('selections','var')
  gdata = guidata(gcf);
  start = gdata.start; stop = gdata.stop; selections = gdata.selections;
end
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
param_season = cell(length(recdata),1); legnow = {};
colors2choose = {[0 0 255],  [0 255 0], [255 255 0], [7 252 255]}; chid =2;
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
  seg_chkbx = findall(gcf,'Style','checkbox','Tag',recdata(recid).id);
  if any(frm_log_vec) && ~seg_chkbx.Value
    %Reset chkbx text
    set(seg_chkbx,'String',recdata(recid).id)
    %Modulate segtmp based on logical vector
    plttmp = struct('Longitude',segtmp.Longitude(frm_log_vec),'Latitude',segtmp.Latitude(frm_log_vec),'Elevation',segtmp.Elevation(frm_log_vec),...
      'GPS_time',segtmp.GPS_time(frm_log_vec),'Frame_IDs',floor(segtmp.Frame_IDs(frm_log_vec)),'param',segtmp.param,'filename',segtmp.filename,'id',segtmp.id);

    %Plot each segment
    gn = plot(gca,plttmp.Longitude,plttmp.Latitude,'.');
    set(gn,'ButtonDownFcn',@(s,e)checkclick(s,e),'UserData',plttmp,...
      'Tag','Segment','LineWidth',2,'HandleVisibility','off')
    hold on
    
    %Check for adding legend entry and season color
    if isempty(legnow) || ~any(cellfun(@(c)strcmp(param_season{recid},c),legnow))
      set(gn,'HandleVisibility','on');
      legnow{end+1} = param_season{recid};
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
%Make the legend
for lid = 1:length(legnow)
  legnow{lid} = strrep(legnow{lid},'_',' ');
end
if isempty(legnow)
  wid = 'MATLAB:legend:PlotEmpty';
  warning('off',wid)
end
legend(legnow,'Visible','off')
%Post the frames annotation
frames_annot
end

function [] = plotrefs(refs,axlim)
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
  %Insert ref in legend
  legnow = legend;
  if isempty(legnow)
    legnow = legend(ref_leg);
  else
    legnow_str = legnow.String;
    for r_id = 1:length(ref_leg)
      legnow_str{end+1} = ref_leg{r_id};
    end
    [legnow,legobjs] = legend(legnow_str);
    %Adjust children text size
    ch = [];
    for l_id = 1:length(legnow.String(1:end-length(ref_leg)))
      chnow = findobj(legobjs,'Type','line','Tag',legnow.String{l_id});
      ch = [ch chnow];
    end
    %Adjust param marker size
    set(ch,'Markersize',35);
  end
end
set(findall(gcf,'Type','legend'),'Fontsize',12,'Visible','on')
end

%% Annotation Functions

function [] = frames_annot()
%Position settings
annot_pos = [0.0108 0.2460 0.1578 0.1649];
%Load the guidata
gdata = guidata(gcf);
%Get the num of selected frames
numsel = length(findall(gca,'Type','Line','Tag','Select'));
%Load the frames
selections = gdata.selections;
%Make the pretty string cell array
if ~isempty(selections) && ~isempty(selections.seg_id{1})
  %Print the unique params
  unqparam = unique(selections.param);
  anstr = {sprintf('Selected Frames[%d](%d):     {Right-Click to Clear All}',length(unqparam),numsel)};
  %Print by param
  for uid = 1:length(unqparam)
    anstr{end+1}=sprintf('\t%s:',unqparam{uid});
    %Print by seg_id
    for sid = 1:length(selections.seg_id)
      %Confirm the param
      if strcmp(unqparam{uid},selections.param{sid})
        anstr{end+1} = sprintf('\t\t%s',selections.seg_id{sid});
        %Add the selected frames
        fstr = frameids2str(selections.frames{sid});
        anstr{end+1} = sprintf('\t\t\t%s',fstr);
      end
    end
  end
end
%Look for the annotation
a = findall(gcf,'Type','textbox','Tag','Frames');
if isempty(selections.seg_id) || isempty(selections.seg_id{1})
  %Delete the annotation
  delete(a)
elseif isempty(a)
  %Plot the annotation
  annotation('textbox', annot_pos, 'String', anstr, ...
    'BackgroundColor','w','interpreter','none','UserData',selections,...
    'ButtonDownFcn',@(s,e)checkclick(s,e),'Tag','Frames','FitBoxToText','on')
else
  %Update the annotation
  set(a,'String',anstr,'FitBoxToText','on')
end
end

function [] = win_annot(selections,start,stop)
%Settings
annot_pos = [0.0141 0.8425  0.1728 0.1524];
%Make the pretty string cell array
  winstr = {'Selected Window:'};
  winstr{end+1} = sprintf('Start:\n\tLat = %.2f^o\tLon = %.2f^o',start.lat,start.lon);
  winstr{end+1} = sprintf('Stop:\n\tLat = %.2f^o\tLon = %.2f^o',stop.lat,stop.lon);
  %Look for the annotation
  a = findall(gcf,'Type','textbox','Tag','Window');
  if isempty(selections)
    %Delete the annotation
    delete(a)
  elseif isempty(a)
    %Plot the annotation
    annotation('textbox', annot_pos, 'String', winstr, ...
      'BackgroundColor','w','UserData',selections,'FitBoxToText','on',...
      'ButtonDownFcn',@(s,e)checkclick(s,e),'Tag','Window')
  else
    %Update the annotation
    set(a,'String',winstr,'FitBoxToText','on')
  end
end


%% Clicking effects
function [] = checkclick(src, evt)
%Get source type
srctype = get(src,'type');
% Switch between selection functions
if strcmp(srctype, 'line')
  if ~isempty(get(src,'UserData'))
    clicksegment(src,evt)
  else
    [start,stop] = clickwindow(src,evt);
    %Update the window
    updatewindow(src,evt,start,stop)
  end
elseif strcmp(srctype, 'textboxshape')
  clickannot(src,evt)
elseif strcmp(srctype, 'image')
  %Warn user to click the line
  warning('Nothing happened... Please select either the red window or the records lines')
else
  keyboard
end
%Update the title
updatetitle
end

function [] = clicksegment(src,evt)
data = guidata(gcf);
selections = data.selections;
%Check the button click
ebut = evt.Button;
%Determine closest frame and all displayed frames
rec_id = near_id(src,evt);
seg_frm_ids = unique(src.UserData.Frame_IDs);
pickfrm_id = src.UserData.Frame_IDs(rec_id);
seg_select = src.UserData.id; param_select = src.UserData.param.season_name;
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
%Determine if segment already selected
if isfield(selections,'param') && ~isempty(selections.param{1})
  seg_log_vec = cellfun(@(c)strcmp(c,seg_select),selections.seg_id);
  if any(seg_log_vec) %== Already selected segment
    %Determine if frame already selected
    frms_vec = selections.frames{seg_log_vec};
    frm_log_vec = zeros(1,length(frms_vec));
    for frm_now = frames_select
      frm_log_vec = frm_log_vec | frms_vec == frm_now;
    end
    %Modulate based on button press
    if ebut == 1
      if sum(frm_log_vec) == length(frames_select) %== All frms selected
        warning('Frame %s has already been selected.',full_id);
        return
      else
        % Add all newly selected frms
        frms_vec = sort(unique([frms_vec frames_select]));
        selections.frames{seg_log_vec} = frms_vec;
      end
    elseif ebut == 3
      if all(~frm_log_vec) %== Frame not selected
        warning('Frame %s is not selected.',full_id);
        return
      else
        %Deselect frm
        frms_vec = frms_vec(~frm_log_vec);
        if isempty(frms_vec) || ~perframe %Deselect entire segment
          if sum(~seg_log_vec) == 0 %Everything is deleted
            newsegs = []; newfrms = []; newprms = [];
          else
            newsegs = selections.seg_id(~seg_log_vec);
            newfrms = selections.frames(~seg_log_vec);
            newprms = selections.param(~seg_log_vec);
          end
          if isempty(newsegs) || ischar(newsegs)
            %Reinitialize
            selections.seg_id = {newsegs};
            selections.frames = {newfrms};
            selections.param = {newprms};
          elseif iscell(newsegs)
            selections.seg_id = newsegs;
            selections.frames = newfrms;
            selections.param = newprms;
          end
        else %Deselect single frame
          selections.frames{seg_log_vec} = sort(unique(frms_vec));
        end          
      end
    else
      warning('Line not eligible for this operation.')
      return
    end
        
  else %== Newly selected segment
    if ebut ==1
      %Add new selected segment and frames
      selections.seg_id{end+1} = seg_select;
      selections.frames{end+1} = frames_select;
      selections.param{end+1} = param_select;
    elseif ebut == 3
      warning('Frame %s is not selected.',full_id);
      return
    else
      warning('Line not eligible for this operation.')
    return
    end
  end
else %No current selections
  if ebut == 1
    selections.seg_id = {seg_select};
    selections.frames = {frames_select};
    selections.param= {param_select};
  elseif ebut == 3
    warning('Frame %s is not selected.',full_id);
    return
  else
    warning('Line not eligible for this operation.')
    return
  end
end
% Output framescheck
data.selections = selections;
guidata(gcf,data);
% Update Figure
select_frames;
end

function [start, stop] = clickwindow(src,evt)
gdata = guidata(gcf);
%Find closest index
[datselec_id, Xdata, Ydata] = near_id(src,evt);
%Plot the corner selected
hold on
gs = plot(Xdata(datselec_id),Ydata(datselec_id),'sk',...
  'MarkerSize',16,'MarkerFaceColor','k');
hold off
%Tell the user to click at the new window location
t = title('Select the new corner location:');
%Make frames unpickable
set(findall(gca,'Type','Line'),'Pickableparts','none')
mapim = findall(gca,'Type','Image');
set(mapim,'Pickableparts','all')
%Do the selection
[newx, newy] = ginput(1);
%Delete unecessary elements
delete(gs)
delete(t)
%Fix the lines and image
set(findall(gca,'Type','Line'),'Pickableparts','all')
set(mapim,'Pickableparts','none')
%Get the unique lat and lon values
lonvec = unique(Xdata(Xdata~=Xdata(datselec_id)));
latvec = unique(Ydata(Ydata~=Ydata(datselec_id)));
%Output new start and stop
start.lon = lonvec; start.lat = latvec;
stop.lon = newx; stop.lat = newy;
%Save guidata
gdata.start = start; gdata.stop = stop;
guidata(gcf, gdata);
end

function [] = clickannot(src,evt)
%Get all the annotations
a = findall(gcf,'Type','textbox','Tag',src.Tag);
if evt.Button == 1 
  %Print the text in the current annotation to the command window
  clc
  for i = 1:length(src.String)
    disp(src.String{i})
  end
elseif evt.Button == 3
  %Delete annotation and reset frames
  delete(a)
  %Get current guidata
  data = guidata(gcf);
  selections.seg_id = cell(1); selections.frames = cell(1); selections.param = cell(1);
  data.selections = selections;
  guidata(gcf,data);
  %Change line color and tag
  set(findobj(gcf,'Type','Line','Tag','Select'),'Color','b','Tag','No');
end
%Update figure
select_frames(src,evt);
end

%% Update functions
function [closeid, Xdata, Ydata] = near_id(src,evt)
cornx = evt.IntersectionPoint(1); corny = evt.IntersectionPoint(2);
%Get the x and y data
Xdata = get(src,'Xdata'); Ydata = get(src,'Ydata');
%Get the closest corner
distdat = sqrt((cornx-Xdata).^2+(corny-Ydata).^2);
[~, closeid] = min(distdat);
end

function [] = updatewindow(src,evt,start,stop)
  %Get frames check
  frames = get(gcf,'UserData');
  %Window annotation
  win_annot(frames,start,stop)
  %Change the axis size
  axchng = true;
  %Replot everything if the window changes
  select_frames(axchng);
end

function [] = updatetitle()
% Make title string
ldata = get(findall(gca,'Type','Line','Tag','Segment'),'UserData');
if ~iscell(ldata)
  ldata = {ldata};
end
if ~isempty(ldata{1})
  numfrms = 0;
  for lid = 1:length(ldata)
    linnow = ldata{lid};
    numfrms = numfrms + length(unique(linnow.Frame_IDs));
  end
else
  numfrms = 0;
end
numsel = length(findobj(allchild(gca),'Type','line','Tag','Select'));
titlstr = {sprintf('Frames Shown: %d',numfrms), sprintf('Frames Selected: %d',numsel)};
title(titlstr)
end