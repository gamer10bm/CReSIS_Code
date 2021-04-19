function app = view_setup(app, recdata, start, stop, refs)
% app = view_setup(app, params, start, stop, refs)
% 
% This function generates the setup figure for the run_view_frames function
%


%Check inputs
if nargin == 4
  %% First opening of the funcxtion
  refs = stop; stop = start; start = recdata; recdata = app;
  app = [];
  %Initialize the setup figure
  set = init_setup(recdata,start,stop,refs);
elseif nargin == 2
  %% Callback
  %Initialize the setup figure
  set = init_setup(recdata);
end
app.set = set;
%Save to guidata
guidata(set.fig,set)

%Check if view frames is already open
vwfrm_fig = findall(0,'type','figure','Tag','View Frames');
if ~isempty(vwfrm_fig)
  vwfrm_dat = guidata(vwfrm_fig);
  if isfield(vwfrm_dat,'recdata')
    %Determine the segment ids and frame ids
    seg_ids = recdata2segids(vwfrm_dat.recdata);
    %Make segchk_selec call
    for sid = 1:length(seg_ids)
      frms = [];
      frms = unique(floor(vwfrm_dat.recdata(sid).Frame_IDs));
      segchk_selec([],[],seg_ids{sid}, frms);
    end
  end
end
end

function set = init_setup(recdata,start,stop,refs)
%% Input Checking
if ~exist('start','var') || isempty(start)
  [start, stopper] = recdata2startstop(recdata);
end
if ~exist('stop','var') || isempty(stop)
  stop = stopper;
end
if ~exist('refs','var')
  refs = [];
end
%Make the base figure
fig_pos = [97    99   607   597];
main_pos = [0 0 .3 1];
selec_pos = [.3 0 .7 1];
set.fig = figure('Position',fig_pos,'Name','Setup View Frames',...
  'Numbertitle','off','Tag','Setup');

%% Setup the main panel and selected frames panel
set.mainp = uipanel('Parent',set.fig,'Title','Setup Window',...
  'FontSize',14,'Background','white','Units','Normalized','Position',main_pos);
set.selecp = uipanel('Parent',set.fig,'Title','Selection Window',...
  'FontSize',14,'Background','white','Units','Normalized','Position',selec_pos);

  %% Make selection buttons
  selbut_pos = [0 .2 1 .8];
  set.seg_selec = uipanel('Parent',set.mainp,'Units','Normalized','Position',selbut_pos,...
    'Title','Select Segments...','FontSize', 12);
  %Select all segments
  all_pos = [.05 .55 .9 .4];
  uicontrol(set.seg_selec,'Style','Pushbutton','String','Select All','Units','Normalized',...
    'Callback',@(s,e)select_all(s,e),'Position',all_pos);
  %Select frames with map
  map_pos = [.05 .05 .9 .4];
  selec_fcn = @(src,evt)select_plot(src,evt);
  uicontrol(set.seg_selec,'Style','Pushbutton','String','Use Map...','Units','Normalized',...
    'Callback',@(s,e)select_frames(s,e,recdata,start,stop,refs,selec_fcn),'Position',map_pos,'Tag','Setup');
    
  %% Make end/ enter buttons
  endbut_pos = [ .05 .025 .4 .15];
  end_fcn = @(han,event)close(set.fig);
  start_out.end_but = uicontrol(set.mainp, 'Style', 'pushbutton', ...
    'String', 'Quit','Units','normalized','Position',endbut_pos,'Callback',end_fcn);

  strtbut_pos = [ .55 .025 .4 .15];
  plot_fcn = @(han, event)plot_start(han, event,start,stop,refs);
  start_out.strt_but = uicontrol(set.mainp, 'Style', 'pushbutton', ...
    'String', 'Plot...','Units','normalized','Position',strtbut_pos,'Callback',plot_fcn);
  
  %% Make Selection Panel
  %Make tabgroup
  tabgrp = uitabgroup(set.selecp);
  %Make a mini panel for each unique season
  param_season = cell(length(recdata),1);
  for r_id = 1:length(recdata)
    param_season{r_id} = recdata(r_id).param.season_name;
  end
  unq_seasons = unique(param_season);
  unq_pans = []; season_strs = cell(1,length(unq_seasons));
  season_frms = cell(1,length(unq_seasons));
  for uid = 1:length(unq_seasons)
    pos_now = [.01 (1/length(unq_seasons)*(uid-1)+.025) .98 1/length(unq_seasons)-.05];
    u = uitab(tabgrp,'Units','Normalized','Title',strrep(unq_seasons{uid},'_',' '));
    unq_pans = [unq_pans u];
    %Find segments for each season
    season_strs{uid} = {};
    season_frms{uid} = {};
    for r_id = 1:length(recdata)

      if strcmp(recdata(r_id).param.season_name,unq_seasons{uid})
        season_strs{uid}{r_id} = recdata(r_id).id;
        season_frms{uid}{r_id} = recdata(r_id).Frame_IDs;
      end
    end
  end
  %Make checkboxes for each season
  numcols = 2;
  for uid = 1:length(season_strs)
    [seg_strs, seg_logvec] = sort(season_strs{uid});
    seg_frms = season_frms{uid}(seg_logvec);
    for colid = 1:numcols
      seg_str_now = seg_strs(colid:numcols:end);
      seg_frm_now = seg_frms(colid:numcols:end);
      for cid = 1:length(seg_str_now)
        %Make unique position
        posnow = [(.025+1/numcols*(colid-1)) (.95-1/length(seg_str_now)*(cid-1))  (1/numcols) .05];
        %Make checkbox
        cbnow = uicontrol(unq_pans(uid),'Style','checkbox','String',seg_str_now{cid},'Units','Normalized',...
          'Position',posnow,'Value',0,'Callback',@(s,e)segchk_selec(s,e),'Tag',seg_str_now{cid});
        %Make frames string below checkbox
        frmtxt = frameids2str(seg_frm_now{cid});
        ftbnow = uicontrol(unq_pans(uid),'Style','edit','String',frmtxt,'FontSize',8,...
          'Units','Normalized','Position',posnow + [.03 -.025 -.04 -.02], 'Tag',seg_str_now{cid},...
          'Callback',@(s,e)editfrms_call(s,e));
      end
    end
  end
  %Load set structure
  set.recdata = recdata; set.start = start; set.stop = stop; set.refs = refs;
end

function frmsvec = editfrms_call(sorfrmstr,e)
%Get the new string
if nargin == 1
  newfrms = sorfrmstr;
else
  newfrms = sorfrmstr.String;
end
%Parse through the newstring
frmsvec = [];
newfrms = strrep(newfrms,'[','');
newfrms = strrep(newfrms,']','');
newfrms = strrep(newfrms,'- ','-');
newfrms = strrep(newfrms,' -','-');
newfrms = strrep(newfrms,' ',',');
%Split the string
spltfrms = strsplit(newfrms,',');
for s_id = 1:length(spltfrms)
  if ~isempty(spltfrms{s_id})
    if isempty(strfind(spltfrms{s_id},'-'))
      %Add single frame
      frmsvec(end+1) = str2num(spltfrms{s_id});
    elseif ~isempty(spltfrms{s_id})
      %Add series
      spltnow = strsplit(spltfrms{s_id},'-');
      frmsvec = [frmsvec str2num(spltnow{1}):str2num(spltnow{2})];
    end
  end
end
%Change edit string
if nargin == 2
  sorfrmstr.String = frameids2str(frmsvec);
end
end

function select_all(s,e)
%Findall the checkboxes
all_cbs = findall(findall(0,'Type','Figure','Tag','Setup'),'Style','Checkbox');
%Modulate based on string
sel1 = 'Select All';
sel2 = 'Deselect All';
switch s.String
  case sel1
    %Select each segment
    valset = 1;
    strset = sel2;
  case sel2
    %Deselect each
    valset = 0;
    strset = sel1;
end
for cb_id = 1:length(all_cbs)
  %Run segchk_selec
  if all_cbs(cb_id).Value ~= valset
    all_cbs(cb_id).Value = valset;
    segchk_selec([],[],all_cbs(cb_id).String)
  end
end
%Reset button string
s.String = strset;
end

function segchk_selec(s,e,seg_id,frms)
if ~exist('seg_id','var') && strcmp(s.Style,'checkbox')
  seg_id = s.String;
end
%Make setup the current figure
setfig = findall(0,'Type','Figure','Tag','Setup');
%Make text below the checkbox that notes the frms selected
cbnow = findall(setfig,'Style','checkbox','Tag',seg_id);
ftbnow = findall(setfig,'Style','edit','Tag',seg_id);
if isempty(ftbnow)
  ftbnow = findall(setfig,'Style','text','Tag',seg_id);
end

if ~exist('frms','var') || isempty(frms)
  frms = unique(editfrms_call(ftbnow.String));
else
  %If frms exist then it is a selections
  cbnow.Value = 1;
end

%Get allowed frames
%Pull from Setup figure
set = guidata(setfig);
recdata = set.recdata;
%Find recdata id to get all the frames
for r_id = 1:length(recdata)
  if strcmp(recdata(r_id).id,seg_id)
    break
  end
end
allfrms = unique(floor(recdata(r_id).Frame_IDs));
flogvec = [];%zeros(1,length(frms));
for fid = 1:length(frms)
  flogvec(fid) = any(frms(fid)==allfrms);
end
frms = frms(logical(flogvec));

%Change ftbnow to style and replace frms string
ftbnow.String = frameids2str(frms);
if isempty(frms) || cbnow.Value == 0
  ftbnow.Style = 'edit';
else
  ftbnow.Style = 'text';
end

end

%% Finishing Action of Use Map
function [] = select_plot(src,evt)
%For each selection run segchk_selec
data = guidata(gcf);
selections = data.selections;
%Close the select_frames window
delete(gcf)
%Select setup window
set = guidata(findall(0,'Type','Figure','Tag','Setup'));
%Change all edit types to []
ftbs = findall(set.fig,'Style','edit');
for ftid = 1:length(ftbs)
  ftbs(ftid).String = '';
end
for s_id= 1:length(selections.seg_id)
  %Run 
  segchk_selec([],[],selections.seg_id{s_id},selections.frames{s_id})
end
end

%% Plot Button Function
function [] = plot_start(src,evt,start,stop,refs)
setfig = findall(0,'Type','Figure','Tag','Setup');
%Get all checkboxes turned on
cbs_on = findall(setfig,'Style','checkbox','Value',1);
%Get the guidata
set = guidata(setfig);
%Get the recdata ids
rec_segs = recdata2segids(set.recdata);
%Get records fields
rec_fields = fieldnames(set.recdata);
%Parse through each turned on cb
rec_data_new  = set.recdata(1);
for cb_id = 1:length(cbs_on)
  seg_now = cbs_on(cb_id).Tag;
  %Get the frames for each
  frms_str = get(findall(setfig,'Style','text','Tag',seg_now),'String');
  %Grab each frame entry
  frmcells = strsplit(frms_str(2:end-1),',');
  frm_vec = [];
  for fc_id = 1:length(frmcells)
    frm = frmcells{fc_id};
    %Determine if range
    frm_pars = strsplit(frm,'-');
    if length(frm_pars)==2
      frm_rng = str2num(frm_pars{1}):str2num(frm_pars{2});
    else
      frm_rng = str2num(frm);
    end
    if ~isempty(frm_rng)
      frm_vec = [frm_vec frm_rng];
    end
  end
  %Get the valid rec entries based on frm_vec
  rec_logvec= strcmp(seg_now,rec_segs);
  rec_now = set.recdata(rec_logvec);
  %Parse by frm_vec
  parse_logvec = zeros(size(rec_now.Frame_IDs));
  floor_frms = floor(rec_now.Frame_IDs);
  for f = frm_vec
    parse_logvec = parse_logvec + (floor_frms==f);
  end
  parse_logvec = logical(parse_logvec);
  %Parse for each nonsingluar field
  for f_id = 1:length(rec_fields)
    if length(rec_now.(rec_fields{f_id}))==length(parse_logvec) && ~ischar(rec_now.(rec_fields{f_id}))
      %Parse based on logvec
      rec_now.(rec_fields{f_id}) = rec_now.(rec_fields{f_id})(parse_logvec);
    end
  end
  %Load to new rec_data
  rec_data_new(cb_id) = rec_now;
end
%Get start and stop based on selected recdata
[startnow, stopnow] = recdata2startstop(rec_data_new);
%Delete any view frames figures
vwfrm_figs = findall(0,'type','figure','tag','View Frames');
if ~isempty(vwfrm_figs)
  delete(vwfrm_figs);
end
%Send to view_frames
view_frames(rec_data_new,startnow,stopnow,refs);
end

%% Recdata parsing functions
function [start, stop] = recdata2startstop(recdata)
%Function makes start and stop structures based on the recdata selected
  start.lat = nan; start.lon = nan;
  stop.lat = nan; stop.lon = nan;
  for r_id = 1:length(recdata)
    minlatnow = min(recdata(r_id).Latitude); minlonnow = min(recdata(r_id).Longitude);
    maxlatnow = max(recdata(r_id).Latitude); maxlonnow = max(recdata(r_id).Longitude);
    if minlatnow < start.lat || isnan(start.lat)
      start.lat = minlatnow;
    end
    if minlonnow < start.lon || isnan(start.lon)
      start.lon = minlonnow;
    end
    if maxlatnow > stop.lat || isnan(stop.lat)
      stop.lat = maxlatnow;
    end
    if maxlonnow > stop.lon || isnan(stop.lon)
      stop.lon = maxlonnow;
    end
  end
end
function [rec_segs] = recdata2segids(recdata)
  %Function grabs the unique segment ids from recdata
  rec_cells = struct2cell(recdata);
  seg_fieldname = 'id';  rec_fields = fieldnames(recdata);
  rec_segs = rec_cells(strcmp(rec_fields,seg_fieldname),:,:);
  rec_segs = reshape(rec_segs,[1,length(rec_segs)]);
end