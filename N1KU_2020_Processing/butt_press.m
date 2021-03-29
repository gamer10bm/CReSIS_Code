function [] = butt_press(handle, event, plotfcn)
%Load data
dat = guidata(handle);

switch event.Key
  case 'rightarrow'
    dat.fn_id = dat.fn_id +1;
  case 'leftarrow'
    dat.fn_id = dat.fn_id -1;
end

%update plot
dat.fn_id = plotfcn(dat.fn_id);
guidata(handle,dat)
end
