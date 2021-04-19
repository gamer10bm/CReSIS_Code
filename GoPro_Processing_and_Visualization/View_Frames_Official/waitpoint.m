function [strt_point] = waitpoint(fig_handle, oldpoint)
strt_point = get(fig_handle,'pointer');
%Modulate based on oldpoint
if ~exist('oldpoint','var') || isempty(oldpoint) 
  if ~strcmp(strt_point,'watch')
    %Turn to watch
    pointstr = 'watch';
  else
    pointstr = 'arrow';
  end
else
  %Turn to given pointer
  pointstr = oldpoint;
end
%Set pointer
set(fig_handle,'pointer',pointstr)
%Draw now
drawnow
end
  