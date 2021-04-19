function frmstr = frameids2str(seg_frms)
%Make string of all the segment frms
frmstr = '[';
unq_frms = unique(floor(seg_frms));
frm_id = 0;
%Add first element
frm_id = frm_id + 1;
frmstr = [frmstr num2str(unq_frms(frm_id))];
ser_strt = unq_frms(frm_id); 
ser_end = ser_strt;
%Iterate through the rest
while frm_id < length(unq_frms)
  frm_id = frm_id + 1;
  frm_now = unq_frms(frm_id);
  if (frm_now-ser_end) > 1 || frm_id == length(unq_frms)
    %Add series end and series strt
    if frm_now-ser_end == 1 && ser_end ~= ser_strt
      %Add last of series
      frmstr = sprintf('%s-%.0f',frmstr,frm_now);
    elseif ser_end-ser_strt > 1
      %Add series end
      frmstr = sprintf('%s-%.0f, %.0f',frmstr,ser_end,frm_now);
      ser_strt = frm_now;
      ser_end = frm_now;
    else
      % Add last two frames
      if ser_strt~=ser_end
        frmstr = sprintf('%s, %.0f',frmstr, ser_end);
      end
      frmstr = sprintf('%s, %.0f',frmstr,frm_now);
      ser_strt = frm_now;
      ser_end = frm_now;
    end
  else
    %Set new series end
    ser_end = frm_now;
  end

end
%Add the end bracket
frmstr = [frmstr ']'];
end