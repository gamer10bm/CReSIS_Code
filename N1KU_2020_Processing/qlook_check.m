function [fn_id] = qlook_check(param, qlook_suffs, fn_id)
%grab figure
f = gcf;
%% Inspect and load basic qlook
    qlook={};
    qlook_dir = ct_filename_out(param,'qlook');
    qlook_fns = dir(qlook_dir);
    
    %Get the qlook_id
    while any(strcmp(qlook_fns(fn_id).name,{'.','..'}))
      fn_id = fn_id +1;
    end
    fn = qlook_fns(fn_id).name;
    
    qlook{end+1} = load(fullfile(qlook_dir,fn));
    
    %% Load other qlooks
    for id_s = 1:length(qlook_suffs)
      qlook_dir = ct_filename_out(param,sprintf('qlook_%s',qlook_suffs{id_s}));
      qlook{end+1} = load(fullfile(qlook_dir,fn));
    end
    
    %% Plot the things
    set(f,'name',fn)
    ql = length(qlook);
    ax = [];
    for qid = 1:ql
      ax(qid) = subplot(ql,1,qid);
      imagesc(lp(qlook{qid}(1).Data))
      if qid == 1
        title('Qlook')
      else
        title(qlook_suffs{qid-1})
      end
    end
    linkaxes(ax);
end