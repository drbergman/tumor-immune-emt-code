function [tracked_names,per_Nums] = track_ind2names(tracker_index)

names = {'Num',...
         'NK',...
         'CTL',...
         'Treg',...
         'Mut_types',...
         'Mut',...
         'TGFB',...
         'EMT_val_quarts',...
         'EMT_bincounts',...
         'Mes'};
     
tracked_names = names(tracker_index);

if ~isempty(find(tracker_index==1,1)) % tracking number of cells, so can normalize by this
    not_normed = {'EMT_val_quarts'}; % can be appended if I decide other tracked values should not be normalized
    for i = 1:length(tracker_index)
        ti = tracker_index(i);
        s = strcmp(names{ti},not_normed);
        if ~any(s)
            per_Nums.(names{ti}) = true;
        else
            per_Nums.(names{ti}) = false;
        end
    end
else
    per_Nums = [];
end
            