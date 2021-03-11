function EEG = extended_source_atlas(EEG,sourcemodel,musclegrid)
%%
try
    EEG = eeg_compatlas(EEG);
catch
    warning('compatlas not calculated due to error')
end
% warning('Adjust path in "extended_source_atlas"!')
% % load('C:\repos\muscle-head\durchbruch_scripts\extended_BEM\musclemodel\muscle2580grid_in_mm.mat');
% % load('C:\repos\muscle-head\durchbruch_scripts\extended_BEM\musclemodel\muscle10648grid_in_mm.mat');
% load('C:\repos\muscle-head\durchbruch_scripts\extended_BEM\musclemodel\muscle3540grid_in_mm.mat');
%
% % load('D:\muscle-head\durchbruch_scripts\extended_BEM\sourcemodel_standard.mat');
% load('C:\repos\muscle-head\durchbruch_scripts\extended_BEM\sourcemodel_brain_only.mat');
% % load('C:\repos\muscle-head\durchbruch_scripts\extended_BEM\sourcemodel_brain_only_fine.mat');

%% coarse fit
% sourcemodel.pos = [sourcemodel.pos; musclegrid.pos];
inside = is_inside(cell2mat({EEG.dipfit.model.posxyzcoarse}'),  EEG.etc.bnd(end));

% EEG.dipfit.model(inside).tissuecoarse = ['Brain_' EEG.dipfit.model(inside).areadk];
for i_dip = 1:length(EEG.dipfit.model)
    if inside(i_dip)
        EEG.dipfit.model(i_dip).tissuecoarse = ['Brain_' EEG.dipfit.model(i_dip).areadk];
        distances = sourcemodel.pos - EEG.dipfit.model(i_dip).posxyzcoarse;
        [EEG.dipfit.model(i_dip).closestdipoledistancecoarse, ~] = min(sqrt(sum(distances.^2,2)));
    else
        distances = musclegrid.pos - EEG.dipfit.model(i_dip).posxyzcoarse;
        
        [EEG.dipfit.model(i_dip).closestdipoledistancecoarse, idx] = min(sqrt(sum(distances.^2,2)));
        % is muscle
        EEG.dipfit.model(i_dip).tissuecoarse = strtrim(musclegrid.tissue(idx,:));
        
    end
    
end

%% fine fit
inside = bounding_mesh(cell2mat({EEG.dipfit.model.posxyz}'),  EEG.etc.bnd(end) .pos,  EEG.etc.bnd(end).tri)==1;

% EEG.dipfit.model(inside).tissuecoarse = strcat('Brain_', EEG.dipfit.model(inside).areadk);
for i_dip = 1:length(EEG.dipfit.model)
    if inside(i_dip)
        EEG.dipfit.model(i_dip).tissue = ['Brain_' EEG.dipfit.model(i_dip).areadk];
        distances = sourcemodel.pos - EEG.dipfit.model(i_dip).posxyz;
        [EEG.dipfit.model(i_dip).closestdipoledistance, ~] = min(sqrt(sum(distances.^2,2)));
    else
        distances = musclegrid.pos - EEG.dipfit.model(i_dip).posxyz;
        
        [EEG.dipfit.model(i_dip).closestdipoledistance, idx] = min(sqrt(sum(distances.^2,2)));
        % is muscle
        EEG.dipfit.model(i_dip).tissue = strtrim(musclegrid.tissue(idx,:));
        
    end
    
end


