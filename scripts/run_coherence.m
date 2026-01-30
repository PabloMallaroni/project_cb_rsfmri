%% Generate multivariate input + CSV outputs 
close all
clear all

%% Paths and setup

datastream    = 'noGSR_schaef200';
datadir   = '/Users/administrator/Documents/MATLAB/project_cb_rest/analysis/nbs_drug_model';
nbsdir    = fullfile(datadir, ['10kperm_', datastream]);
regiondir = fullfile(datadir, ['regional_', datastream]);
savedir   = fullfile(datadir, ['multivariate_', datastream]);

if ~exist(savedir,'dir')
    mkdir(savedir)
end

measures   = {'sfc','dfc','efc_mse'};
conditions = {'cb','pla','psil'};
n_ses  = numel(conditions);
n_sub  = 20;
n_scale = 1;

raw_measures = struct;

%% Load measures


% MSE entropy
load(fullfile(regiondir,'data_regional_mse_entropy.mat'))
raw_measures.cb.efc_mse   = cb.efc;
raw_measures.psil.efc_mse = psil.efc;
raw_measures.pla.efc_mse  = pla.efc;
clear cb psil pla

% sFC / dFC
load(fullfile(nbsdir,'nbs_allFC.mat'))
raw_measures.cb.sfc   = cb.sfc;
raw_measures.psil.sfc = psil.sfc;
raw_measures.pla.sfc  = pla.sfc;
clear cb psil pla

raw_measures.cb.dfc   = cb.dfc;
raw_measures.psil.dfc = psil.dfc;
raw_measures.pla.dfc  = pla.dfc;
clear cb psil pla

%% Build multivariate regional feature matrices (z-scored)


nroi     = size(raw_measures.(conditions{1}).(measures{1}),1);
n_feat   = numel(measures);
mat_diag = eye(nroi);

similarity_dat = cell(n_sub, n_ses);

for ses = 1:n_ses
    for sub = 1:n_sub
        for type = 1:n_feat

            if strcmp(measures{type},'sfc') || strcmp(measures{type},'dfc')

                mat = raw_measures.(conditions{ses}).(measures{type})(:,:,sub);
                mat(mat_diag==1) = 0;

                similarity_dat{sub,ses}(:,type) = ...
                    zscore(sum(mat,2,'omitmissing'),0,1);

            elseif strcmp(measures{type},'efc_mse')

                similarity_dat{sub,ses}(:,type) = ...
                    zscore(raw_measures.(conditions{ses}).(measures{type})(:,sub,n_scale),0,1);

            end
        end
    end
end

save(fullfile(savedir,'similarity_dat.mat'), ...
     'similarity_dat','measures','conditions','-v7.3')

%% Feature–feature correlation plots


cb_feat   = nanmean(cat(3, similarity_dat{:,1}), 3);
pla_feat  = nanmean(cat(3, similarity_dat{:,2}), 3);
psil_feat = nanmean(cat(3, similarity_dat{:,3}), 3);
all_feat  = nanmean(cat(3, cb_feat, pla_feat, psil_feat), 3);

figure;

subplot(1,4,1)
imagesc(corr(cb_feat),[-1 1]); axis square
xticks(1:n_feat); yticks(1:n_feat)
xticklabels(measures); yticklabels(measures)
title('cb')

subplot(1,4,2)
imagesc(corr(pla_feat),[-1 1]); axis square
xticks(1:n_feat); yticks(1:n_feat)
xticklabels(measures); yticklabels(measures)
title('pla')

subplot(1,4,3)
imagesc(corr(psil_feat),[-1 1]); axis square
xticks(1:n_feat); yticks(1:n_feat)
xticklabels(measures); yticklabels(measures)
title('psil')

subplot(1,4,4)
imagesc(corr(all_feat),[-1 1]); axis square
xticks(1:n_feat); yticks(1:n_feat)
xticklabels(measures); yticklabels(measures)
title('all')

saveas(gcf, fullfile(savedir,'feature_feature_correlations.png'))


%% Region-wise multivariate similarity (raw + Fisher)


roi_mat_raw    = cell(n_sub, n_ses);
roi_mat_fisher = cell(n_sub, n_ses);

for sub = 1:n_sub
    for ses = 1:n_ses

        X = similarity_dat{sub,ses};

        R = corr(X','type','Pearson');
        R(isnan(R)) = 0;
        R(mat_diag==1) = 0;

        Z = atanh(R);
        Z(isinf(Z)) = 0;
        Z(mat_diag==1) = 0;

        roi_mat_raw{sub,ses}    = R;
        roi_mat_fisher{sub,ses} = Z;
    end
end

save(fullfile(savedir,'roi_multivariate_similarity.mat'), ...
     'roi_mat_raw','roi_mat_fisher','conditions')

%% Multivariate nodal strength + edge vectors (CSV outputs)


uppertri_mask = triu(true(nroi),1);

cb_multivariate   = zeros(n_sub,nroi);
psil_multivariate = zeros(n_sub,nroi);
pla_multivariate  = zeros(n_sub,nroi);

cb_mean   = zeros(n_sub,sum(uppertri_mask(:)));
psil_mean = zeros(n_sub,sum(uppertri_mask(:)));
pla_mean  = zeros(n_sub,sum(uppertri_mask(:)));

for sub = 1:n_sub

    tmp = roi_mat_raw{sub,1};
    tmp(isinf(tmp)) = 0; tmp(isnan(tmp)) = 0;
    cb_multivariate(sub,:) = sum(tmp,2,'omitnan')';
    cb_mean(sub,:) = tmp(uppertri_mask);

    tmp = roi_mat_raw{sub,3};
    tmp(isinf(tmp)) = 0; tmp(isnan(tmp)) = 0;
    psil_multivariate(sub,:) = sum(tmp,2,'omitnan')';
    psil_mean(sub,:) = tmp(uppertri_mask);

    tmp = roi_mat_raw{sub,2};
    tmp(isinf(tmp)) = 0; tmp(isnan(tmp)) = 0;
    pla_multivariate(sub,:) = sum(tmp,2,'omitnan')';
    pla_mean(sub,:) = tmp(uppertri_mask);

end

writematrix(cb_multivariate,   fullfile(savedir,'cb_nodal_multivariate.csv'))
writematrix(psil_multivariate, fullfile(savedir,'psil_nodal_multivariate.csv'))
writematrix(pla_multivariate,  fullfile(savedir,'pla_nodal_multivariate.csv'))

writematrix(cb_mean,   fullfile(savedir,'cb_edges_uppertri.csv'))
writematrix(psil_mean, fullfile(savedir,'psil_edges_uppertri.csv'))
writematrix(pla_mean,  fullfile(savedir,'pla_edges_uppertri.csv'))

