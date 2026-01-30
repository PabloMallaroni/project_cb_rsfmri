%% NBS DRUG FULL
close all
clear all

%% Paths
maindir = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(maindir)

datadir   = fullfile(maindir,'data');
resultdir = fullfile(maindir,'analysis','nbs_drug_model');

addpath(genpath(maindir))
addpath('/Users/administrator/Documents/MATLAB/project_cb_rest/toolbox/NBS1.2');
addpath('/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/BrainNetViewer');
addpath('/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/NIfTI+toolbox');

%% Config
configs.parc = 'schaefer232';

%% NBS parameters
NBS_param.size    = 'Intensity';
NBS_param.thresh  = 4.1; %Based on apriori effect size
NBS_param.test    = 'F-test';
NBS_param.permute = 10000;
alpha = 0.05;

NBS_type = '10kperm_noGSR_schaef200';

%% Load data
load(fullfile(datadir,'ROI_timeseries_PSC_noGSR_scrub.mat'))

sfc_file = fullfile( ...
    maindir, ...
    'analysis','nbs_drug_model', ...
    ['10kperm_',stream], ...
    ['sFC_edge_',stream,'.mat'] );

load(sfc_file, 'roi_mat_static');



dcc_file = fullfile(result_dir, ['10kperm_',stream],...
    ['DCC_entropy_variance_', stream, '.mat']);

load(dcc_file, 'entropy_mat', 'variance_mat')

roi_mat_entropic = entropy_mat;
roi_mat_dynamic  = variance_mat;



n_sub = size(roi_mat_static,1);
n_ses = size(roi_mat_static,2);
n_roi = size(roi_mat_static{1,1},1);

roi_upper = triu(true(n_roi),1);

%% Prepare NBS input (subjects × conditions → edges)
dat_nbs_sfc = cat(3, roi_mat_static{:});
dat_nbs_dfc = cat(3, roi_mat_dynamic{:});
dat_nbs_efc = cat(3, roi_mat_entropic{:});

%% Load NBS results
nbs_path = fullfile(resultdir,NBS_type);

nbs_sfc = load(fullfile(nbs_path,'result_nbs_sfc_fd.mat'));
nbs_dfc = load(fullfile(nbs_path,'result_nbs_dfc_fd.mat'));
nbs_efc = load(fullfile(nbs_path,'result_nbs_efc.mat'));

nbs_array = {nbs_sfc, nbs_dfc, nbs_efc};
nbs_names = {'drug × sFC','drug × dFC','drug × eFC'};

%% Visualise NBS-restricted ROI matrices (no network averaging)

figure
for m = 1:numel(nbs_array)

    stat = nbs_array{m}.nbs.NBS.test_stat;
    pval = nbs_array{m}.nbs.NBS.pval;

    stat(abs(stat) < NBS_param.thresh) = 0;

    subplot(1,3,m)
    imagesc(stat)
    axis square
    colorbar
    title([nbs_names{m}, '  p = ', num2str(pval)])

end

%% Post-hoc: ROI-level edge masking based on NBS component

%     All subsequent analyses are strictly restricted to edges belonging
%     to the NBS-identified component. Edge-wise t-statistics and
%     uncorrected p-values (p < 0.05) are used *only* to characterize the
%     internal structure, sign, and spatial distribution of effects
%     within the NBS component.
%
%     These p-values are NOT interpreted as independent hypothesis tests
%     and are NOT subjected to additional multiple-comparisons correction,
%     because:
%       - the search space has already been constrained by NBS,
%       - edges are not statistically independent,
%       - NBS already provides family-wise error control at the component
%         level.


clear idx
idx = nbs_sfc.nbs.NBS.test_stat;
idx(abs(idx) < NBS_param.thresh) = 0;
idx(abs(idx) >= NBS_param.thresh) = 1;

for sub = 1:n_sub
    for ses = 1:n_ses
        tmp = roi_mat_static{sub,ses};
        tmp(~idx) = 0;
        thresh_roi_mat_static{sub,ses} = tmp;
    end
end

clear idx
idx = nbs_dfc.nbs.NBS.test_stat;
idx(abs(idx) < NBS_param.thresh) = 0;
idx(abs(idx) >= NBS_param.thresh) = 1;

for sub = 1:n_sub
    for ses = 1:n_ses
        tmp = roi_mat_dynamic{sub,ses};
        tmp(~idx) = 0;
        thresh_roi_mat_dynamic{sub,ses} = tmp;
    end
end

%% Pairwise ROI-level t-tests

for sub = 1:n_sub
    cb.sfc(:,:,sub)   = thresh_roi_mat_static{sub,1};
    pla.sfc(:,:,sub)  = thresh_roi_mat_static{sub,2};
    psil.sfc(:,:,sub) = thresh_roi_mat_static{sub,3};

    cb.dfc(:,:,sub)   = thresh_roi_mat_dynamic{sub,1};
    pla.dfc(:,:,sub)  = thresh_roi_mat_dynamic{sub,2};
    psil.dfc(:,:,sub) = thresh_roi_mat_dynamic{sub,3};
end

%% ROI-level contrasts 

[~,p,~,stats] = ttest(psil.sfc, pla.sfc, 'dim',3,'tail','both');
t_sfc_psil = stats.tstat; p_sfc_psil = p;

[~,p,~,stats] = ttest(cb.sfc, pla.sfc, 'dim',3,'tail','both');
t_sfc_cb = stats.tstat; p_sfc_cb = p;

[~,p,~,stats] = ttest(psil.dfc, pla.dfc, 'dim',3,'tail','both');
t_dfc_psil = stats.tstat; p_dfc_psil = p;

[~,p,~,stats] = ttest(cb.dfc, pla.dfc, 'dim',3,'tail','both');
t_dfc_cb = stats.tstat; p_dfc_cb = p;

%% Mask by p < 0.05 (descriptive, within NBS)

t_sfc_psil(p_sfc_psil >= alpha) = 0;
t_sfc_cb(p_sfc_cb >= alpha)     = 0;
t_dfc_psil(p_dfc_psil >= alpha) = 0;
t_dfc_cb(p_dfc_cb >= alpha)     = 0;

%% Save ROI-level outputs

writematrix(t_sfc_psil, fullfile(nbs_path,'roi_tmat_sfc_psil_nbs_p05.csv'))
writematrix(t_sfc_cb,   fullfile(nbs_path,'roi_tmat_sfc_cb_nbs_p05.csv'))
writematrix(t_dfc_psil, fullfile(nbs_path,'roi_tmat_dfc_psil_nbs_p05.csv'))
writematrix(t_dfc_cb,   fullfile(nbs_path,'roi_tmat_dfc_cb_nbs_p05.csv'))

