%% Generate static FC (edge-level) and global FC (regional-level)

close all
clear all

%% Paths
basedir = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(basedir)
addpath(genpath(basedir))

datadir = fullfile(basedir,'data');

%% Datastream
stream = 'PSC_noGSR_sch232';

savedir_edge   = fullfile(basedir,'analysis','nbs_drug_model',['10kperm_',stream]);
savedir_region = fullfile(basedir,'analysis','nbs_drug_model',['regional_',stream]);

if ~exist(savedir_edge,'dir'),   mkdir(savedir_edge);   end
if ~exist(savedir_region,'dir'), mkdir(savedir_region); end

%% Load data
load(fullfile(datadir,'ROI_timeseries_PSC_noGSR_scrub.mat'))

%% Config
parc = 'schaefer232';
in_series = roi_timeseries.(parc);

[nsub,ncond] = size(in_series);
nroi = size(in_series{1,1},2);

mat_diag = eye(nroi);

%% Initialise outputs
roi_mat_static = cell(nsub,ncond);   % ROI × ROI (edge-level)
roi_gfc        = cell(nsub,ncond);   % ROI × 1  (regional-level)

%% Compute static FC and global FC

for sub = 1:nsub
    for ses = 1:ncond

        if ~isempty(in_series{sub,ses})

            ts = double(in_series{sub,ses});   % TR × ROI

            % Pearson correlation
            fc = corr(ts);

            % Handle NaNs
            fc(isnan(fc)) = 0;

            % Fisher z-transform
            fc = atanh(fc);

            % Remove diagonal
            fc(mat_diag==1) = NaN;

            % Save edge-level sFC
            roi_mat_static{sub,ses} = fc;

            % Save regional gFC (mean nodal strength)
            roi_gfc{sub,ses} = mean(fc,2,'omitnan');

        else
            roi_mat_static{sub,ses} = NaN(nroi,nroi);
            roi_gfc{sub,ses}        = NaN(nroi,1);
        end
    end
end


save(fullfile(savedir_edge, ...
    ['sFC_edge_',stream,'.mat']), ...
    'roi_mat_static');

save(fullfile(savedir_region, ...
    ['gFC_regional_',stream,'.mat']), ...
    'roi_gfc');

