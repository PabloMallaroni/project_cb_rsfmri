%% get regional and whole-brain SampEn (scale = 1)

close all
clear all

%% paths
maindir = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(maindir)

addpath(genpath('/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/CopBET-master'))

stream    = 'noGSR_schaef200';
datadir   = '/Users/administrator/Documents/MATLAB/project_cb_rest/analysis/nbs_drug_model';
resultdir = fullfile(datadir, ['regional_', stream]);
savedir   = resultdir;

in_data = fullfile(datadir,'ROI_timeseries_PSC_noGSR_scrub.mat');
load(in_data)


%% configs
configs.parc = 'schaefer232';

%% restructure data
series_cell = {};
i = 1;

for cond = 1:3
    for sub = 1:20
        series_cell{i,1} = roi_timeseries.(configs.parc){sub,cond};
        i = i + 1;
    end
end

[~, num_rois] = size(series_cell{1});

%% SampEn parameters
r     = 0.3;
m     = 2;
scale = 1;

%% regional SampEn
entropy = cell(length(series_cell),1);

for ses = 1:length(series_cell)
    disp(['session ', num2str(ses)])
    for roi = 1:num_rois
        ts = series_cell{ses}(:,roi);
        r_val = r * std(double(ts));
        mse = sample_entropy(m, r_val, ts, scale);
        entropy{ses}(1,roi) = mse(1);
    end
end

entropy_out = entropy';
save(fullfile(savedir, ['data_regional_mse_entropy_', datastream, '.mat']), ...
     'entropy_out')

%% whole-brain SampEn
clear entropy

entropy = nan(length(series_cell),1);

for ses = 1:length(series_cell)
    disp(['session ', num2str(ses)])
    ts = nanmean(series_cell{ses}, 2);
    r_val = r * std(double(ts));
    mse = sample_entropy(m, r_val, ts, scale);
    entropy(ses,1) = mse(1);
end

entropy_out = entropy;
save(fullfile(savedir, ['data_wholebrain_mse_entropy_', stream, '.mat']), ...
     'entropy_out')
