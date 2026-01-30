%% Generate timeseries temporal complexity measures (Lempel–Ziv)

close all;
clear all;

%% Paths
basedir = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(basedir)
addpath(genpath(basedir));

BETdir = '/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/CopBET-master';
addpath(genpath(BETdir));

stream = 'noGSR_schaef200';

datadir = fullfile(basedir,'analysis','nbs_drug_model');
savedir = fullfile(datadir, ['regional_', stream]);

if ~exist(savedir,'dir')
    mkdir(savedir)
end

in_data = fullfile(basedir,'data','ROI_timeseries_PSC_noGSR_scrub.mat');
load(in_data)

in_parc     = 'schaefer232';
in_analysis = 'LZ76temporal';

%% Restructure data (scan-wise)
in = reshape(roi_timeseries.(in_parc), [], 1);

%% Run Lempel–Ziv temporal complexity
out = CopBET_time_series_complexity( ...
    in, ...
    in_analysis, ...
    'keepdata', false, ...
    'parallel', false);

results_lz = out.entropy;   % one value per scan

%% Save output (same data stream)
save(fullfile(savedir, ['data_lz76_temporal_entropy_', stream, '.mat']), ...
     'results_lz', )

