%% Generate degree distribution entropy

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

in_parc = 'schaefer232';

%% Restructure data (scan-wise)
in = reshape(roi_timeseries.(in_parc), [], 1);

%% Run degree distribution entropy
out = CopBET_degree_distribution_entropy( ...
    in, ...
    'keepdata', false, ...
    'parallel', false);

results_degree = out.entropy;   % one value per scan

%% Save output (s
save(fullfile(savedir, ['data_degree_entropy_', stream, '.mat']), ...
     'results_degree')
