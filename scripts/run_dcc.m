%% Generate DCC-based measures (long-running, use a cluster - seriously!!!)

close all;
clear all;

%% Paths
basedir = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(basedir)
addpath(genpath(basedir));

BETdir = '/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/CopBET-master';
addpath(genpath(BETdir));

stream = 'PSC_noGSR_sch232';

datadir = fullfile(basedir,'analysis','nbs_drug_model');
savedir = fullfile(datadir, ['10kperm_', stream]);

if ~exist(savedir,'dir')
    mkdir(savedir)
end

in_data = fullfile(basedir,'data','ROI_timeseries_PSC_noGSR_scrub.mat');
load(in_data)

%% Parameters
in_parc = 'schaefer232';
iter    = 5000;

%% Restructure data (scan-wise)
series_cell = reshape(roi_timeseries.(in_parc), [], 1);
n_ses = numel(series_cell);

%% Start parallel pool
num_workers = 28;
p = gcp('nocreate');
if isempty(p)
    parpool(num_workers);
end

%% Run DCC per scan (writes one file per scan)

parfor sesidx = 1:n_ses
    disp(['DCC scan ', num2str(sesidx)])
    Barrett_analysis_no_correction_sps253( ...
        series_cell{sesidx}, ...
        sesidx, ...
        savedir, ...
        stream, ...
        iter);
end


%% Aggregate DCC outputs

disp('Aggregating DCC outputs...')

n_sub  = size(roi_timeseries.(in_parc),1);
n_cond = size(roi_timeseries.(in_parc),2);

entropy_mat  = cell(n_sub, n_cond);
variance_mat = cell(n_sub, n_cond);

sub  = 1;
cond = 1;

for ses = 1:n_ses

    load(fullfile(savedir, ...
        ['DCC_',stream,'_ses',num2str(ses),'.mat']), ...
        'entropy_out','variance_out');

    entropy_mat{sub,cond}  = entropy_out;
    variance_mat{sub,cond} = variance_out;

    if sub == n_sub
        sub = 1;
        cond = cond + 1;
    else
        sub = sub + 1;
    end
end

%% Save clean, final outputs

save(fullfile(savedir, ...
    ['DCC_entropy_variance_',stream,'.mat']), ...
    'entropy_mat','variance_mat','in_parc','stream','-v7.3');

disp('DCC aggregation complete and saved.')
