%% Setup paths

close all
clear all

addpath(genpath('/Users/administrator/Documents/MATLAB/imaging_tools/toolboxes/rotate_parcellation-master'));

parcdir = '/Users/administrator/Documents/MATLAB/project_cb_rest/mni_atlas/parc_networks';

parc_cort = fullfile(parcdir, ...
    'space-MNI152NLin2009cAsym_atlas-Schaefer217_res-2_dseg.nii');

%example approach
plot_cort = readmatrix(fullfile(parcdir, 'sydnor2021_SAaxis_sch200_17net.txt'));

%% Load parcellation and map values to volume

GM_parc = load_nii(parc_cort);
nReg = max(GM_parc.img(:));

RC_surface = zeros(size(GM_parc.img));

for roi = 1:nReg
    RC_surface(GM_parc.img == roi) = plot_cort(roi);
end

GM_parc.img = RC_surface;

niftiPath  = 'cort_gradient.nii';
save_nii(GM_parc, niftiPath)

%% Visualize with BrainNet

[brainnet_path, ~, ~] = fileparts(which('BrainNet.m'));
in_surface = fullfile(brainnet_path, ...
    'Data','SurfTemplate','BrainMesh_ICBM152_smoothed.nv');

BrainNet_MapCfg(in_surface, niftiPath);

%% Spin permutation testing

spin_coord = readmatrix(fullfile(parcdir, 'schaefer232_coords.txt'));
spin_coord = spin_coord(1:200,:);   % cortex only

perm_id = rotate_parcellation( ...
    spin_coord(1:100,:), ...
    spin_coord(101:200,:), ...
    10000);

%% Correlation + spin test

pls_psych = readmatrix('loadings_z_overall.csv');
spin_psych = pls_psych(1:200);

spin_gradient = zscore(plot_cort);

[r,p] = corr(spin_psych, spin_gradient, 'type','Pearson');

pspin = perm_sphere_p( ...
    spin_psych, ...
    spin_gradient, ...
    perm_id, ...
    'Pearson');

disp(['r = ', num2str(r), ', p = ', num2str(p), ', spin-p = ', num2str(pspin)])
