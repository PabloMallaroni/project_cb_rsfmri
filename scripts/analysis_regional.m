%% Regional analysis: main effect + contrasts
close all
clear all

maindir   = '/Users/administrator/Documents/MATLAB/project_cb_rest';
cd(maindir)

resultdir = fullfile(maindir,'analysis','nbs_drug_model','regional_noGSR_schaef200');
addpath(genpath(maindir))

configs.parc = 'schaefer232';
conditions   = {'cb','pla','psil'};
n_cond       = numel(conditions);
alpha_desc   = 0.05;

netfile = xlsread( ...
    '/Users/administrator/Documents/MATLAB/project_cb_rest/mni_atlas/parc_networks/schaefer18networks_idx.xlsx', ...
    configs.parc);

yeoROI = netfile(:,1);
n_roi  = max(yeoROI);

measures = {'efc_mse','gfc'};

for m = 1:numel(measures)

    measure = measures{m};
    disp(['Running regional analysis for ',measure])

    switch measure

        case 'efc_mse'
            load(fullfile(resultdir,'values_regional_mse_entropy_noGSR.mat'))
            n_ses = 3;
            n_sub = length(entropy_out)/n_ses;

            data = cell(n_sub,n_ses);
            data(:,1) = entropy_out(1:n_sub);
            data(:,2) = entropy_out(n_sub+1:2*n_sub);
            data(:,3) = entropy_out(2*n_sub+1:3*n_sub);

        case 'gfc'
            load(fullfile(resultdir,'data_regional_gFC_nbs.mat'))
            data = cell(n_sub,3);
            for s = 1:n_sub
                data{s,1} = cb.gfc(:,s)';
                data{s,2} = pla.gfc(:,s)';
                data{s,3} = psil.gfc(:,s)';
            end
            n_ses = 3;
    end

    idx_sub  = repmat((1:n_sub)',n_ses,1);
    idx_cond = [ones(n_sub,1); zeros(n_sub,1); ones(n_sub,1)];

    p_main = zeros(n_roi,1);
    f_main = zeros(n_roi,1);

    for roi = 1:n_roi
        Y = [];
        for ses = 1:n_ses
            tmp = cellfun(@(x) x(roi), data(:,ses));
            Y   = [Y; tmp];
        end

        tbl = table(Y, idx_sub, idx_cond, ...
            'VariableNames', {'y','subject','condition'});

        lme = fitlme(tbl,'y ~ condition + (1|subject)', ...
                     'DummyVarCoding','effects');

        a = anova(lme);
        p_main(roi) = a.pValue(2);
        f_main(roi) = a.FStat(2);
    end

    p_main_fdr = mafdr(p_main,'BHFDR',true);
    sig_roi    = find(p_main_fdr < 0.05);

    t_cb   = zeros(n_roi,1);
    p_cb   = ones(n_roi,1);

    t_psil = zeros(n_roi,1);
    p_psil = ones(n_roi,1);

    t_diff = zeros(n_roi,1);
    p_diff = ones(n_roi,1);

    for r = sig_roi'

        cb_vals   = cellfun(@(x) x(r), data(:,1));
        pla_vals  = cellfun(@(x) x(r), data(:,2));
        psil_vals = cellfun(@(x) x(r), data(:,3));

        [~,p,~,st] = ttest(cb_vals,pla_vals);
        t_cb(r) = st.tstat; p_cb(r) = p;

        [~,p,~,st] = ttest(psil_vals,pla_vals);
        t_psil(r) = st.tstat; p_psil(r) = p;

        [~,p,~,st] = ttest(psil_vals,cb_vals);
        t_diff(r) = st.tstat; p_diff(r) = p;
    end

    t_cb(p_cb >= alpha_desc)     = 0;
    t_psil(p_psil >= alpha_desc) = 0;
    t_diff(p_diff >= alpha_desc) = 0;

    save(fullfile(resultdir, ...
        ['stats_regional_',measure,'_noGSR.mat']), ...
        't_cb','p_cb', ...
        't_psil','p_psil', ...
        't_diff','p_diff', ...
        'f_main','p_main','p_main_fdr','sig_roi')

end
