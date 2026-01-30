# ============================================================
# Imports
# ============================================================

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

from scipy.stats import zscore
from scipy.spatial import distance

from dominance_analysis import Dominance
from brainspace.null_models import MoranRandomization

# ============================================================
# Paths
# ============================================================

contrast_path = "/Users/administrator/Documents/MATLAB/project_cb_rest/analysis/nbs_drug_model/"
recep_path    = "/Users/administrator/Documents/MATLAB/project_cb_rest/mni_atlas/parc_receptors/schaefer232/"
parc_path     = "/Users/administrator/Documents/MATLAB/project_cb_rest/mni_atlas/parc_networks/"
plot_path     = "/Users/administrator/Documents/MATLAB/project_cb_rest/analysis/nbs_drug_model/figures/"

# ============================================================
# Load nodal outcome maps (NO gFC, only sFC / dFC / SampEn)
# ============================================================

# sFC
sfc_cb   = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_sfc_zerothresh_cb.csv',   header=None, names=['sfc_cb'])
sfc_psil = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_sfc_zerothresh_psil.csv', header=None, names=['sfc_psil'])
sfc_diff = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_sfc_zerothresh_diff.csv', header=None, names=['sfc_diff'])

# dFC
dfc_cb   = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_dfc_zerothresh_cb.csv',   header=None, names=['dfc_cb'])
dfc_psil = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_dfc_zerothresh_psil.csv', header=None, names=['dfc_psil'])
dfc_diff = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nodal_tstat_dfc_zerothresh_diff.csv', header=None, names=['dfc_diff'])

# SampEn / eFC
efc_cb   = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nothresh_tstat_roi_cb_efc.csv',   header=None, names=['sampEn_cb'])
efc_psil = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nothresh_tstat_roi_psil_efc.csv', header=None, names=['sampEn_psil'])
efc_diff = pd.read_csv(contrast_path+'regional_noGSR_schaef200/nothresh_tstat_roi_diff_efc.csv', header=None, names=['sampEn_diff'])

# ============================================================
# Combine outcomes
# ============================================================

outcomes_df = pd.concat(
    [
        sfc_psil, sfc_cb, sfc_diff,
        dfc_psil, dfc_cb, dfc_diff,
        efc_psil, efc_cb, efc_diff
    ],
    axis=1
)

outcomes_df = outcomes_df.apply(zscore)

# ============================================================
# Load receptor data
# ============================================================

with open(recep_path+'receptor_names.csv') as f:
    r_names = [x.strip() for x in f.readlines()]

r_dat = pd.read_csv(recep_path+'receptor_data.csv', header=None, names=r_names)

# z-score cortex and subcortex separately
r_dat.iloc[:200] = zscore(r_dat.iloc[:200], axis=0)
r_dat.iloc[200:] = zscore(r_dat.iloc[200:], axis=0)

# drop unused receptors
r_dat = r_dat.drop(columns=['5HT6','5HT4'])

# ============================================================
# Moran null model
# ============================================================

coords = pd.read_table(os.path.join(parc_path,'schaefer232_coords.txt'), header=None)

dist_mat = distance.squareform(distance.pdist(coords, 'euclidean'))
W = np.exp(-dist_mat)
np.fill_diagonal(W, 0)

msr = MoranRandomization(
    n_rep=10000,
    procedure='singleton',
    tol=1e-6,
    random_state=0
)
msr.fit(W)

# ============================================================
# Dominance analysis per contrast
# ============================================================

predictors = r_dat
dominance_results = {}
top_pairs_with_stats = {}
data_frames = {}

for contrast in ['psil', 'cb', 'diff']:

    cols = [c for c in outcomes_df.columns if c.endswith(contrast)]
    dominance_results[contrast] = {}
    top_pairs_with_stats[contrast] = []

    for col in cols:
        dom = Dominance(
            data=pd.concat([predictors, outcomes_df[col]], axis=1),
            target=col,
            objective=1
        )
        dom.incremental_rsquare()
        dom.dominance_stats()
        dominance_results[contrast][col] = dom

    # Build radar DataFrame
    radar_df = pd.DataFrame({
        k: v.dominance_stats()['Percentage Relative Importance']
        for k,v in dominance_results[contrast].items()
    })
    data_frames[contrast] = radar_df

# ============================================================
# Radar + scatter plotting
# ============================================================

labelsize = 18
ticksize  = 16
textsize  = 14

def plot_radar(ax, df, title):
    cats = df.index.tolist()
    N = len(cats)
    angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(cats, fontsize=labelsize)
    ax.set_ylim(0,100)
    ax.set_yticks([25,50,75])
    ax.set_yticklabels(['25','50','75'], fontsize=labelsize)

    colors = sns.color_palette('tab10', df.shape[1])
    cmap = {}

    for i,col in enumerate(df.columns):
        vals = df[col].tolist() + [df[col].iloc[0]]
        ax.plot(angles, vals, lw=2, color=colors[i], label=col)
        ax.fill(angles, vals, alpha=0.15, color=colors[i])
        cmap[col] = colors[i]

    ax.set_title(title, fontsize=labelsize)
    ax.legend(frameon=False, fontsize=labelsize, ncol=df.shape[1],
              loc='upper center', bbox_to_anchor=(0.5,-0.15))
    return cmap

fig = plt.figure(figsize=(24, 5*len(data_frames)))
outer = gridspec.GridSpec(len(data_frames), 4, wspace=0.25)

for r,(contrast,radar_df) in enumerate(data_frames.items()):

    ax_radar = fig.add_subplot(outer[r,0], polar=True)
    color_map = plot_radar(ax_radar, radar_df, contrast)

    for j,(metric) in enumerate(radar_df.columns):
        ax = fig.add_subplot(outer[r,j+1])

        y = outcomes_df[metric]
        best_receptor = radar_df[metric].idxmax()
        x = predictors[best_receptor]

        r_obs = np.corrcoef(x, y)[0,1]
        r_null = msr.randomize(x.values)
        p_moran = (np.abs(np.corrcoef(r_null.T, y)[0,1:]) >= abs(r_obs)).mean()

        sns.regplot(
            x=x, y=y, ax=ax,
            scatter_kws={'s':20,'alpha':0.5,'color':color_map[metric]},
            line_kws={'color':color_map[metric]}
        )

        ax.set_xlabel('Normalised receptor density', fontsize=labelsize)
        ax.set_ylabel('Normalised nodal t-stat', fontsize=labelsize)
        ax.tick_params(labelsize=ticksize)
        ax.set_title(f'{best_receptor}', fontsize=labelsize)

        ax.text(
            0.05,0.95,
            f'r = {r_obs:.2f}\nMoran p = {p_moran:.3f}',
            transform=ax.transAxes,
            va='top', fontsize=textsize
        )

plt.savefig(
    os.path.join(plot_path,'dominance_spider_scatter_sfc_dfc_sampen.png'),
    dpi=600, bbox_inches='tight'
)
plt.show()
