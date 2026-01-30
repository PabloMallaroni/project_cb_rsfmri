# ============================================================
# Standalone Behavioral PLS: Regional Drug–Placebo Contrasts
# ============================================================

import os
import random
import numpy as np
import pandas as pd
import pyls
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt

from scipy.stats import zscore
from pyls import utils, compute
from matplotlib import cm

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

base_dir = '/Users/administrator/Documents/MATLAB/project_cb_rest'
parc_dir = os.path.join(base_dir, 'mni_atlas', 'parc_networks')
out_dir  = os.path.join(base_dir, 'analysis', 'nbs_drug_model')

labels_path = os.path.join(parc_dir, 'schaefer232_labels.txt')

# ------------------------------------------------------------
# Load ROI labels
# ------------------------------------------------------------

with open(labels_path, 'r') as f:
    labels = f.read().splitlines()

n_roi = len(labels)

# ------------------------------------------------------------
# Load regional multivariate data (X)
# ------------------------------------------------------------

pla_multivariate  = pd.read_csv(os.path.join(out_dir,'pla_multivariate.csv'), header=None)
cb_multivariate   = pd.read_csv(os.path.join(out_dir,'cb_multivariate.csv'), header=None)
psil_multivariate = pd.read_csv(os.path.join(out_dir,'psil_multivariate.csv'), header=None)

pla_multivariate.columns  = labels
cb_multivariate.columns   = labels
psil_multivariate.columns = labels

# ------------------------------------------------------------
# Load behavioral effects (Y)
# ------------------------------------------------------------

pla_effects  = pd.read_csv(os.path.join(out_dir,'pla_effects.csv'))
cb_effects   = pd.read_csv(os.path.join(out_dir,'cb_effects.csv'))
psil_effects = pd.read_csv(os.path.join(out_dir,'psil_effects.csv'))

# ------------------------------------------------------------
# Construct drug–placebo contrasts
# ------------------------------------------------------------

cbpla_multivariate   = cb_multivariate   - pla_multivariate
psilpla_multivariate = psil_multivariate - pla_multivariate

cbpla_effects   = cb_effects   - pla_effects
psilpla_effects = psil_effects - pla_effects

# Stack contrasts (CB–PLA followed by PSIL–PLA)
X = pd.concat([cbpla_multivariate, psilpla_multivariate], axis=0)
Y = pd.concat([cbpla_effects, psilpla_effects], axis=0)

# ------------------------------------------------------------
# PLS configuration
# ------------------------------------------------------------

n_sub  = cbpla_multivariate.shape[0]
groups = [n_sub]          # one group, two conditions
n_cond = 2

random.seed(8)

# ------------------------------------------------------------
# Run behavioral PLS
# ------------------------------------------------------------

out = pyls.behavioral_pls(
    X, Y,
    groups=groups,
    n_cond=n_cond,
    seed=8,
    n_perm=10000,
    n_boot=10000
)

print('Permutation p-values:', out.permres.pvals)
print('Variance explained:', out.varexp)

# ============================================================
# BOOTSTRAP LOADINGS (regional stability)
# ============================================================

# CB loadings
cb_boot = out.bootres.x_loadings_boot[:, :n_roi, 0].T
cb_bsr  = np.mean(cb_boot, axis=1) / np.std(cb_boot, axis=1)
cb_bsr  = pd.DataFrame(cb_bsr, index=labels, columns=['BSR'])

# PSIL loadings
psil_boot = out.bootres.x_loadings_boot[:, n_roi:, 0].T
psil_bsr  = np.mean(psil_boot, axis=1) / np.std(psil_boot, axis=1)
psil_bsr  = pd.DataFrame(psil_bsr, index=labels, columns=['BSR'])

# Overall loadings
overall_bsr = pd.DataFrame(
    out.bootres.x_weights_normed[:, 0],
    index=labels,
    columns=['BSR']
)

# ------------------------------------------------------------
# Display strongest contributors
# ------------------------------------------------------------

print('\nTop CB regions:')
print(cb_bsr.nlargest(10, 'BSR'))

print('\nTop PSIL regions:')
print(psil_bsr.nlargest(10, 'BSR'))

print('\nTop overall regions:')
print(overall_bsr.nlargest(10, 'BSR'))

# ============================================================
# Plot: Top 20 overall regions
# ============================================================

top_regions = overall_bsr.reindex(
    overall_bsr['BSR'].abs().sort_values(ascending=False).index[:20]
)

colors = [
    cm.coolwarm(0.1) if v < 0 else cm.coolwarm(0.9)
    for v in top_regions['BSR']
]

plt.figure(figsize=(12, 8))
top_regions['BSR'].plot(kind='barh', color=colors)
plt.xlabel('Bootstrap ratio')
plt.title('Top 20 regions (overall PLS)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.show()

# ============================================================
# Variance explained vs permutation null
# ============================================================

null_sing = out.permres.permsingular
real_exp  = out.varexp * 100
real_p    = out.permres.pvals

null_exp = []
for i in range(null_sing.shape[1]):
    d = np.diag(null_sing[:, i])
    null_exp.append(np.diag(compute.varexp(d)))

null_exp = np.array(null_exp) * 100

plt.figure(figsize=(7, 7))
sns.boxplot(
    data=null_exp[:, :6],
    color='gray',
    showfliers=False,
    width=0.4,
    orient='h'
)

for i in range(min(len(real_exp), 6)):
    plt.scatter(
        real_exp[i], i,
        color='red' if real_p[i] < 0.05 else 'gray',
        s=80,
        zorder=5
    )

plt.yticks(range(6), [f'LV{i+1}' for i in range(6)])
plt.xlabel('Variance explained (%)')
plt.ylabel('Latent variable')
plt.tight_layout()
plt.show()

# ============================================================
# Network-level aggregation of raw loadings (descriptive)
# ============================================================

net_idx = pd.read_excel(
    os.path.join(parc_dir, 'schaefer18networks_idx.xlsx'),
    sheet_name='schaefer232'
)

net_labels = pd.read_excel(
    os.path.join(parc_dir, 'schaefer18networks_idx.xlsx'),
    sheet_name='net_label'
)

groups_dummy = utils.dummy_code(out.inputs.groups, out.inputs.n_cond)
Xcov = np.row_stack([
    compute.xcorr(X[g], out['x_scores'][g], covariance=out.inputs.covariance)
    for g in groups_dummy.T.astype(bool)
])

loadings = pd.DataFrame(Xcov[:n_roi, 0], columns=['loading'])
loadings['roi'] = np.arange(n_roi)

merged = net_idx.merge(loadings, on='roi')
grouped = merged.groupby('net')['loading']

mean_load = grouped.mean()
sem_load  = grouped.apply(stats.sem)

plt.figure(figsize=(10, 8))
plt.errorbar(
    range(len(mean_load)),
    mean_load,
    yerr=1.5 * sem_load,
    fmt='o',
    color='black'
)

plt.xticks(
    range(len(mean_load)),
    net_labels.set_index('net').loc[mean_load.index, 'label'],
    rotation=45,
    ha='right'
)

plt.ylabel('Mean loading')
plt.title('Network-level PLS loadings')
plt.tight_layout()
plt.show()
