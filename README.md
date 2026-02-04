## Code and derivatives accompanying the paper

**Spatiotemporal mapping of brain organisation following the administration of 2C-B and psilocybin**  
*Molecular Psychiatry* 
https://www.nature.com/articles/s41380-026-03447-0

## Data availability

This repository does **not** contain raw neuroimaging or behavioural data. This work was finished almost a year and a half ago, but much of that time was spent waiting at the editorial decision stages rather than in peer review. Glad to be able to share, but note that as a result, the example derivatives used in the paper provided under `/derivatives` need to be extended and cleared up.

To reproduce the full analyses, users will need:

- Preprocessed, parcellated resting-state fMRI time series  
- PET-derived receptor and transporter density maps (parcellated)  
- Behavioural and questionnaire data  
- Matlab/python
  
For raw data, you can contact the corresponding author to set up a DSA upon reasonable request.



## Usage notes

- Scripts are study-specific and require manual adjustment of paths.
- Analyses are provided individually rather than as a single pipeline so that they can be used in a modular way.
- The code is shared for transparency and reproducibility rather than as a general-purpose toolbox.



## What this code does

This repository contains the analysis code used in the paper to quantify drug-induced changes in resting-state functional brain organisation under 2C-B and psilocybin. The code implements the following analyses:

### Functional connectivity

- **Static functional connectivity (sFC)**  
  Computes pairwise Pearson correlations between parcellated resting-state fMRI time series, followed by Fisher z-transformation as well as its nodal derivative global FC.

- **Dynamic functional connectivity (dFC)**  
  Estimates time-varying connectivity using framewise conditional correlations and quantifies variance of connectivity fluctuations across time.


### Signal complexity

- **Sample entropy (sampEn)**  
  Computes regional BOLD signal complexity using fixed parameters (m = 2, r = 0.3 × SD).

- **Lempel–Ziv complexity (zivEn)**  
  Computes whole-brain signal complexity from binarised time series.

- **Degree-distribution entropy (degreeEn)**  
  Quantifies entropy of graph degree distributions as a complementary complexity metric.

- **Dynamic connectivity entropy (dccEn)**  
  Computes Shannon entropy of dynamic conditional correlation time series, providing an edge-wise measure of variability in functional connectivity dynamics.

### Multivariate integration

- **Regional multivariate coherence**  
  Integrates static connectivity, dynamic connectivity variance, and entropy measures to quantify how consistently regions respond across metrics.

### Statistical inference

- **Network-Based Statistics (NBS)**  
  Identifies connected components of edges showing significant drug effects in a repeated-measures design.

- **Regional statistical models**  
  Performs linear mixed-effects analyses and post-hoc contrasts on nodal measures with multiple-comparison correction.

- **Multilevel partial least squares (PLS)**  
  Relates drug-induced changes in brain organisation to subjective experience measures while accounting for within-subject drug contrasts.

### Pharmacology-informed mapping

- **Dominance analysis with PET maps**  
  Quantifies the relative contribution of PET-derived receptor and transporter density maps to the spatial distribution of drug effects.


## Additional dependencies
- `NBS` — https://www.nitrc.org/projects/nbs/
- `pyls` — https://github.com/rmarkello/pyls  
- `COPBET` — https://github.com/anders-s-olsen/CopBET  
- `rotate_parcellation` — https://github.com/frantisekvasa/rotate_parcellation  


## Citation

If you use this code, please cite:

Mallaroni, P., Singleton, S.P., Mason, N.L. et al. 
Spatiotemporal mapping of brain organisation following the administration of 2C-B and psilocybin. Mol Psychiatry (2026). 
https://doi.org/10.1038/s41380-026-03447-0


---

## Contact

For questions regarding the code or analyses:

**Pablo Mallaroni**  
p.mallaroni@ic.ac.uk  
p.mallaroni@maastrichtuniversity.nl
