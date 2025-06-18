# GBM-PTEN Model Codes

This repository contains code and results for the multiscale modeling framework used in the research article:

> **"Multiscale modeling of tumor-macrophage interactions in PTEN-Null glioblastoma predicts synergistic combination therapies"**

The modeling framework includes:
- A **Mean-field PDEs-ODEs model** that simulates tumor growth, signaling pathways, and drug responses at the population level.
- A **Multiscale agent-based model (MSABM)** that captures tumor-immune interactions and therapy effects at the single-cell level.

---

## Directory Structure

### `MSABM/`

MSABM for PTEN-Null and PTEN-WT glioblastoma.

- `PTEN_WT_fit_codes & results`: Code and results for fitting tumor growth in PTEN-WT models.
- `PTEN_combination_therapies_codes & results`: Simulations of drug combination treatments under PTEN-Null and PTEN-WT gliomas.
- `PTEN_growth_fit_codes & results`: Tumor growth fitting under PTEN-Null and PTEN-WT conditions.
- `PTEN_synergy_codes & results`: Synergy analysis for drug combination treatments using MSABM simulations.

### `Mean-field PDEs-ODEs model/`

Mean-field simulations of tumor and immune signaling interactions.

- `1_signaling_pathway_fit`: Fitting of intracellular signaling pathways.
- `2_glioma_growth_volume_fit`: PTEN-Null and PTEN-WT glioma growth model calibration using volume data.
- `3_sensitivity_analysis`: Global and local sensitivity analysis of model parameters.
- `4_KM_curves_fit`: Fitting model predictions to Kaplan-Meier survival curves.
- `5_KM_curves_plot`: Code for plotting survival curves under various treatments.
- `6_CSF1RI_only_as_control`: Model outputs for CSF1R inhibitor monotherapy.
- `7_drug_combination`: Simulations of dual-drug combination treatments.

---

For questions or collaboration inquiries, please contact the corresponding author listed in the paper.
