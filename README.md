# Orthogonal Policy Learning with Ordinal Outcomes

Code accompanying the paper *Orthogonal Policy Learning with Ordinal Outcomes* by Yue Zhang, Shanshan Luo, and Yangbo He.

## Files

- `data_gen.R`: Simulation settings, data generation, and shared helper functions.
- `sim_Rsup.R`: Worst-case regret simulation.
- `sim_Rsup_tab.R`: Tables and figures for the regret simulation.
- `sim_trPr.R`: Utility-threshold sensitivity simulation.
- `sim_trPr_plot.R`: Treatment-proportion figures for the sensitivity simulation.
- `generate_sipp_december.R`: Optional preparation of the original SIPP data.
- `sipp.r`: SIPP bootstrap analysis and fitted policy trees.
- `sipp_plot.r`: SIPP figures, workbook, and summary CSV.
- `sipp_data_dec.rds`: Prepared December sample containing 8,319 distinct respondents.

## Usage

Keep the eight R scripts in the same directory. Before running them, fill in the blank paths marked `"USER PATH"` and install the R packages required by the scripts.

Run the scripts in the following order:

1. **Regret simulation:** `data_gen.R` → `sim_Rsup.R` → `sim_Rsup_tab.R`
2. **Sensitivity simulation:** `data_gen.R` → `sim_trPr.R` → `sim_trPr_plot.R`
3. **SIPP application:** `sipp.r` → `sipp_plot.r`

The included SIPP data can be used directly; regeneration is optional.
