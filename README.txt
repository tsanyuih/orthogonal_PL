Orthogonal Policy Learning with Ordinal Outcomes
================================================
Code accompanying the manuscript by Yue Zhang, Shanshan Luo, and Yangbo He.

Files
-----
data_gen.R                Simulation data, settings, and shared helpers.
sim_Rsup.R                Worst-case regret simulation.
sim_Rsup_tab.R            Simulation regret tables and figures.
sim_trPr.R                Utility-threshold sensitivity simulation.
sim_trPr_plot.R           Simulation treatment-proportion figures.
generate_sipp_december.R  Optional preparation from the original SIPP data.
sipp.r                    SIPP bootstrap analysis and fitted policy trees.
sipp_plot.r               SIPP figures, workbook, and summary CSV.
sipp_data_dec.rds         Prepared December sample of 8,319 distinct respondents.

Usage
-----
Keep the eight R scripts together. Before running, fill in the blank paths
marked "USER PATH" and install the R packages required by the scripts.
Run the scripts in this order:
  Regret simulation: data_gen.R -> sim_Rsup.R -> sim_Rsup_tab.R
  Sensitivity:       data_gen.R -> sim_trPr.R -> sim_trPr_plot.R
  SIPP application:  sipp.r -> sipp_plot.r
The included SIPP data can be used directly; regeneration is optional.
