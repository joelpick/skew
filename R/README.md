# 00_functions.R

- Contains the functions that are commonly used in other scripts.


# 00_analysis_options.R

- Input analysis options which are used for the downstream analysis, such as whether the analysis is using a reduced set of fixed effects and what variables are included in which analyses.
- Produces "analysis_options.Rdata"


# 01_data_prep.R

- Extract day 15 body size and associated data
- Produces "Data/Intermediate/chick_data.Rdata" which contains data.frames:
    - 'TBHW' - all day 15 data
    - 'TBHW_noRep' - day 15 data without repeated measures
    - 'ped' - pedigree
 

# 02_stan_data.R

- Makes data for trait models in stan
- Produces "Data/Intermediate/stan_data_reduced.Rdata" and "Data/Intermediate/stan_data.Rdata" with reduced and full fixed effect structures respectively. These contain the lists:
    - 'stan_data_DS' - data for stan skew t dam sire analysis with repeated measures
    - 'stan_data_DS_noRep' - data for stan skew t dam sire analysis without repeated measures
    - 'stan_data_ped' - data for stan skew t animal model analysis with repeated measures
	- 'stan_data_ped_reduced' - data for stan gaussian reduced animal model analysis with repeated measures
    - 'stan_data_ped_noRep' - data for stan skew t animal model analysis without repeated measures
    - 'stan_data_ped_noRep' - data for stan skew t animal model analysis without repeated measures
    - 'stan_data_ped_reduced_noRep' - data for stan gaussian reduced animal model analysis without repeated measures


# 03_generate_starting_values.R

- Generates starting values for trait models in stan using asreml
- Runs similarly structured asreml models assuming normality
- Produces "Data/Intermediate/starting_values.Rdata" and "Data/Intermediate/starting_values_reduced.Rdata", which contains lists:
    - modT - tarsus asreml model results
    - modHB - head-bill asreml model results
    - modW - wing asreml model results
    - modM - weight asreml model results


# 04_damSire_analysis.R

- Runs dam-sire multi-membership models with skew-t dam/sire, between and within nest effects and normal residuals for models with repeated measures. 
- Also has sensitivity analyses for wing length models with alternate priors
- For each of the four traits this produces
    - Data/Intermediate/stanMod_DS_reduced_trait_time.Rdata
    - Data/Intermediate/stanMod_DS_reduced_trait_alpha10_time.Rdata
    - Data/Intermediate/stanMod_DS_reduced_trait_alpha1_time.Rdata
  which all contain 'model_z' containing the posterior distributions of the traits models


# 05_animalModel_analysis.R

- Runs animal models with normal breeding values and skew t between and within nest effects and normal residuals for models with repeated measures. 
- Also runs Gaussian animal models.
- For each of the four traits this produces
    - Data/Intermediate/stanMod_pedN_reduced_trait_time.Rdata
    - Data/Intermediate/stanModNormal2_pedN_reduced_trait_time.Rdata
  which all contain 'model_z' containing the posterior distributions of the traits models


# 06_survival_analysis.R

- Runs survival models for all traits accounting for measurement error
- Produces
    - Data/Intermediate/day15_survival_ME_trait_time.Rdata


# 07_selection_gradients.R

- Runs selection gradient code
- Produces
    - selection_gradients_ME_trait.Rdata
    - internal_optimum_ME.Rdata
- Plots
	- figure4_ME.pdf


# 08_h2.R
- Runs heritability code
- Produces
    - h2_trait.Rdata


# 09_POreg.R
- Calculates non-linear parent offspring regression
- Produces
    - nonLinearPO_trait.Rdata
    - h2_results.Rdata
- Plots 
    - figure5.pdf


# 10_sim_ped.R
- Runs simulations to test the impact of skew and cross fostering on the estimation of Va and h2
- Produces
    - sim_ped_weight_g_N_XR.Rdata
    - sim_ped_weight_g_ST_XR.Rdata
    - sim_ped_weight_g_NR.Rdata
    - sim_ped_weight_g_NR.Rdata
- Plots
    - figure6.pdf


# 11_MA_data.R
- Collates data for the meta-analysis of skew across species
- Produces
    - MA_data.csv


# 12_meta_analysis.R
- Runs meta-analysis on across species skew data
- Produces
    - meta_analysis.Rdata
- Plots
    - figure2.pdf


# 13_response.R
- Calcaultes response to selection
- Produces
    - response.Rdata
- Plots 
    - figure_SM_heywood.pdf


# figure1.R
- Code for producing figure1.pdf


# figure3.R
- Code for producing figure3.pdf and figure_SM_skewt_Gaussian.pdf


# SM_dam_sire_egg_analysis.R
- Runs gaussian models of blue tit data in asreml to compare dam-sire and animal models, models with and without egg size, and models with two different pedigrees
- Produces
	- dam_sire_egg.Rdata (model results with polygamous pedigree)
    - dam_sire_egg_M.Rdata  (model results with monogamous pedigree)


# SM_dam_sire_egg_plots.R
- Plots results from SM_dam_sire_egg_analysis.R and produces 
    - figure_SM_eggSize.pdf, 
    - figure_SM_DS_Animal.pdf
    - figure_SM_pedigree.pdf


# SM_fix_eff_sim.R
- Simulates data based on wing length dam-sire model, and analyses it with a Gaussian LMM
- Produces
    - lmer_sims.Rdata
- Plots
    - figure_SM_lmer_sim.pdf


# SM_heterosked_sim.R
- Simulations f different scenarios, including skewed data, heteroskedasticly and unbalanced level, analysed with gaussian, skew, skew t and t models.
- Produces
    - MP_sims_N.Rdata
    - MP_sims_S.Rdata
    - MP_sims_S_hom.Rdata
    - MP_sims_S_homM.Rdata
    - MP_sims_S_homN.Rdata
    - MP_sims_S_P.Rdata
- Plots
    - figure_SM_fixed_sim.pdf


# SM_post_check.R
- Posterior predictive checks of Gaussian and skew t models for all traits.
- Plots
    - figure_SM_post_check.pdf


# SM_prior_post.R
- Plots posterior distributions and priors of dam-sire models for all traits
    - figure_SM_prior_post_trait.pdf


# SM_prior_sensitivity.R
- Plots comparison of dam-sire models of wing length with different priors. 
    - figure_SM_prior_sensitivity1.pdf
    - figure_SM_prior_sensitivity2.pdf


# SM_priors.R
- Plots induced priors on gammas given different priors on alpha and delta
    - figure_SM_skew_priors.pdf


# SM_tables.R
- Produces results tables, saved as
    - table_data.Rdata


# SM_wing_skew.R
- Plots the disparity in fixed effect estimation between skew t and Gaussian models of wing length 
    - figure_SM_wing_skew.pdf

