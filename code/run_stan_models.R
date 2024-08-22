# Script for running SEIR model 

pacman::p_load(cmdstanr, posterior, bayesplot,
               readr, tidyverse, zoo, RcppRoll,
               here, readxl)

set_cmdstan_path(path = NULL)

## SEIR 2 vaccination groups
source(file.path(here::here(), "data_processing", "prep_data_v2.R"))

stan_mod <- "seir_cases_2vac_age"

file <- file.path(here::here(), "stan", "final", str_c(stan_mod,".stan"))
mod <- cmdstan_model(file)
mod_res <- mod$sample(data =  seir_prep_dat_v2(date_end = "2021-12-31",
                                               theta = 1/4.6,
                                               gamma = 1/2.1,
                                               ve_inf = 0.5,
                                               ve_susc1 =0.5,
                                               ve_susc2 =0.8),
                      adapt_delta = 0.95,
                      parallel_chains = 4,
                      max_treedepth = 10,
                      iter_warmup = 400,
                      iter_sampling = 400,
                      seed = 24590
)


temp_rds_file <- file.path(here::here(), "results", "stan_fit", 
                           str_c(stan_mod,"_", today(), "_a8.rds"))

mod_res$save_object(file = temp_rds_file)

