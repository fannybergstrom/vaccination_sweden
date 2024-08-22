// Age structured SEIR model. Two vaccination groups. Counting removed chaning vaccination status and direct and indirect effect.

data {
  // Selected parameters of SEIR model and data for fitting
  // time and population data 
  int<lower=1> n_days; // number of days
  int<lower=1> N; // population size
  int<lower=1> n_a; // number age groups
  row_vector[n_a] N_a; // population size per age group
  matrix[n_a, n_a] c; //contact matrix with entries c[a',a]: 
  // "av. nuber of contacts an individual from age-group a' has with age-group a"
  
  // SARS-Cov-2 related data
  matrix[n_days, n_a] n_vac1; // No. vaccinated w/ 1st dose
  matrix[n_days, n_a] n_vac2; // No. vaccinated w/ 2nd dose
  array[n_days, n_a] int<lower=0> rep_cases_local; // reported cases due to local infection
  matrix[n_days, n_a] rep_cases_import_u; // reported cases, imported 
  matrix[n_days, n_a] rep_cases_import_v1; 
  matrix[n_days, n_a] rep_cases_import_v2;
  
  // Initialization 
  row_vector[n_a] E_u_1; // exposed at day 1    
  row_vector[n_a] I_u_1; // infectious at day 1                                                                                                                                                         
  row_vector[n_a] R_u_1; // recovered at day 1
  row_vector[n_a] new_inf_1; // newly infected at day 1                                                                                                                                                
  
  // Fixed Parameters of SEIR model
  real<lower=0> theta; // rate Exposed->Infectious SEIR model  
  real<lower=0> gamma; // recovery rate SEIR model 
  
  // vaccine efficacy
  real<lower=0, upper=1> ve_inf;   // ve against infectiousness
  real<lower=0, upper=1> ve_susc1; // ve against susceptibility 1st dose
  real<lower=0, upper=1> ve_susc2; // ve against susceptibility 2nd dose
  
  // time- and age-varying probability of reporting newly infected as case
  matrix[n_days, n_a] frac_obs;
}

parameters {
  vector[n_days] log_beta;
  
  real<lower=0> sd_beta;
}

transformed parameters {
  // Compartments
  matrix[n_days, n_a] S_u; // Unvaccinated susceptibles
  matrix[n_days, n_a] S_v1; // Vaccinated susceptibles
  matrix[n_days, n_a] S_v2; // Vaccinated susceptibles
  matrix[n_days, n_a] E_u; // Unvaccinated exposed
  matrix[n_days, n_a] E_v1; // Vaccinated exposed
  matrix[n_days, n_a] E_v2; // Vaccinated exposed
  matrix[n_days, n_a] I_u; // Unvaccinated infectious
  matrix[n_days, n_a] I_v1; // Vaccinated infectious 
  matrix[n_days, n_a] I_v2; // Vaccinated infectious 
  matrix[n_days, n_a] R_u; // Unvaccinated recovered
  matrix[n_days, n_a] R_v1; // Vaccinated recovered 
  matrix[n_days, n_a] R_v2; // Vaccinated recovered 
  
  // Helper compartments
  matrix[n_days, n_a] new_inf; // New infections
  matrix[n_days, n_a] new_imp_u; // New imported infections unvaccinated   
  matrix[n_days, n_a] new_imp_v1; // New imported infections vaccinated   
  matrix[n_days, n_a] new_imp_v2; // New imported infections vaccinated   
  
  // Time-varying transmission rate
  vector<lower=0>[n_days] beta;
  
  // Initialization
  row_vector[n_a] Init_v_1;
  matrix[n_days, n_a] Init_m;
  // Vaccinated at time 1 ()
  Init_v_1 = rep_row_vector(0, n_a);
  // Exposed at time 1
  Init_m = rep_matrix(0.0000001, n_days, n_a);
  // Exposed 
  // Unvaccinated
  E_u = Init_m;
  E_u[1,  : ] = E_u_1; // Initialize based on Input data                                                                                                                                                                               
  
  E_v1 = Init_m; // No vaccinations at time 1
  E_v2 = Init_m; // No vaccinations at time 1
  
  // Infectious
  I_u = Init_m;
  I_u[1,  : ] = I_u_1; // Initialize based on Input data      
  I_v1 = Init_m; // No vaccinations at time 1
  I_v2 = Init_m; // No vaccinations at time 1
  
  // Susceptible
  S_u = Init_m;
  S_u[1,  : ] = N_a - R_u_1 - I_u_1 - E_u_1; // Susceptible = N-I-E-R, Initialize based in Input data
  S_v1 = Init_m; // No vaccinations at time 1
  S_v2 = Init_m; // No vaccinations at time 1
  
  // Recovered
  R_u = Init_m;
  R_u[1,  : ] = R_u_1; // Initialize based on Input data 
  R_v1 = Init_m; // No vaccinations at time 1
  R_v2 = Init_m; // No vaccinations at time 1
  
  // Helper compartment newly infected
  new_inf = Init_m;
  new_inf[1,  : ] = new_inf_1; // Initialize based on Input data
  
  // Imported cases per day and age-group (input data, upscaled by inverse of assumed reporting fraction)
  for (t in 1 : n_days) {
    for (a in 1 : n_a) {
      new_imp_u[t, a] = 1 / frac_obs[t, a] * rep_cases_import_u[t, a];
      new_imp_v1[t, a] = 1 / frac_obs[t, a] * rep_cases_import_v1[t, a];
      new_imp_v2[t, a] = 1 / frac_obs[t, a] * rep_cases_import_v2[t, a];
    }
  }
  
  // The SEIR model
  beta = exp(log_beta);
  
  for (t in 2 : n_days) {
    for (a in 1 : n_a) {
      // Unv. susceptibles at time t: Susceptibles at t-1 minus vaccinations among susceptibles   
      // minus new infections
      
      S_u[t, a] = S_u[t - 1, a]
                  - n_vac1[t - 1, a] * S_u[t - 1, a]
                    / (S_u[t - 1, a] + R_u[t - 1, a])
                  - beta[t]
                    * sum(c[ : , a]'
                          .* (I_u[t - 1,  : ]
                              + (1 - ve_inf)
                                * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                    * S_u[t - 1, a] / N_a[a];
      // 1-dose vacc susceptibles at time t: 1-dose vacc susceptibles at t-1 plus 
      // vaccination among unvacinated susceptibles at t-1 minus 
      // new infections mins
      // vaccinations among 1-dose susceptibles   
      
      S_v1[t, a] = S_v1[t - 1, a]
                   + n_vac1[t - 1, a] * S_u[t - 1, a]
                     / (S_u[t - 1, a] + R_u[t - 1, a])
                   - beta[t]
                     * sum(c[ : , a]'
                           .* (I_u[t - 1,  : ]
                               + (1 - ve_inf)
                                 * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                     * S_v1[t - 1, a] / N_a[a] * (1 - ve_susc1)
                   - n_vac2[t - 1, a] * S_v1[t - 1, a]
                     / (S_v1[t - 1, a] + R_v1[t - 1, a] + 0.00000001);
      // 2-dose vacc susceptibles at time t: 2-dose vacc susceptibles at t-1 plus 
      // vaccinations among 1-dose susceptibles at t-1 minus 
      // new infections  
      S_v2[t, a] = S_v2[t - 1, a]
                   + n_vac2[t - 1, a] * S_v1[t - 1, a]
                     / (S_v1[t - 1, a] + R_v1[t - 1, a] + 0.00000001)
                   - beta[t]
                     * sum(c[ : , a]'
                           .* (I_u[t - 1,  : ]
                               + (1 - ve_inf)
                                 * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                     * S_v2[t - 1, a] / N_a[a] * (1 - ve_susc2);
      
      // Exposed at time t: Exposed at t-1 + new infected, minus new infectious   
      E_u[t, a] = E_u[t - 1, a]
                  + beta[t]
                    * sum(c[ : , a]'
                          .* (I_u[t - 1,  : ]
                              + (1 - ve_inf)
                                * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                    * S_u[t - 1, a] / N_a[a]
                  - theta * E_u[t - 1, a];
      
      E_v1[t, a] = E_v1[t - 1, a]
                   + beta[t]
                     * sum(c[ : , a]'
                           .* (I_u[t - 1,  : ]
                               + (1 - ve_inf)
                                 * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                     * S_v1[t - 1, a] / N_a[a] * (1 - ve_susc1)
                   - theta * E_v1[t - 1, a];
      
      E_v2[t, a] = E_v2[t - 1, a]
                   + beta[t]
                     * sum(c[ : , a]'
                           .* (I_u[t - 1,  : ]
                               + (1 - ve_inf)
                                 * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                     * S_v2[t - 1, a] / N_a[a] * (1 - ve_susc2)
                   - theta * E_v2[t - 1, a];
      
      // Infectious at time t: Infectious at t-1 + new infectious, minus recovered, plus \"imported\" infectious                                                                                            
      
      I_u[t, a] = I_u[t - 1, a] - gamma * I_u[t - 1, a]
                  + theta * E_u[t - 1, a] + new_imp_u[t, a];
      
      I_v1[t, a] = I_v1[t - 1, a] - gamma * I_v1[t - 1, a]
                   + theta * E_v1[t - 1, a] + new_imp_v1[t, a];
      
      I_v2[t, a] = I_v2[t - 1, a] - gamma * I_v2[t - 1, a]
                   + theta * E_v2[t - 1, a] + new_imp_v2[t, a];
      
      // Recovered at t: Recovered at t-1 + new recovered, minus/plus vaccinations among recovered                                                                                                             
      
      R_u[t, a] = R_u[t - 1, a] + gamma * I_u[t - 1, a]
                  - (n_vac1[t - 1, a] * R_u[t - 1, a]
                     / (S_u[t - 1, a] + R_u[t - 1, a] + 0.0000001));
      
      R_v1[t, a] = R_v1[t - 1, a] + gamma * I_v1[t - 1, a]
                   + (n_vac1[t - 1, a] * R_u[t - 1, a]
                      / (S_u[t - 1, a] + R_u[t - 1, a]))
                   - (n_vac2[t - 1, a] * R_v1[t - 1, a]
                      / (S_v1[t - 1, a] + R_v1[t - 1, a]) + 0.000001);
      
      R_v2[t, a] = R_v2[t - 1, a] + gamma * I_v2[t - 1, a]
                   + n_vac2[t - 1, a] * R_v1[t - 1, a]
                     / (S_v1[t - 1, a] + R_v1[t - 1, a] + 0.0000001);
      
      // newly infected
      new_inf[t, a] = beta[t]
                      * (sum(c[ : , a]'
                             .* (I_u[t - 1,  : ]
                                 + (1 - ve_inf)
                                   * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                         * S_u[t - 1, a] / N_a[a]
                         + sum(c[ : , a]'
                               .* (I_u[t - 1,  : ]
                                   + (1 - ve_inf)
                                     * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                           * S_v1[t - 1, a] / N_a[a] * (1 - ve_susc1)
                         + sum(c[ : , a]'
                               .* (I_u[t - 1,  : ]
                                   + (1 - ve_inf)
                                     * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                           * S_v2[t - 1, a] / N_a[a] * (1 - ve_susc2));
    }
  }
}
model {
  // First-order random walk for log-beta
  sd_beta ~ normal(1, 0.2);
  log_beta[1] ~ normal(0.66, 0.3);
  
  for (t in 2 : n_days) {
    target += normal_lpdf(log_beta[t] | log_beta[t - 1], sd_beta);
  }
  // Poisson-Likelihood observed cases
  for (t in 1 : n_days) {
    for (a in 1 : n_a) {
      target += poisson_lpmf(rep_cases_local[t, a] | frac_obs[t, a]
                                                     * new_inf[t, a]);
    }
  }
}
generated quantities {
  matrix[n_days, n_a] R_e;
  
  matrix[n_days, n_a] R_e_count;
  
  matrix[n_days, n_a] S_count;
  
  matrix[n_days, n_a] E_count;
  
  matrix[n_days, n_a] I_count;
  
  matrix[n_days, n_a] R_count;
  
  matrix[n_days, n_a] new_inf_count;
  
  matrix[n_days, n_a] new_inf_u;
  
  matrix[n_days, n_a] new_inf_v1;
  
  matrix[n_days, n_a] new_inf_v2;
  
  matrix[n_days, n_a] S_dir; // Unvaccinated susceptibles                                                                                                                                                 
  
  matrix[n_days, n_a] E_dir; // Unvaccinated exposed                                                                                                                                                   
  
  matrix[n_days, n_a] I_dir; // Vaccinated infectious                                                                                                                                                      
  
  matrix[n_days, n_a] R_dir; // Unvaccinated recovered                                                                                                                                                                               
  
  matrix[n_days, n_a] new_inf_dir; // New infections
  
  matrix[n_days, n_a] R_u_v1; // Removed moving from unvaccinated to vaccinated
  
  matrix[n_days, n_a] R_v1_v2; // Removed moving from one dose to two doses
  
  row_vector[n_a] diff_inf;
  
  E_count[1,  : ] = E_u_1; //                                                                                                                                                                               
  
  I_count[1,  : ] = I_u_1; //   
  
  R_count[1,  : ] = R_u_1; //
  
  S_count[1,  : ] = N_a - R_u_1 - I_u_1 - E_u_1;
  
  new_inf_count[1,  : ] = new_inf_1;
  
  new_inf_u[1,  : ] = new_inf_1;
  
  new_inf_v1[1,  : ] = Init_v_1;
  
  new_inf_v2[1,  : ] = Init_v_1;
  
  new_inf_dir[1,  : ] = Init_v_1;
  
  R_u_v1[1,  : ] = Init_v_1;
  
  R_v1_v2[1,  : ] = Init_v_1;
  
  for (t in 2 : n_days) {
    for (a in 1 : n_a) {
      // Counterfactual susceptibles at time t: Susceptibles at t-1 minus new infections
      S_count[t, a] = S_count[t - 1, a]
                      - beta[t] * sum(c[ : , a]' .* I_count[t - 1,  : ])
                        * S_count[t - 1, a] / N_a[a];
      // Counterfactual exposed at t: exposed at t-1 + new infected - new infectious 
      E_count[t, a] = E_count[t - 1, a] * (1 - theta)
                      + beta[t] * sum(c[ : , a]' .* I_count[t - 1,  : ])
                        * S_count[t - 1, a] / N_a[a];
      // Counterfactual infectious at t: infectious at t-1 + new infectious + imported infectious - new recovered
      I_count[t, a] = I_count[t - 1, a] * (1 - gamma)
                      + theta * E_count[t - 1, a] + new_imp_u[t, a]
                      + new_imp_v1[t, a] + new_imp_v2[t, a];
      // Counterfactual Recovered at t: recovered at t-1 + new recovered 
      R_count[t, a] = R_count[t - 1, a] + gamma * I_count[t - 1, a];
      // New infected counterfactual
      new_inf_count[t, a] = beta[t] * sum(c[ : , a]' .* I_count[t - 1,  : ])
                            * S_count[t - 1, a] / N_a[a];
      
      new_inf_u[t, a] = beta[t]
                        * sum(c[ : , a]'
                              .* (I_u[t - 1,  : ]
                                  + (1 - ve_inf)
                                    * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                        * S_u[t - 1, a] / N_a[a];
      
      new_inf_v1[t, a] = beta[t]
                         * sum(c[ : , a]'
                               .* (I_u[t - 1,  : ]
                                   + (1 - ve_inf)
                                     * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                         * S_v1[t - 1, a] / N_a[a] * (1 - ve_susc1);
      
      // Removed changing compartment from unvaccinated to one dose
      
      R_u_v1[t, a] = (n_vac1[t - 1, a] * R_u[t - 1, a]
                      / (S_u[t - 1, a] + R_u[t - 1, a]));
      
      // Removed changing compartment from to one dose to two doses
      
      R_v1_v2[t, a] = (n_vac2[t - 1, a] * R_v1[t - 1, a]
                       / (S_v1[t - 1, a] + R_v1[t - 1, a]));
    }
  }
  
  for (a in 1 : n_a) {
    diff_inf[a] = (S_count[1, a] - S_count[n_days, a])
                  - (S_u[1, a] - S_u[n_days, a] - S_v1[n_days, a]);
  }
  
  // Estimating direct effect        
  E_dir[1,  : ] = E_u_1; //                                                                                                                                                                               
  I_dir[1,  : ] = I_u_1; //
  R_dir[1,  : ] = R_u_1; //
  S_dir[1,  : ] = N_a - R_u_1 - I_u_1 - E_u_1;
  new_inf_dir[1,  : ] = new_inf_1;
  
  for (t in 2 : n_days) {
    for (a in 1 : n_a) {
      // Susceptibles at time t: Susceptibles at t-1 minus new infections                                                                                  
      S_dir[t, a] = S_dir[t - 1, a]
                    - beta[t]
                      * sum(c[ : , a]'
                            .* (I_u[t - 1,  : ]
                                + (1 - ve_inf)
                                  * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                      * S_dir[t - 1, a] / N_a[a];
      
      // Exposed at time t: Exposed at t-1 + new infected, minus recovered - new infectious                                                                                          
      E_dir[t, a] = E_dir[t - 1, a]
                    + beta[t]
                      * sum(c[ : , a]'
                            .* (I_u[t - 1,  : ]
                                + (1 - ve_inf)
                                  * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                      * S_dir[t - 1, a] / N_a[a]
                    - theta * E_dir[t - 1, a];
      
      // Infectious at time t: Infectious at t-1 + new infectous, minus recovered, plus \"imported\" infectious                                                                                            
      I_dir[t, a] = I_dir[t - 1, a] - gamma * I_dir[t - 1, a]
                    + theta * E_dir[t - 1, a] + new_imp_u[t, a]
                    + new_imp_v1[t, a] + new_imp_v2[t, a];
      
      // Recovered at t: Recovered at t-1 + new recovered                                                                                                    
      R_dir[t, a] = R_dir[t - 1, a] + gamma * I_dir[t - 1, a];
      
      // newly infected
      new_inf_dir[t, a] = beta[t]
                          * sum(c[ : , a]'
                                .* (I_u[t - 1,  : ]
                                    + (1 - ve_inf)
                                      * (I_v1[t - 1,  : ] + I_v2[t - 1,  : ])))
                          * S_dir[t - 1, a] / N_a[a];
    }
  }
}
