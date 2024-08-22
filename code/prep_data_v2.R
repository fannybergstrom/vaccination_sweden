seir_prep_dat_v2 <- function(date_end = "2021-12-31",
                             theta = 1 / 4.6,
                             gamma = 1 / 2.1,
                             ve_inf = 0.5,
                             ve_susc1 = 0.5,
                             ve_susc2 = 0.8,
                             n_age_groups = 8) {
  
  length_time <- (date_end %>% as.Date() - "2021-01-01" %>% as.Date()) %>% as.numeric() + 1
  
  # Fraction observed
  frac_obs <- readRDS(str_c(here::here(), "/data/processed/frac_rep_age_modfree.rds")) %>%
    head(length_time)
  
  drac_frac_2020 <- 1 / readRDS(str_c(here::here(), "/data/processed/frac_rep_age_modfree_drakfrac.rds")) %>%
    select(-date)
  
  
  # Case data
  cases <- read_rds(file.path(here::here(), "data", "fohm", "case_daily_age_import.Rdata")) %>%
    select(-contains("age_0_9"), -contains("age_10_19"))
  
  # 7-day average of cases
  cases <- cases %>%
    mutate_if(is.numeric, roll_mean, n = 14, fill = NA, align = "right", na.rm = T) %>%
    filter(date <= date_end) %>%
    na.omit()
  
  cases_2020 <- cases %>% filter(date < ymd("2021-01-01")) 
  
  cases_mod <- cases %>% filter(date >= ymd("2021-01-01"), date <= ymd(date_end))
  
  # Intital values SEIR
  
  new_inf_1 <- cases_2020 %>%
    slice_tail(n = 1) %>%
    select(contains("change_local")) %>%
    slice_tail(n = 1) %>%
    as.matrix() /
    (frac_obs %>% slice_tail(n = 1))
  
  E_u_1 <- ((cases_2020 %>% slice_tail(n = 26) %>% select(contains("change_local")) * (1 - theta)^(25:0) * (1 - gamma)^(25:0)) %>% as.matrix() *
              (drac_frac_2020 %>% slice_tail(n = 26))) %>%
    apply(2, sum)
  
  I_u_1 <- (((cases_2020 %>% slice_tail(n = 26) %>% select(contains("change_local")) * (1 - gamma)^(25:0)) %>% as.matrix() *
               (drac_frac_2020 %>% slice_tail(n = 26))) %>%
              apply(2, sum) - E_u_1) %>% round()
  
  R_u_1 <- (cases_2020 %>% select(contains("change_local")) %>% as.matrix() * drac_frac_2020) %>%
    apply(2, sum) %>%
    round() - E_u_1 - I_u_1
  
  # Vaccination
  v <- readRDS(file.path(here::here(), "data", "processed", "vacc_frac_swe_age.rds")) %>%
    filter(date >= ymd("2021-01-01"), date <= ymd(date_end)) %>%
    select(-contains("age_0_9"), -contains("age_10_19")) %>%
    select(contains("vacc"))
  
  vacc <- readRDS(file.path(here::here(), "data", "processed", "vacc_data_daily_swe_age.rds")) %>%
    select(-contains("age_0_9"), -contains("age_10_19"))
  
  vacc$date.x <- seq(as.Date("2020-12-05"), as.Date("2022-04-10"), by = 1)
  
  # Contact matrix
  contact <- readRDS(str_c(here::here(), "/data/processed/contact_matrix.rds"))
  
  # Age distribution Sweden
  N_a <- read_excel("data/age_dist_swe.xlsx") %>%
    na.omit() %>%
    pull(N)
  N_a <- N_a[2:9]
  
  # Return data list
  list(
    n_days = as.numeric(ymd(date_end) - ymd("2021-01-01")) + 1,
    N_a = N_a,
    N = sum(N_a),
    n_a = N_a %>% length(),
    rep_cases_local = cases_mod %>% filter(date <= date_end) %>% 
      select(contains("change_local")) %>% as.matrix() %>% round(),
    rep_cases_import_u = cases_mod %>% filter(date <= date_end) %>%
      select(contains("change_import")) %>% as.matrix() *
      (1 - v) / ((1 - v) + v * (1 - ve_inf)),
    rep_cases_import_v1 = cases_mod %>% filter(date <= date_end) %>%
      select(contains("change_import")) %>% as.matrix() *
      v * (1 - ve_inf) / ((1 - v) + v * (1 - ve_inf)) ,
    rep_cases_import_v2 = cases_mod %>% filter(date <= date_end) %>%
      select(contains("change_import")) %>% as.matrix() *
      v * (1 - ve_inf) / ((1 - v) + v * (1 - ve_inf)),
    n_vac1 = vacc %>% filter(date.x >= as.Date("2021-01-01") + 20, date.x <= as.Date(date_end) + 20) %>% select(contains("vacc")) %>% as.matrix(),
    n_vac2 = vacc %>% filter(date.x >= as.Date("2021-01-01") - 10, date.x <= as.Date(date_end) - 10) %>% select(contains("vacc")) %>% as.matrix(),
    E_u_1 = E_u_1 %>% unlist(),
    I_u_1 = I_u_1 %>% unlist(),
    R_u_1 = R_u_1 %>% unlist(),
    new_inf_1 = new_inf_1 %>% unlist(),
    c = contact,
    theta = theta,
    gamma = gamma,
    ve_inf = ve_inf,
    ve_susc1 = ve_susc1,
    ve_susc2 = ve_susc2,
    frac_obs = frac_obs
  )
}
