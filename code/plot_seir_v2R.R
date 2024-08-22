seir_v2R_main_plot <- function(stan_fit = res){

N <- 8086929
seir_sum = stan_fit$summary(c("S_u", "S_v1", "S_v2", "E_u", "E_v1", "E_v2", 
                              "I_u", "I_v1", "I_v2", "R_u", "R_v1", "R_v2"))
seir_sum = seir_sum %>% select(variable, median, q5, q95) %>% 
  mutate(comp = factor(substr(variable, 1, 1), levels = c("S", "E", "I", "R")),
         vacc_stat = factor(substr(variable, 3, 3), levels = c("u", "v")),
         vacc_stat_2 = factor(substr(variable, 3, 4)),
         frac = median/N,
         day = str_extract(variable, "\\d+,") %>% str_remove(",") %>% as.numeric() + ymd("2021-01-01")-1) 


seir_vacc_stat_sum = stan_fit$summary(c("R_u_v1")) 
seir_vacc_stat_sum_1 = seir_vacc_stat_sum %>% mutate(comp = factor(substr(variable, 1, 1), levels = c("S", "E", "I", "R")),
                                                     vacc_stat = factor(substr(variable, 3, 3), levels = c("u", "v")),
                                                     vacc_stat_2 = factor(substr(variable, 3, 4)),
                                                     frac = median/N,
                                                     day = str_extract(variable, "\\d+,") %>% str_remove(",") %>% as.numeric() + ymd("2021-01-01")-1) %>% 
  group_by(day,comp) %>% 
  summarise(frac = sum(median, na.rm = T)/(N)) %>% 
  as.data.frame() %>% 
  mutate(f = cumsum(frac)) 

seir_res = seir_sum %>% 
  group_by(day,vacc_stat,comp) %>% 
  summarise(frac = sum(median)/(N)) %>% 
  ggplot() +
  geom_area(aes(day, frac, fill = vacc_stat)) +
   geom_area(data = seir_vacc_stat_sum_1 %>% 
               mutate(vacc_stat = "v1"), aes(day, f,fill = vacc_stat)) +
  facet_wrap(~comp, scales = "free_y", ncol = 4) +
  theme(legend.position = "bottom") +
  ylab("Pop. share")+ 
  xlab("Date") +
  scale_x_date(breaks = ymd(c("2021-01-01", "2021-07-01", "2021-12-01")), 
               date_labels ="%b") +  #%Y") + 
  theme(legend.title = element_text()) +
  guides(fill = guide_legend(title = "Vaccination status")) +
  scale_fill_manual(values =c('#42B540FF',cols[5], cols[4]), 
                    labels = c('Unvaccinated',"Vaccinated prior infection", "Vaccinated post infection"))


seir_res# +
 # theme(legend.margin = 0)

}