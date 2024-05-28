
# This script saves the code for some additional figures
# produced for presentation or reports

#### NHANES example ####

# 2011-2014 binary

df_nhanes <- read_rds(here("Data/nhanes_bi.rds"))
df_nhanes %>% 
  select(SEQN, Z, sind) %>% 
  filter(SEQN %in% sample(unique(df_nhanes), size = 100)) %>% head()
  pivot_wider(names_from = sind, values_from = Z) %>% 
  select(-SEQN) %>% distinct(.)


pred_appl_gfosr_l5 %>% 
  filter(id %in% rand_id) %>%
  select(id, sind, pred360) %>%
  pivot_wider(names_from = sind, values_from = pred360) %>% 
  select(-id) %>% distinct(.)