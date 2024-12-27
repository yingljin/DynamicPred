# This script subsets the NHANES data 
# so that is can be stored on GitHub

library(tidyverse)

df <- read_rds(here("DataRaw/nhanes_bi.rds"))

# take 1000 subjects
subset_id <- sample(unique(df$SEQN), size = 1000)

df <- df %>% filter(SEQN %in% subset_id)

saveRDS(df, here("DataRaw/nhanes_bi_sub.rds"))
