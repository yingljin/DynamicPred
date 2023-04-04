# data application of NHANES

set.seed(404)

library(here)
library(tidyverse)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)

df <- read_rds(here("Data/nhanes_bi.rds"))

N <- length(unique(df$SEQN)) # sample size 8763
J <- max(df$sind) # 1440 measures for each subject

# visulization
rand_id <- sample(unique(df$SEQN), size = 4)

df %>% 
  filter(SEQN %in% rand_id) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Z))+
  facet_wrap(~SEQN)
    
