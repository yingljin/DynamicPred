
# This script includes a motivating example created for the manuscript
# model after leroux2017 
# use NHANES binary indicator

library(here)
library(tidyverse)
library(ggpubr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(kableExtra)
library(mvtnorm)
library(mgcv)
library(splines)
library(LaplacesDemon)
library(arsenal)
library(RColorBrewer)
theme_set(theme_minimal())

#### data ####

df <- read_rds(here("Data/nhanes_bi.rds"))
load(here("Data/ApplOutput_fGFPCA.RData"))
length(unique(df$SEQN)) # 8763 subjects
df %>% select(SEQN, age_years_interview, BMI, gender) %>% distinct(.) %>% 
  select(BMI) %>% 
  lapply(mean)

# age: 49.53 (17.52)
# BMI: 29.08 (7.04)
# gender: 4202 male, 4561 female


#### Figure ####

# illustration of dynamic prediction 
# one male (62161) and one female (62164)
pred_nhanes_fgfpca %>% filter(id==62161 | id == 62164) %>% 
  mutate(Y = ifelse(sind<=720, Y, NA)) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y))+
  geom_line(aes(x=sind, y = exp(pred720)/(1+exp(pred720))))+
  facet_wrap(~id)

df_plot <- pred_nhanes_fgfpca %>% filter(id == 62164) %>% 
  pivot_longer(starts_with("pred"))
df_plot$Y[df_plot$name=="pred360" & df_plot$sind>360] <- NA
df_plot$Y[df_plot$name=="pred720" & df_plot$sind>720] <- NA
df_plot$Y[df_plot$name=="pred1080" & df_plot$sind>1080] <- NA


# plot 
display.brewer.pal(4, "Set2")
brewer.pal(5, "Set2")

df_plot %>% 
  mutate(value = exp(value)/(1+exp(value))) %>%
  mutate(name=as.numeric(gsub("pred", "", name))) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.5)+
  geom_line(aes(x=sind, y = value), col = "red")+
  geom_vline(aes(xintercept = name), linetype = "dashed")+
  facet_grid(cols = vars(name))+
  scale_x_continuous(breaks = c(360, 720, 1080), 
                     labels = c("6am", "12pm", "6pm"))+
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend = element_blank())+
  labs(x="", y = "Y")
ggsave(here("Images/MotiveExp.jpeg"),
       width=10, height=3)
