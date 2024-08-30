
# This script saves the code for figures in the manuscript
# produced for presentation or reports

#### set up ####

set.seed(1114)


library(here)
library(tidyverse)
library(ggpubr)
library(refund)
library(lme4)
library(ROCR)
library(GLMMadaptive)
library(kableExtra)
library(knitr)
library(mvtnorm)
library(mgcv)
library(splines)
library(boot)
library(RColorBrewer)

# plot settings
theme_set(theme_minimal())
cols <- c(brewer.pal(3, "Set2"), "#000000") # define a color palette 
names(cols) <- c("6am", "12pm", "6pm",  "True")
# display.brewer.all()
# display.brewer.pal(4, "Set2")
# brewer.pal(5, "Set2")

## prediction window 
window <- seq(0, 1440, by = 360)
# M <- 5

# code 
source(here("Code/PredEval.R"))


#### load and format data ####
# load data
load(here("Data/DataAppl/ApplOutput_fGFPCA.RData"))
load(here("Data/DataAppl/ApplOutput_GLMMadaptive.RData"))
load(here("Data/DataAppl/Appl_GFOSR_L5.RData"))
load(here("Data/DataAppl/Appl_GFOSR_L1.RData"))

head(pred_appl_gfosr_l1)
head(pred_nhanes_fgfpca)

# format
## FGFPCA
pred_nhanes_fgfpca <- pred_nhanes_fgfpca %>% 
  select(!ends_with("b")) 

pred_nhanes_fgfpca$pred360[pred_nhanes_fgfpca$sind <= 360] <- NA
pred_nhanes_fgfpca$pred720[pred_nhanes_fgfpca$sind <= 720] <- NA
pred_nhanes_fgfpca$pred1080[pred_nhanes_fgfpca$sind <= 1080] <- NA

head(pred_nhanes_fgfpca)

# GFOSR
pred_appl_gfosr_l5 <- pred_appl_gfosr_l5 %>% 
  select(-window) %>% 
  rename(pred360=pred_w1, pred720=pred_w2, pred1080=pred_w3)

pred_appl_gfosr_l1 <- pred_appl_gfosr_l1 %>% 
  select(-window) %>% 
  rename(pred360=pred_w1, pred720=pred_w2, pred1080=pred_w3)

# GLMM
pred_nhanes_adglmm <-  pred_nhanes_adglmm %>% 
  select(-t)
head(pred_nhanes_adglmm)

#### Figure ####

# choose 4 subjects to plot
rand_id <- sample(unique(pred_nhanes_fgfpca$id), 4)

# individual prediction tracks
bind_rows(
  pred_nhanes_fgfpca %>% filter(id %in% rand_id) %>% 
    mutate(method="fGFPCA"),
  pred_nhanes_adglmm %>% filter(id %in% rand_id) %>%
    mutate(method = "GLMMadaptive"),
  pred_appl_gfosr_l1 %>% filter(id %in% rand_id) %>%
    mutate(method = "GFOSR (L=1)"),
  pred_appl_gfosr_l5 %>% filter(id %in% rand_id) %>%
    mutate(method = "GFOSR (L=5)") 
) %>% 
  mutate_at(vars(starts_with("pred")), function(x)exp(x)/(1+exp(x))) %>% 
  mutate(id = paste("ID","=", id)) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.1, alpha = 0.3)+
  geom_line(aes(x=sind, y = pred360, col = "6am"), linewidth=1.1)+
  geom_line(aes(x=sind, y = pred720, col = "12pm"), linewidth=1.1)+
  geom_line(aes(x=sind, y = pred1080, col = "6pm"), linewidth=1.1)+
  scale_x_continuous(name = "",
                     breaks = c(360, 720, 1080), 
                     labels = c("6am", "12pm", "6pm"))+
  facet_grid(cols = vars(method), rows = vars(id))+
  labs(col = "Maximum observation time", x = "Time", y="")+
  theme(legend.position = "bottom")+
  scale_color_manual(values = cols)
ggsave(filename = here("Images/DataAppl.pdf"), width=10, height=10)

#### AUC ####
## fGFPCA
pred_nhanes_fgfpca %>% 
  mutate(window = cut(sind, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>% 
  group_by(window) %>%
  summarise(auc1 = get_auc(Y, pred360),
            auc2 = get_auc(Y, pred720),
            auc3 = get_auc(Y, pred1080)) 

## GFOSR
pred_appl_gfosr_l1 %>% 
  mutate(window = cut(sind, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>% 
  group_by(window) %>%
  summarise(auc1 = get_auc(Y, pred360),
            auc2 = get_auc(Y, pred720),
            auc3 = get_auc(Y, pred1080)) 

pred_appl_gfosr_l5 %>% 
  mutate(window = cut(sind, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>% 
  group_by(window) %>%
  summarise(auc1 = get_auc(Y, pred360),
            auc2 = get_auc(Y, pred720),
            auc3 = get_auc(Y, pred1080)) 

## GLMM

pred_nhanes_adglmm %>% 
  mutate(window = cut(sind, breaks = window, include.lowest = T)) %>% 
  select(Y, starts_with("pred"), window) %>% 
  group_by(window) %>%
  summarise(auc1 = get_auc(Y, pred360),
            auc2 = get_auc(Y, pred720),
            auc3 = get_auc(Y, pred1080)) 
