
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
cols <- c(brewer.pal(4, "Set2"), "#000000") # define a color palette 
names(cols) <- c("0.2", "0.4", "0.6", "0.8", "True")
# display.brewer.all()
# display.brewer.pal(4, "Set2")
# brewer.pal(5, "Set2")

## prediction window 
window <- seq(0, 1, by = 0.2)
# M <- 5

# code 
source(here("Code/PredEval.R"))



#### Simulation output ####
# load data
## N = 500
load(here("Data/SimN500/SimOutput_fGFPCA.RData"))
load(here("Data/SimN500/SimOutput_GLMMadaptive.RData"))
load(here("Data/SimN500/SimOutput_GFOSR_L1.RData"))
load(here("Data/SimN500/SimOutput_GFOSR_L5.RData"))

## N = 100
load(here("Data/SimN100/SubSimOutput_fGFPCA.RData"))
load(here("Data/SimN100/SubSimOutput_GLMMadaptive.RData"))
load(here("Data/SimN100/SubSimOutput_GFOSR_L1.RData"))
load(here("Data/SimN100/SubSimOutput_GFOSR_L5.RData"))

# choose 4 subjects to plot
rand_id <- sample(501:600, size = 4)
rand_id <- sample(101:200, size = 4)

#### Format ####
##### N = 500 #####
pred_list_fGFPCA <- lapply(pred_list_fGFPCA, function(x){x %>% select(-bin)})
name_vec <- colnames(pred_list_fGFPCA[[1]])
pred_list_fGFPCA <- lapply(pred_list_fGFPCA,
                           function(x){
                             x[x$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
                             x[x$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
                             x[x$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
                             x[x$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA
                             return(x)
                           })


pred_list_GLMMad <- lapply(pred_list_GLMMad, function(x){
  colnames(x) <- name_vec
  return(x)
})

pred_list_gfofr_l1 <- lapply(pred_list_gfofr_l1, function(x){
  x <- x %>% select(-window)
  colnames(x) <- name_vec
  return(x)
})

pred_list_gfofr_l5 <- lapply(pred_list_gfofr_l5, function(x){
  x <- x %>% select(-window)
  colnames(x) <- name_vec
  return(x)
})


##### N = 100 #####
# N = 100
name_vec <- colnames(pred_subset_fGFPCA[[1]])
pred_list_fGFPCA <- lapply(pred_subset_fGFPCA, 
                           function(x){
                             x[x$t<=0.2, c("pred0.2", "pred0.2_lb", "pred0.2_ub")] <- NA
                             x[x$t<=0.4, c("pred0.4", "pred0.4_lb", "pred0.4_ub")] <- NA
                             x[x$t<=0.6, c("pred0.6", "pred0.6_lb", "pred0.6_ub")] <- NA
                             x[x$t<=0.8, c("pred0.8", "pred0.8_lb", "pred0.8_ub")] <- NA
                             return(x)    
                           })
head(pred_list_fGFPCA[[1]])


head(pred_subset_adglmm[[1]])
which(num_probs==1)
pred_subset_adglmm <- pred_subset_adglmm[!num_probs] # remove data with numeric problems
length(pred_subset_adglmm) # 6 nonconvergence problem
pred_list_GLMMad <- lapply(pred_subset_adglmm, function(x){
  colnames(x) <- name_vec
  return(x)
})

head(pred_subset_gfofr_l1[[1]])
pred_list_gfofr_l1 <- lapply(pred_subset_gfofr_l1, function(x){
  x <- x %>% select(-window)
  colnames(x) <- name_vec
  return(x)
})

pred_list_gfofr_l5 <- lapply(pred_subset_gfofr_l5, function(x){
  x <- x %>% select(-window)
  colnames(x) <- name_vec
  return(x)
})


# individual prediction tracks
bind_rows(
  pred_list_fGFPCA[[1]] %>% filter(id %in% rand_id) %>% 
    mutate(method="fGFPCA"),
  pred_list_GLMMad[[1]] %>% filter(id %in% rand_id) %>%
    mutate(method = "GLMMadaptive"),
  pred_list_gfofr_l1[[1]] %>% filter(id %in% rand_id) %>%
    mutate(method = "GFOSR (L=1)"),
  pred_list_gfofr_l5[[1]] %>% filter(id %in% rand_id) %>%
    mutate(method = "GFOSR (L=5)") 
) %>% 
  mutate_at(vars(eta_i, starts_with("pred")), 
            .funs = function(x){exp(x)/(1+exp(x))}) %>%
  mutate(method=factor(method, 
                       levels = c("fGFPCA", "GFOSR (L=5)", "GFOSR (L=1)","GLMMadaptive"))) %>%
  mutate(id = factor(id, levels = unique(id), labels = paste0("Subject ", 1:4))) %>%
  ggplot()+
  geom_point(aes(x=t, y=Y), size = 0.1, alpha = 0.3)+
  geom_line(aes(x=t, y=eta_i, col = "True"), linetype = "dashed", linewidth=1.1)+
  geom_line(aes(x=t, y=pred0.2, col = "0.2"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=pred0.4, col = "0.4"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=pred0.6, col = "0.6"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=pred0.8, col = "0.8"), na.rm = T, linewidth=1.1)+
  guides(alpha = "none")+
  facet_grid(rows = vars(id), cols = vars(method))+
  labs(col = "Maximum observation time", x = "t", y="")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = cols)
ggsave(filename = here("Images/simN500.pdf"), width=10, height=10)
ggsave(filename = here("Images/simN100.pdf"), width=10, height=10)


# coverage rate of prediction interval
bind_rows(calc_predint_cover(pred_list_fGFPCA) %>% mutate(method = "fGFPCA"),
          calc_predint_cover(pred_list_gfofr_l5) %>% mutate(method = "GFOSR (L=5)"),
          calc_predint_cover(pred_list_gfofr_l1) %>% mutate(method = "GFOSR (L=1)"),
          calc_predint_cover(pred_list_GLMMad) %>% mutate(method = "GLMMadaptive")) %>%
  ggplot()+
  geom_line(aes(x=t, y=cover0.2, col = "0.2"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=cover0.4, col = "0.4"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=cover0.6, col = "0.6"), na.rm = T, linewidth=1.1)+
  geom_line(aes(x=t, y=cover0.8, col = "0.8"), na.rm = T, linewidth=1.1)+
  facet_grid(cols = vars(method))+
  labs(col = "Maximum observation time", x = "t", y="Coverage rate")+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = cols)
ggsave(filename = here("Images/simN500_ci.pdf"), width=10, height=3)
ggsave(filename = here("Images/simN100_ci.pdf"), width=10, height=3)


#### ISE and AUC ####
calc_ISE(pred_list_fGFPCA, window)


calc_AUC(pred_list_fGFPCA, window)

mean(time_subset_fGFPCA)
mean(time_gfofr_l1)
mean(time_gfofr_l5)
mean(time_GLMMad)


