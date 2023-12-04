


#### Example figures for presentation ####

library(here)
library(tidyverse)
library(mgcv)
theme_set(theme_minimal())

# binary function

df <- readRDS(here("Data/nhanes_bi.rds"))
length(unique(df$SEQN))
length(unique(df$sind))

## all subjects

df1 <- df %>% filter(SEQN == 62161)
df2 <- df %>% filter(SEQN == 62209)

df1$eta_hat <- predict(loess(Z~sind, df1))
df2$eta_hat <- predict(loess(Z~sind, df2))

df1$eta_hat[df1$eta_hat<0] <- 0
df2$eta_hat[df2$eta_hat<0] <- 0

bind_rows(df1, df2) %>%
  # mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_point(aes(x=sind, y=Z), size = 0.5)+
  geom_line(aes(x=sind, y=eta_hat), col = "blue")+
  facet_wrap(~gender, ncol = 2, nrow = 1)+
  labs(title = "Binary activity indicator from NHANES data", y = "", x = "Minute")
ggsave(filename = "bifunc.png", path = "Images", width = 12, height = 6, bg = "white")


# try glm
library(lme4)

fit_glm <- glmer(Z ~ sind + (1 | sind), data = df, family = binomial)

# continous function

list.files("Data")
df2 <- readRDS(here("Data/NHANES_AC_processed.rds"))

df2 %>% 
  filter(SEQN==21005) %>% 
  select(WEEKDAY, starts_with("min")) %>%
  group_by(WEEKDAY) %>%
  # summarise_at(-1, mean) %>%
  pivot_longer(-1, names_to = "Minute") %>% 
  mutate(Minute = gsub("MIN", "", Minute)) %>%
  mutate(Minute = as.numeric(Minute),
         WEEKDAY = as.factor(WEEKDAY)) %>%
  ggplot(aes(x=Minute, y=value, col=WEEKDAY, group = WEEKDAY))+
  geom_point(size = 0.5, alpha=0.5)+
  geom_smooth(se=F, size = 0.5)+
  labs(y="", col = "Weekday", title = "Total activitiy count from one subject")
ggsave(filename = "contfunc_multilevel.png", path = "Images", width = 6, height = 5, bg = "white")

# count function
library(refund)
data(cd4)

cd4[1:5, ] %>% data.frame() %>%
  pivot_longer(1:61, names_to = "month", values_to = "count") %>%
  mutate(month = rep(-18:42, 5), 
         id = rep(1:5, each=61)) %>%
  filter(complete.cases(.)) %>%
  ggplot(aes(x=month, y=count, col=as.factor(id), group = id))+
  geom_point()+
  geom_line()+
  labs(x="Month", y="CD4 cell count", col = "ID", title = "Observed CD4 cell counts")
ggsave(filename = "countfunc.png", path = "Images", width = 6, height = 6, bg = "white") 

# image
data(DTI) # ID = 2001

index <- which(DTI$ID==2001)
DTI$cca[index, ] %>% 
  data.frame() %>%
  # mutate(gender = c(rep("Female", 5), rep("Male", 6))) %>%
  # group_by(gender) %>%
  # summarize_at(1:93, mean) %>%
  pivot_longer(1:93, names_to = "index", values_to = "val") %>% 
  mutate(index = rep(1:93, 5), visit = rep(1:5, each = 93)) %>% 
  mutate(visit = as.factor(visit)) %>%
  ggplot(aes(x=index, y=val, col=visit, group= visit))+
  geom_point(size = 0.5)+
  geom_line()+
  labs(title = "Fractional anisotropy tract profiles for the corpus callosum",
       x = "Position", y = "Intensity", col = "Visit")
ggsave(filename = "DTI_multilevel.png", path = "Images", width = 6, height = 5, bg = "white") 
  
#### Local GLMMs ####
load(here("Data/SimOutput_fGFPCA.RData"))
load(here("Data/Sim_data.RData"))
library(lme4)

head(sim_data[[1]])

id_vec <- c(223, 292, 322, 375)

df <- sim_data[[1]] %>% filter(id %in% id_vec) %>% filter(sind_inx <= 50) #%>%

# observation  
df %>% 
  ggplot(aes(x=sind_inx, y=Y, col = id))+
  geom_point()+
  labs(x="Time index", y="Outcome", col = "ID")
ggsave(filename = "exp_Y.png", path = "Images", width = 6, height = 6, bg = "white") 

# bins
df %>% 
  ggplot(aes(x=sind_inx, y=Y, col = id))+
  geom_point()+
  geom_vline(xintercept = c(10, 20, 30, 40, 50), linetype="dashed")+
  labs(x="Time index", y="Outcome", col = "ID")
ggsave(filename = "exp_bin.png", path = "Images", width = 6, height = 6, bg = "white") 

# estimate latent function

source(here("Code/GLMM-FPCA.R"))


df$bin <- cut(df$sind_inx, breaks = c(0, 10, 20, 30, 40, 50), include.lowest = T, 
              labels = c(5, 15, 25, 35, 45))
df$bin <- as.numeric(as.character(df$bin))

df_bin_lst <- split(df, f = df$bin)

df_est_latent <- lapply(df_bin_lst, function(x){pred_latent(x)}) 
df_est_latent <- bind_rows(df_est_latent)

df_est_latent %>% 
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_point(aes(x=sind_inx, y=Y, col = id))+
  geom_vline(xintercept = c(10, 20, 30, 40, 50), linetype="dashed")+
  geom_point(aes(x=bin, y=eta_hat, col = id), shape = 5)+
  labs(x="Time index", y="Outcome", col = "ID")
ggsave(filename = "exp_est_eta.png", path = "Images", width = 6, height = 6, bg = "white") 


df_est_latent %>% 
  mutate(eta_hat = exp(eta_hat)/(1+exp(eta_hat))) %>%
  ggplot()+
  geom_point(aes(x=sind_inx, y=Y, col = id))+
  geom_vline(xintercept = c(10, 20, 30, 40, 50), linetype="dashed")+
  geom_point(aes(x=bin, y=eta_hat, col = id), shape = 5)+
  geom_line(aes(x=bin, y=eta_hat, col=id))+
  labs(x="Time index", y="Outcome", col = "ID")
ggsave(filename = "exp_est_eta_line.png", path = "Images", width = 6, height = 6, bg = "white") 

#### fpca #####

df <- data.frame(
  time = 1:ncol(fpca_mod$Y), 
  Y = fpca_mod$Y[223, ],
  Mean = fpca_mod$mu, 
  PC1 = fpca_mod$efunctions[, 1],
  PC2 = fpca_mod$efunctions[, 2],
  PC3 = fpca_mod$efunctions[, 3],
  PC4 = fpca_mod$efunctions[, 4]
)

df %>% pivot_longer(2:7, names_to = "Function") %>%
  mutate(Function = relevel(as.factor(Function), ref = "Y")) %>%
  ggplot()+
  geom_line(aes(x=time, y=value))+
  facet_wrap(~Function, scales = "free", ncol = 3, dir="v")+
  labs(x="Time index", y="")
ggsave(filename = "fpca.png", path = "Images", width = 9, height = 6, bg = "white") 
