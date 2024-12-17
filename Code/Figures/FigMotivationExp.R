
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
load(here("Data//DataAppl/ApplOutput_fGFPCA.RData"))
length(unique(df$SEQN)) # 8763 subjects
df %>% select(SEQN, age_years_interview, BMI, gender) %>% distinct(.) %>% 
  select(BMI) %>% 
  lapply(summary)



# age: 49.53 (17.52)
# BMI: 29.08 (7.04)
# gender: 4202 male, 4561 female

4561/8763

#### Figure ####

# motivation example a
## data structure
id_vec <- unique(pred_nhanes_fgfpca$id)[1:6]

p1 <- pred_nhanes_fgfpca %>% 
  filter(id %in% id_vec) %>%
  select(id, Y, sind, starts_with("pred1080")) %>% 
  mutate_at(vars(starts_with("pred1080")), function(x)exp(x)/(1+exp(x))) %>%
  mutate(id = paste("ID", "=", id )) %>% 
  ggplot()+
  geom_line(aes(x=sind, y = pred1080), col = "red", size = 1.1)+
  geom_point(aes(x=sind, y = Y), size = 0.1, alpha = 0.3)+
  facet_wrap(~id)+
  scale_x_continuous(name = "",
                     breaks = c(360, 720, 1080), 
                     labels = c("6am", "12pm", "6pm"))+
  scale_y_continuous(name = "Observed indicator",
                     breaks = c(0, 1),
                     sec.axis = dup_axis(name="Predicted probablity",
                                         breaks = seq(0, 1, by = 0.25)))+
  theme(legend = element_blank(),
        axis.title.y.right = element_text(color="red"),
        axis.text.y.right = element_text(color="red"))
ggsave(here("Images/MotiveExp_a.pdf"),
       width=10, height=5)



# motivation example b
# dynamic prediction illustration
## one male (62161) and one female (62164)
pred_nhanes_fgfpca %>% select(id, gender) %>% distinct(.) %>% head()

df_plot <- pred_nhanes_fgfpca %>% filter(id==62164 | id == 62184) %>% 
  select(id, Y, sind, pred360, pred720, pred1080) %>% 
  pivot_longer(starts_with("pred")) 

df_plot$Y[df_plot$name=="pred360" & df_plot$sind>360] <- NA
df_plot$Y[df_plot$name=="pred720" & df_plot$sind>720] <- NA
df_plot$Y[df_plot$name=="pred1080" & df_plot$sind>1080] <- NA

df_plot$value[df_plot$name=="pred360" & df_plot$sind <= 360] <- NA
df_plot$value[df_plot$name=="pred720" & df_plot$sind <= 720] <- NA
df_plot$value[df_plot$name=="pred1080" & df_plot$sind <= 1080] <- NA



p2 <- df_plot %>% 
  mutate_at("value", function(x)exp(x)/(1+exp(x))) %>% 
  mutate(name = as.numeric(gsub("pred", "", name))) %>%
  mutate(max_t=factor(name, levels = c(360, 720, 1080), 
                      labels = c("Observed up to 6am", 
                                 "Observed up to 12pm", 
                                 "Observed up to 6pm"))) %>% 
  mutate(id = paste("ID", "=", id)) %>% 
  ggplot()+
  geom_point(aes(x=sind, y=Y), size = 0.1, alpha = 0.3)+
  geom_line(aes(x=sind, y = value), col = "red", linewidth=1.1)+
  geom_vline(aes(xintercept = name), linetype = "dashed")+
  facet_grid(cols = vars(max_t), rows = vars(id))+
  scale_x_continuous(name = "",
                     breaks = c(360, 720, 1080), 
                     labels = c("6am", "12pm", "6pm"))+
  scale_y_continuous(name = "Observed indicator",
                     breaks = c(0, 1),
                     sec.axis = dup_axis(name="Predicted probablity",
                                         breaks = seq(0, 1, by = 0.25)))+
  theme(legend = element_blank(),
        axis.title.y.right = element_text(color="red"),
        axis.text.y.right = element_text(color="red"))
ggsave(here("Images/MotiveExp_b.pdf"),
       width=10, height=5)


ggarrange(p1, p2, col=1, 
          )
