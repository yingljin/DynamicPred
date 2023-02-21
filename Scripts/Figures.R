# compare point prediction with in-sample

in_pred_subj1 %>% select(mid, type, eta, id, true)

p1 <- bind_rows(out_pred_subj1 %>% mutate(sample="Out-of-sample prediction"),
          in_pred_subj1 %>% select(mid, type, eta, id, true) %>%
            mutate(sample = "In-sample estimation")) %>%
  ggplot()+
  geom_line(aes(x=mid, y=true), col = "black")+
  geom_line(aes(x=mid, y=eta, col= type, group = type), 
            linetype="dashed", na.rm = T)+
  facet_wrap(~sample)+
  labs(x="", y = "", title = "Subject 1", col = "")

p2 <- bind_rows(out_pred_subj2 %>% mutate(sample="Out-of-sample prediction"),
                in_pred_subj2 %>% select(mid, type, eta, id, true) %>%
                  mutate(sample = "In-sample estimation")) %>%
  ggplot()+
  geom_line(aes(x=mid, y=true), col = "black")+
  geom_line(aes(x=mid, y=eta, col= type, group = type), 
            linetype="dashed", na.rm = T)+
  facet_wrap(~sample)+
  labs(x="", y = "", title = "Subject 2", col = "")

p3 <- bind_rows(out_pred_subj3 %>% mutate(sample="Out-of-sample prediction"),
                in_pred_subj3 %>% select(mid, type, eta, id, true) %>%
                  mutate(sample = "In-sample estimation")) %>%
  ggplot()+
  geom_line(aes(x=mid, y=true), col = "black")+
  geom_line(aes(x=mid, y=eta, col= type, group = type), 
            linetype="dashed", na.rm = T)+
  facet_wrap(~sample)+
  labs(x="", y = "", title = "Subject 3", col = "")

p4 <- bind_rows(out_pred_subj4 %>% mutate(sample="Out-of-sample prediction"),
                in_pred_subj4 %>% select(mid, type, eta, id, true) %>%
                  mutate(sample = "In-sample estimation")) %>%
  ggplot()+
  geom_line(aes(x=mid, y=true), col = "black")+
  geom_line(aes(x=mid, y=eta, col= type, group = type), 
            linetype="dashed", na.rm = T)+
  facet_wrap(~sample)+
  labs(x="", y = "", title = "Subject 4", col = "")

ggpubr::ggarrange(p1, p2, p3, p4, ncol = 1, common.legend = T) %>%
  annotate_figure(bottom = "Time", left = "Value")
ggsave(filename = "DynPred.jpeg", width=8, height=12)

# subject 1
p1 <- grid.arrange(
  out_pred_subj1 %>% 
    group_by(type) %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col = "black")+
    geom_line(aes(x=mid, y=eta, col= type, group = type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y = "", title = "Out-of-sample prediction")+
    theme(title = element_text(size = 10)),
  
  in_pred_subj1 %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col="black")+
    geom_line(aes(x=mid, y=eta, col=type, group=type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y="", title = "In-sample estimation")+
    theme(title = element_text(size = 10)), 
  nrow = 1,
  bottom = text_grob("Subject 1", vjust = -0.5, hjust = 0))

# subject 2
p2 <- grid.arrange(
  out_pred_subj2 %>% 
    group_by("type") %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col = "black")+
    geom_line(aes(x=mid, y=eta, col= type, group = type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y = "", title = "Out-of-sample prediction")+
    theme(title = element_text(size = 10)),
  
  in_pred_subj2 %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col="black")+
    geom_line(aes(x=mid, y=eta, col=type, group=type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y="", title = "In-sample estimation")+
    theme(title = element_text(size = 10)), 
  nrow = 1,
  bottom = text_grob("Subject 2", vjust = -0.5, hjust = 0))


# subject 3
p3 <- grid.arrange(
  out_pred_subj3 %>% 
    group_by("type") %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col = "black")+
    geom_line(aes(x=mid, y=eta, col= type, group = type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y = "", title = "Out-of-sample prediction")+
    theme(title = element_text(size = 10)),
  
  in_pred_subj3 %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col="black")+
    geom_line(aes(x=mid, y=eta, col=type, group=type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y="", title = "In-sample estimation")+
    theme(title = element_text(size = 10)), 
  nrow = 1,
  bottom = text_grob("Subject 3", vjust = -0.5, hjust = 0))


  # subject 4
p4 <- grid.arrange(
  out_pred_subj4 %>% 
    group_by("type") %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col = "black")+
    geom_line(aes(x=mid, y=eta, col= type, group = type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y = "", title = "Out-of-sample prediction")+
    theme(title = element_text(size = 10)),
  
  in_pred_subj4 %>%
    ggplot()+
    geom_line(aes(x=mid, y=true), col="black")+
    geom_line(aes(x=mid, y=eta, col=type, group=type), linetype="dashed", na.rm = T,
              show.legend = FALSE)+
    labs(x="", y="", title = "In-sample estimation")+
    theme(title = element_text(size = 10)), 
  nrow = 1,
  bottom = text_grob("Subject 4", vjust = -0.5, hjust = 0))

grid.arrange(p1, p2, p3, p4, ncol = 1, width=8, heights = 16)
