
#### GLMMadaptive ####

pred_lst_ref <- list()
runtime_ref <- rep(NA, M)



pb = txtProgressBar(min = 0, max = M, initial = 0, style = 3) 

for(m in 1:M){
  df <- sim_data[[m]]
  
  t1 <- Sys.time()
  # GLMM adaptive model
  adglmm <- mixed_model(Y ~ sind, random = ~ sind | id, data =  df, family = binomial())
  
  # container
  df_est_latent <- df
  df_est_latent[, 'pred_t200'] <- df_est_latent[, 'pred_t400'] <- df_est_latent[, 'pred_t600'] <- df_est_latent[, 'pred_t800'] <- NA 
  
  # time up to t
  pred_t200 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 200),
                       newdata2 =  df %>% filter(sind_inx > 200), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t200[df_est_latent$sind_inx>200] <- pred_t200$newdata2$pred
  
  pred_t400 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 400),
                       newdata2 =  df %>% filter(sind_inx > 400), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t400[df_est_latent$sind_inx>400] <- pred_t400$newdata2$pred
  
  pred_t600 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 600),
                       newdata2 =  df %>% filter(sind_inx > 600), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t600[df_est_latent$sind_inx>600] <- pred_t600$newdata2$pred
  
  pred_t800 <- predict(adglmm, newdata = df %>% filter(sind_inx <= 800),
                       newdata2 =  df %>% filter(sind_inx > 800), type = "subject_specific", type_pred = "link",
                       se.fit = TRUE, return_newdata = TRUE)
  df_est_latent$pred_t800[df_est_latent$sind_inx>800] <- pred_t800$newdata2$pred
  
  t2 <- Sys.time()
  
  
  pred_lst_ref[[m]] <- df_est_latent
  runtime_ref[m] <- t2-t1
  
  setTxtProgressBar(pb, m)
  
}

save(pred_lst_ref, runtime_ref, file = here("Data/SimOutput_GLMMadaptive.RData"))

#### Plots #####

par(mfrow=c(2, 2))
plot(mid, fpca_mod$efunctions[,1])
plot(mid, fpca_mod$efunctions[,2])
plot(mid, fpca_mod$efunctions[,3])
plot(mid, fpca_mod$efunctions[,4])



pred_lst_ref[[1]] %>% filter(id==5) %>%
  left_join(df %>% filter(id==5) %>%  select(sind, eta_i)) %>%
  ggplot()+
  geom_line(aes(x=sind, y=eta_i))+
  geom_line(aes(x=sind, y = pred_t200, col = "200"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t400, col = "400"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t600, col = "600"), linetype = "dashed")+
  geom_line(aes(x=sind, y = pred_t800, col = "800"), linetype = "dashed")


mean(runtime_ref)
