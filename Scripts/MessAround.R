head(fpca_fit$mu)

df_pred_latent %>% select(id, bin, eta_hat) %>% distinct(.) %>%
  group_by(bin) %>% summarise_at("eta_hat", mean)
