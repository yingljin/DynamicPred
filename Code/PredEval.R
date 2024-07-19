# This script codes the functions of evaluation of predictive performance

#### A function to calculate ISE ####

## calculation
# pred list: must contain columns: pred.., eta_i, id

calc_ISE <- function(pred_list, window){
  
  M <- length(pred_list)  # number of simulations
  
  # container
  ise_arr <- array(NA, dim = c(length(window)-1, length(window)-2, M))
  # prediction window by max obervation time by iteration
  
  # within simulation, average over subject
  for(m in 1:M){
    ise_tb <- pred_list[[m]] %>%
      mutate(err1 = (pred0.2-eta_i)^2,
             err2 = (pred0.4-eta_i)^2,
             err3 = (pred0.6-eta_i)^2,
             err4 = (pred0.8-eta_i)^2) %>%
      select(id, t, starts_with("err")) %>% 
      mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
      group_by(window, id) %>% 
      summarise_at(vars(err1, err2, err3, err4), sum) %>% 
      group_by(window) %>% 
      summarize_at(vars(err1, err2, err3, err4), mean) %>%
      # filter(window != "[0,0.2]") %>% 
      select(starts_with("err"))
    
    ise_arr[ , , m] <- as.matrix(ise_tb)
   
  }

  # average over simulation
  ise_arr <- apply(ise_arr, c(1, 2), mean, na.rm = T)
  colnames(ise_arr) <- c("0.2", "0.4", "0.6", "0.8")
  
  return(ise_arr)
}


#### Functions for AUC ####

## a function to calculate AUC with NA presence
get_auc <- function(y, pred){
  if(sum(is.na(y))>0 | sum(is.na(pred))>0){
    auc <- NA
  }
  else{
    this_perf <- performance(prediction(pred, y), measure = "auc")
    auc <- this_perf@y.values[[1]]
  }
  return(auc)
}


## calculate AUC for simulation

calc_AUC <- function(pred_list, window){
  
  M <- length(pred_list)  # number of simulations
  
  # container
  auc_arr <- array(NA, dim = c(length(window)-1, length(window)-2, M))
  # prediction window by max obervation time by iteration
  
  # calculate AUC over subjects within one simulation
  for(m in 1:M){
    auc_tb <- pred_list[[m]] %>% 
      mutate(window = cut(t, breaks = window, include.lowest = T)) %>% 
      select(Y, starts_with("pred"), window) %>%
      group_by(window) %>%
      summarise(auc1 = get_auc(Y, pred0.2),
                auc2 = get_auc(Y, pred0.4),
                auc3 = get_auc(Y, pred0.6),
                auc4 = get_auc(Y, pred0.8)) %>%
      select(starts_with("auc"))
    auc_arr[, ,m] <- as.matrix(auc_tb)
  }
  
  # average over simulation
  auc_arr <- apply(auc_arr, c(1, 2), mean, na.rm = T)
  colnames(auc_arr) <- c("0.2", "0.4", "0.6", "0.8")
  
  return(auc_arr)
}



#### Interval coverage rate ####
calc_predint_cover <- function(pred_list){
  
  covrate <- pred_list %>% bind_rows(.id="sim") %>%
    mutate(cover0.2 = pred0.2_lb<=eta_i & pred0.2_ub>=eta_i,
           cover0.4 = pred0.4_lb<=eta_i & pred0.4_ub>=eta_i,
           cover0.6 = pred0.6_lb<=eta_i & pred0.6_ub>=eta_i,
           cover0.8 = pred0.8_lb<=eta_i & pred0.8_ub>=eta_i) %>% 
    group_by(sim, t) %>%
    summarize_at(vars(starts_with("cover")), mean) %>% 
    group_by(t) %>%
    summarize_at(vars(starts_with("cover")), mean)
  
  return(covrate)
}

