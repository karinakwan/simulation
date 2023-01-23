remotes::install_github("karinakwan/casebaseweights@95375b14b1b8876dc01db6ade5bcdc35af714657")
library(survival)
library(parallel)

create_data <- function(n, x, beta, event.scale, event.shape,
                        cens.scale, cens.shape){
  
  true_time <- as.numeric((-log(runif(n))/(exp(x*beta)*event.scale))^(1/event.shape))
  cens_time <- rweibull(n, shape = cens.shape, scale = cens.scale^(-1/cens.shape))
  
  # Censoring 
  obs_time <- pmin(true_time, cens_time)
  event <- ifelse(obs_time < cens_time, 1, 0)
  
  if (n == 400){
    df <- data.frame(time = obs_time, x = x, status = event)
    df$wts <- rep(1, length(x))
    
    return(df)
    
  } else {
    
    # get equal number of cases and controls by randomly sampling controls
    # from n-(number of cases)
    n.ctrls <- sum(event)
    cases.index <- which(event == 1)
    ctrls.index <- sample(which(event == 0), size = n.ctrls, replace = FALSE)
    df.sample <- data.frame(time = obs_time[c(cases.index, ctrls.index)], 
                            x = x[c(cases.index, ctrls.index)], 
                            status = event[c(cases.index, ctrls.index)])
    
    df.sample$wts <- ifelse(df.sample$status==0, 
                            (n-sum(df.sample$status==1))/sum(df.sample$status==0), 
                            1)
    
    df.full <- data.frame(time = obs_time, x = x, status = event)
    
    return(list(df.full, df.sample))
    
  }
  
}

fit_models <- function(full_dataset, sample_dataset, ratio){
  
  sample.wts <- sample_dataset$wts
  
  # Cox full
  temp_cox_full <- coxph(Surv(time, status) ~ x, data = full_dataset)
  cox_out_full <- c(cox_coef_full = coef(summary(temp_cox_full))[1],
                    cox_se_full = coef(summary(temp_cox_full))[3])
  
  # Cox naive
  temp_cox_naive <- coxph(Surv(time, status) ~ x, data = sample_dataset)
  cox_out_naive <- c(cox_coef_naive = coef(summary(temp_cox_naive))[1],
                     cox_se_naive = coef(summary(temp_cox_naive))[3])
  
  # Cox robust 
  temp_cox_robust <- coxph(Surv(time, status) ~ x, data = sample_dataset, 
                           weights = sample.wts,
                           robust = TRUE)
  cox_out_robust <- c(cox_coef_robust = coef(summary(temp_cox_robust))[1],
                      cox_se_robust = coef(summary(temp_cox_robust))[4])
  
  # casebase full
  cb_out_full <- lapply(ratio, FUN = function(ratio_i){
    temp_cb_full <- casebase::fitSmoothHazard(status ~ log(time) + x, 
                                              data = full_dataset, 
                                              ratio = ratio_i)
    cb_out <- c(coef(summary(temp_cb_full))[3, 1],
                coef(summary(temp_cb_full))[3, 2])
    names(cb_out) <- paste0("cb_r_", ratio_i,  c("_coef", "_se"), "_full")
    return(cb_out)
  })
  
  # casebase naive 
  cb_out_naive <- lapply(ratio, FUN = function(ratio_i){
    temp_cb_naive <- casebase::fitSmoothHazard(status ~ log(time) + x, 
                                               data = sample_dataset, 
                                               ratio = ratio_i)
    cb_out <- c(coef(summary(temp_cb_naive))[3, 1],
                coef(summary(temp_cb_naive))[3, 2])
    names(cb_out) <- paste0("cb_r", ratio_i,  c("_coef", "_se"), "_naive")
    return(cb_out)
  })
  
  sample_dataset$wts <- sample.wts
  wts <- sample.wts
  
  # casebase robust 
  cb_out_robust <- lapply(ratio, FUN = function(ratio_i){
    temp_cb_robust <- casebaseweights::fitSmoothHazard(status ~ log(time) + x,
                                                       data = sample_dataset, 
                                                       ratio = ratio_i, wts = wts)
    cb_out <- c(coef(summary(temp_cb_robust))[3, 1],
                coef(summary(temp_cb_robust))[3, 2])
    names(cb_out) <- paste0("cb_r", ratio_i,  c("_coef", "_se"), "_robust")
    return(cb_out)
  })
  
  return(c(cox_out_full, cox_out_naive, cox_out_robust, unlist(cb_out_full),
           unlist(cb_out_naive), unlist(cb_out_robust)))
  
}


sim_run <- function(n, x, beta, event.scale, event.shape,
                    cens.scale, cens.shape, iternb, path, ratio){
  
  df <- create_data(n = n, x = x,
                    beta = beta, event.scale = event.scale, 
                    event.shape = event.shape, cens.scale = cens.scale,
                    cens.shape = cens.shape)
  
  modelfit <- fit_models(full_dataset = df[[1]], sample_dataset = df[[2]], 
                         ratio = ratio)
  filename <- paste0(path, "modelout_", iternb, ".rds")
  saveRDS(modelfit, filename)
  
}

# n=10000
set.seed(251)
n <- 10000
x.mean <- 0
x.sd <- 1
x <- rlnorm(n, meanlog = x.mean, sdlog = x.sd)
beta <- 1.5
event.scale <- 10^(-8)
event.shape <- 4
cens.scale <- 10^(-5)
cens.shape <- 8

sim_10000 <- mclapply(1:1000, sim_run, n = n, x = x,
                      beta = beta, event.scale = event.scale, 
                      event.shape = event.shape, cens.scale = cens.scale,
                      cens.shape = cens.shape, 
                      path = "path_to_save_model_output/",
                      ratio = c(100,200,500), mc.cores = 16)
