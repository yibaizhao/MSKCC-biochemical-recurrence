library(here)
library(dplyr)
library(survival)
library(survminer)
library(gtsummary)
library(icenReg)

load(here("MSK/data", "rec_mskcc.RData"))
source(here("MSK/src", "analysis_functions.R"))

data_processing <- function(dset, digit){
  dset <- dset %>% mutate(
    time_l = ifelse(is.na(time_l), time_f, time_l),
    time_r = ifelse(is.na(time_r), Inf, time_r))
  # dset <- dset %>% mutate(
  #   time_l = ifelse(is.na(time_l), round(time_f, digit), floor(time_l)),
  #   time_r = ifelse(is.na(time_r), Inf, ceiling(time_r)))
  dset <- dset %>% mutate_at(vars(matches("dxage|dxpsa|txage|txpsa")), scale)
  # dset <- dset %>% mutate(status = ifelse(time_r==Inf, 0, 3))
  # categorize grade
  dset <- dset %>% mutate(txpathgrade_cat = ifelse(txpathgrade<=6, "<=6", ">6"))
  dset <- dset %>% filter(!(time_l==0 & time_r==Inf))
  dset <- dset %>% dplyr::select(-time_f)
  return(dset)
}
dset <- data_processing(dset = rec_mskcc, digit = 0)
dset_diff <- dset %>% filter(status==1) %>% mutate(diff=time_r-time_l)
hist(dset_diff$diff, freq = FALSE)

analysis <- function(dset){
  dset_sub <- dset %>% dplyr::select(c(ptid, time_l, time_r, txgrade, txpathgrade_cat, dxpsa, dxage)) %>% na.omit()
  npmle_fit <- ic_np(cbind(time_l, time_r) ~ txpathgrade_cat, data = dset_sub)
  lab <- names(npmle_fit$nGrp)
  npmle_data_0 <- data.frame(time = c(0, npmle_fit$scurves[[1]]$Tbull_ints[,1]),
                             surv = c(1, npmle_fit$scurves[[1]]$S_curves$baseline),
                             group = lab[1])
  npmle_data_1 <- data.frame(time = c(0, npmle_fit$scurves[[2]]$Tbull_ints[,1]),
                             surv = c(1, npmle_fit$scurves[[2]]$S_curves$baseline),
                             group = lab[2])
  npmle_data <- rbind(npmle_data_0, npmle_data_1) %>% filter(surv!=0)
  g <- ggplot(data = npmle_data, aes(x = time, y = surv, color = group)) +
    geom_step() +
    labs(x = "Time", y = "Survival Probability", color = "txpathgrade_cat") +
    theme_minimal()
  g
  # mixcure(Surv(time=time_l, time2=time_r, type = "interval2") ~ txgrade, 
  #         ~ txgrade,
  #         data = dset_sub)
  # fit_cure <- smcure(Surv(time=time_l, time2=time_r, type = "interval2")~txgrade+dxpsa+dxage,
  #        cureform=~dxpsa+dxage,
  #        data=dset_sub, model="ph", nboot = 100)
  # printsmcure(fit_cure)
  # predm <- predictsmcure(fit_cure,newX=c(0,1),newZ=c(0,1), model = "ph")
  # plotpredictsmcure(predm)
  # flexsurvcure(Surv(time=time_l, time2=time_r, type = "interval2")~txgrade+dxpsa+dxage,
  #              data = dset_sub, dist="weibullPH", link = "logistic", mixture = TRUE,
  #              anc = list(scale=~txgrade+dxpsa+dxage))
  # 
  # # The Generalized Odds Rate Mixture Cure model is fitted for interval censored survival data. 
  # library(GORCure)
  # fit<-GORMC(survfun=Surv(time=time_l, time2=time_r, type = "interval2")~txgrade+dxpsa+dxage,
  #            curefun=~dxpsa+dxage,
  #            data=dset_sub,r=0)
  # ggsurvplot(fit, data = dset_sub)
  
  
  survfun <- Surv(time_l, time_r, type = 'interval2')~txpathgrade_cat+txage
  curefun <- ~txpathgrade_cat+txage
  ests <- mxcureF(survfun, curefun, dset, bs_iter = NULL, dist = "weibull")
  point_ests <- ests$par
  var <- round(diag(solve(-ests$hessian)), 6)
  # ests <- mxcureF(survfun, curefun, dset, bs_iter = 100, dist = "npmle")
  # data <- mxcureF(survfun, curefun, dset)
  # L=data$sdata$L
  # R=data$sdata$R
  # delta=data$sdata$status
  # x=data$X
  # z=data$Z
  
  


  
}

