############################
# Use bootstrap to calculate variance of cumulative incident function
############################
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(splines)
library(nnet)
library(survival)
library(purrrlyr)
library(foreach)
library(doParallel)

datestamp <- "2022-11-29"

filename <- "data_2022-09-13.csv"

read_data <- function(filename){
  dset <- read.csv(here("data", filename))
  # 0. Filter out "Chemo" and patients having adjuvant treatment
  dset <- dset %>% filter(salvage_tx_type != "Chemo" | is.na(salvage_tx_type))
  dset <- dset %>% filter(is.na(tt_adjuvant_tx))
  # 1. Categorization
  dset <- dset %>%
    mutate(tx_grade_gp = case_when(
      tx_grade == "6" ~ "<=6",
      tx_grade %in% c("8", "9", "10") ~ ">7",
      TRUE ~ tx_grade
    ),
    tx_grade_gp = factor(tx_grade_gp, levels = c("<=6", "3+4", "4+3", ">7")))
  
  dset <- dset %>% 
    mutate(tx_tstage_gp = case_when(
      grepl("T2|T0", tx_tstage) ~ "<=T2",
      grepl("T3", tx_tstage) ~ "T3",
      grepl("T4", tx_tstage) ~ "T4"
    ))
  
  dset <- dset %>%
    mutate(race_gp = case_when(
      grepl("Asian|Other", race) ~ "Other",
      TRUE ~ race),
      # 2. Race="White" as baseline
      race_gp = factor(race_gp, levels = c("White", "Black", "Other")))
  
  # 3. Reorer therapy factor
  dset <- dset %>%
    mutate(across(salvage_tx_type, ~ifelse(is.na(tt_calendar_salvage_tx), "None", .)),
           salvage_tx_type = factor(salvage_tx_type, levels = c("None", "Hormones", "RT", "Hormones + RT")))
  # age group
  dset <- dset %>%
    mutate(tx_age_gp = cut(tx_age, c(0, 55, 70,100)))
  # Categorize calendar year
  dset <- dset %>%
    mutate(
      across(matches("calendar|year"), ~ifelse(.>= 2001 & .<=2003, "2001-2003",
                                               ifelse(.>= 2004 & .<=2006, "2004-2006",
                                                      ifelse(.>= 2007 & .<=2009, "2007-2009",
                                                             ifelse(.>= 2010 & .<=2012, "2010-2012",
                                                                    ifelse(.>= 2013 & .<=2015, "2013-2015",
                                                                           ifelse(.>= 2016 & .<=2018, "2016-2018",
                                                                                  ifelse(.>= 2019 & .<=2021, "2019-2021", NA))))))),
             .names = "{col}_gp"),
      tt_calendar_bcr_gp2 = ifelse(tt_calendar_bcr >= 2001 & tt_calendar_bcr <= 2005, "2001-2005",
                                   ifelse(tt_calendar_bcr > 2005 & tt_calendar_bcr <= 2010, "2006-2010",
                                          ifelse(tt_calendar_bcr > 2010 & tt_calendar_bcr <= 2015, "2011-2015",
                                                 ifelse(tt_calendar_bcr > 2015 & tt_calendar_bcr <= 2021, "2016-2021", NA)))),
      across(matches("calendar|year"), ~factor(.)))
  # if no salvage treatment, then tt_salvage_tx=tt_last_fu
  dset <- dset %>% mutate(across(tt_salvage_tx, ~case_when(salvage_tx_type == "None" ~ tt_last_fu,
                                                           TRUE ~ .)))
  # time(month) from bcr to salvage treatment
  dset <- dset %>% mutate(tt_bcr_to_salvage_tx = tt_salvage_tx-tt_bcr)
  
  # categorize comobility
  dset <- dset %>%
    mutate(cci_gp = case_when(cci == 0 ~ "0",
                              cci <= 2 ~ "1-2",
                              TRUE ~ "3+")
    )
  
  
  return(dset)
}

# CIF
integrand <- function(t, df_new, fun_pi, fun_ft, trt_type){
  df_new <- data.frame(tx_grade_group=df_new$tx_grade_group,
                       tt_bcr_to_salvage_tx_year=t,
                       tx_age=df_new$tx_age,
                       tt_calendar_bcr_gp2=df_new$tt_calendar_bcr_gp2,
                       tt_bcr=df_new$tt_bcr,
                       tx_tstage_gp=df_new$tx_tstage_gp)
  pi <- predict(fun_pi, newdata=df_new,type="probs")
  ft_risk <- predict(fun_ft, newdata=df_new, type='link')
  shape <- 1/fun_ft$scale
  scale <- as.numeric(exp(ft_risk))
  f_t <- dweibull(x=t, shape=shape, scale=scale, log = FALSE)
  (pi*f_t)[,trt_type]
}

fit_vertical_model <- function(dset, last_follow_up_years = seq(0.01, 15, 0.1)){
  dset <- dset %>% mutate(tt_bcr_to_salvage_tx_year = tt_bcr_to_salvage_tx/12,
                          tt_bcr_to_salvage_tx_year = ifelse(tt_bcr_to_salvage_tx_year==0, 0.001, tt_bcr_to_salvage_tx_year))
  dset_trt <- dset %>% dplyr::filter(salvage_tx_type!="None")
  dset_trt <- dset_trt %>% mutate(salvage_tx_type = factor(salvage_tx_type))
  
  # estimate conditional probability of receiving each treatment at time t
  mlog_spline <- multinom(salvage_tx_type ~ bs(tt_bcr_to_salvage_tx_year) + tt_calendar_bcr_gp2+tx_age+tx_grade_group+tt_bcr+tx_tstage_gp, data = dset_trt)
  dset_cov <- dset %>% select(tt_calendar_bcr_gp2, tx_age, tx_grade_group, tt_bcr, tx_tstage_gp)
  newdset <- data.frame(tx_age=63,
                        tt_bcr=13,
                        tx_tstage_gp="T3")
  newdset <- expand_grid(newdset, 
                         tt_bcr_to_salvage_tx_year=last_follow_up_years,
                         tx_grade_group=c("<=6", "3+4", "4+3", ">=8"),
                         tt_calendar_bcr_gp2 = head(levels(dset_cov$tt_calendar_bcr_gp2), 3))
  # vary years based on maximum followups
  newdset <- newdset %>% filter((tt_calendar_bcr_gp2=="2001-2005"&tt_bcr_to_salvage_tx_year<=20)|
                                  (tt_calendar_bcr_gp2=="2006-2010"&tt_bcr_to_salvage_tx_year<=15)|
                                  (tt_calendar_bcr_gp2=="2011-2015"&tt_bcr_to_salvage_tx_year<=10)|
                                  (tt_calendar_bcr_gp2=="2016-2021"&tt_bcr_to_salvage_tx_year<=5))
  probs <- predict(mlog_spline, newdata=newdset, type="probs")
  colnames(probs) <- levels(dset_trt$salvage_tx_type)
  gdset <- data.frame(newdset, probs)
  
  # Assume T~Weibull, estimate density corresponding to the overall failure time distribution
  fit <- survreg(Surv(tt_bcr_to_salvage_tx_year,salvage_tx_type!="None") ~ tt_calendar_bcr_gp2+tx_age+tx_grade_group+tt_bcr+tx_tstage_gp, data=dset)
  # add salvage type
  newdset <- expand_grid(newdset,
                         salvage_tx_type = c("Hormones", "RT", "Hormones + RT"))
  out <- newdset %>% by_row(function(row){
    integrate(integrand, 
              df_new=row, fun_pi=mlog_spline, fun_ft=fit, 
              trt_type = row$salvage_tx_type, 
              lower = 0, upper = row$tt_bcr_to_salvage_tx_year)$value
  }, .collate = "rows", .to = "CIF")
}

bootstrap_analysis <- function(dset, B, numCores, saveit=TRUE){
  registerDoParallel(numCores)  # use multicore, set to the number of our cores
  bout <- foreach (i=1:B, .combine=rbind) %dopar% {
    bset <- dset %>% sample_frac(frac=1, replace = T)
    bout_i <- fit_vertical_model(bset, last_follow_up_years = seq(0.01, 15, 0.1))
    bout_i$B <- i
    bout_i
  }
  if(saveit){
    filename <- str_glue('bs_out_{datestamp}.csv')
    write.csv(bout, here("results", filename))
  }
}

control <- function(filename){
  dset <- read_data(filename)
  bootstrap_analysis(dset, B=20, numCores=32, saveit=TRUE)
}
control(filename = filename)

