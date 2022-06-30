requiredPackages = c('foreach', 'parallel', 'doParallel')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

piF <- function(intercept, slope, x){
  coef <- matrix(c(intercept, slope), ncol = 1)
  linear <- x%*%coef
  exp(linear)/(1+exp(linear))
}

logSurvF <- function(t, tau, surv_param, beta, z){
  lam <- c(unlist(sapply(t, function(i) ifelse(is.infinite(i), i, surv_param[which(i==tau)]))))
  beta <- matrix(beta, ncol = 1)
  -lam*c(exp(z%*%beta))
}

SurvF <- function(t, tau, surv_param, beta, z, dist){
  if(dist == "npmle"){
    log_surv <- logSurvF(t, tau, surv_param, beta, z)
    surv <- exp(log_surv)
  }
  if(dist == "weibull"){
    u <- surv_param[1]
    v <- surv_param[2]
    surv <- exp(-exp(u + exp(v)*log(t) + z%*%beta))
  }
  
  return(surv)
}

loglikF <- function(params, L, R, x, z, delta, dist){
  tau <- sort(unique(c(L, R)))
  tau <- tau[!is.infinite(tau)]
  n_cure <- ncol(x)
  n_beta <- ncol(z)
  
  intercept <- params[1]
  slope <- params[2:n_cure]
  beta <- params[(n_cure+1):(n_cure+n_beta)]
  if(dist == "npmle"){
    surv_param <- tail(params, length(tau))
  }
  if(dist == "weibull"){
    surv_param <- tail(params, 2) # c(u,v): u=log(labmda); v=log(gamma)
  }
  
  pi <- piF(intercept, slope, x)
  Surv_L <- SurvF(L, tau, surv_param, beta, z, dist)
  Surv_R <- SurvF(R, tau, surv_param, beta, z, dist)
  
  likli <- (Surv_L-Surv_R)^delta * (1-pi + pi*Surv_L)^(1-delta)
  loglikli <- log(likli)
  loglikli[is.infinite(loglikli)] <- 0
  sum(loglikli)
}

EiF <- function(ests, L, R, x, z, delta, dist){
  tau <- sort(unique(c(L, R)))
  tau <- tau[!is.infinite(tau)]
  n_cure <- ncol(x)
  n_beta <- ncol(z)
  
  intercept <- ests[1]
  slope <- ests[2:n_cure]
  beta <- ests[(n_cure+1):(n_cure+n_beta)]
  if(dist == "npmle"){
    surv_param <- tail(ests, length(tau))
  }
  if(dist == "weibull"){
    surv_param <- tail(ests, 2) # c(u,v): u=log(labmda); v=log(gamma)
  }
  
  pi <- piF(intercept, slope, x)
  Surv_L <- SurvF(L, tau, surv_param, beta, z, dist)
  # Surv_R <- SurvF(R, tau, surv_param, beta, z, dist)
  
  Ei <- pi / (pi + (1-pi)*Surv_L)
  Ei[delta==1] <- 1
  return(Ei)
}

QloglikF <- function(params, Ei, L, R, x, z, delta, dist){
  tau <- sort(unique(c(L, R)))
  tau <- tau[!is.infinite(tau)]
  n_cure <- ncol(x)
  n_beta <- ncol(z)
  
  intercept <- params[1]
  slope <- params[2:n_cure]
  beta <- params[(n_cure+1):(n_cure+n_beta)]
  if(dist == "npmle"){
    surv_param <- tail(params, length(tau))
  }
  if(dist == "weibull"){
    surv_param <- tail(params, 2) # c(u,v): u=log(labmda); v=log(gamma)
  }
  
  pi <- piF(intercept, slope, x)
  Surv_L <- SurvF(L, tau, surv_param, beta, z, dist)
  Surv_R <- SurvF(R, tau, surv_param, beta, z, dist)
  
  
  qloglik <- Ei*log((1-pi)*(Surv_L-Surv_R)) + (1-Ei)*log(pi)
  qloglik[is.infinite(qloglik)] <- 0
  sum(qloglik, na.rm = T)
}


gradF <- function(params, Ei, L, R, x, z, delta, dist){
  tau <- sort(unique(c(L, R)))
  tau <- tau[!is.infinite(tau)]
  n_cure <- ncol(x)
  n_beta <- ncol(z)
  
  intercept <- params[1]
  slope <- params[2:n_cure]
  beta <- params[(n_cure+1):(n_cure+n_beta)]
  if(dist == "npmle"){
    surv_param <- tail(params, length(tau))
  }
  if(dist == "weibull"){
    surv_param <- tail(params, 2) # c(u,v): u=log(labmda); v=log(gamma)
  }
  
  pi <- piF(intercept, slope, x)
  Surv_L <- SurvF(L, tau, surv_param, beta, z, dist)
  Surv_R <- SurvF(R, tau, surv_param, beta, z, dist)
  
  if(dist == "npmle"){
    mx <- 1-pi + pi*Surv_L
    dL_pi <- (1-delta) * (Surv_L-1)/mx
    dpi_alpha <- apply(x, 2, function(col) (pi - pi^2)*col)
    dL_i_alpha <- apply(dpi_alpha, 2, function(col) dL_pi*col)
    
    dL_i_lambda <- sapply(tau, function(tau_j){
      dL_S_ij <- delta/(Surv_L-Surv_R)*((L==tau_j)*1 - (R==tau_j)*1) +
        (1-delta)/mx*pi*(L==tau_j)
      dS_ij_lamda_j <- -SurvF(tau_j, tau, lambdas, beta, z)*exp(z%*%beta)
      dL_S_ij*dS_ij_lamda_j
    })
    
    lambdas_L <- sapply(L, function(i) lambdas[which(i==tau)])
    lambdas_R <- sapply(R, function(i) ifelse(is.infinite(i), 0, lambdas[which(i==tau)]))
    dL_i_beta <- apply(z, 2, function(col){
      delta/(Surv_L-Surv_R)*(Surv_L*(-lambdas_L) - Surv_R*(-lambdas_R) )*exp(z%*%beta)*col +
        (1-delta)/mx*pi*Surv_L*(-lambdas_L)*exp(z%*%beta)*col
    })
    
    dL_alpha <- colSums(dL_i_alpha)
    dL_lambda <- colSums(dL_i_lambda)
    dL_beta <- colSums(dL_i_beta)
    
    return(c(dL_alpha, dL_beta, dL_lambda))
  }
  if(dist == "weibull"){
    dQ_pi <- -(Ei/(1-pi)) + (1-Ei)/pi
    dpi_alpha <- apply(x, 2, function(col) (pi - pi^2)*col)
    dQ_i_alpha <- apply(dpi_alpha, 2, function(col) dQ_pi*col)
    
    dS_u <- function(t, tau, surv_param, beta, z, dist){
      surv <- SurvF(t, tau, surv_param, beta, z, dist)
      surv2 <- surv*log(surv)
      surv2[is.infinite(t)] <- 0
      surv2
      }

    dS_v <- function(t, tau, surv_param, beta, z, dist){
      surv <- SurvF(t, tau, surv_param, beta, z, dist)
      surv2 <- surv*log(surv)*log(t)*exp(surv_param[2])
      surv2[is.infinite(t)] <- 0
      surv2
    }

    dS_beta <- function(t, tau, surv_param, beta, z, dist){
      surv <- SurvF(t, tau, surv_param, beta, z, dist)
      surv2 <- apply(z, 2, function(col) surv*log(surv)*log(t)*col)
      surv2[is.infinite(t),] <- 0
      surv2
    }
    
    dQ_i_u <- Ei/(Surv_L-Surv_R)*(dS_u(t=L, tau, surv_param, beta, z, dist) - dS_u(t=R, tau, surv_param, beta, z, dist))
    dQ_i_v <- Ei/(Surv_L-Surv_R)*(dS_v(t=L, tau, surv_param, beta, z, dist) - dS_v(t=R, tau, surv_param, beta, z, dist))
    dQ_i_beta <- apply((dS_beta(t=L, tau, surv_param, beta, z, dist) - dS_beta(t=R, tau, surv_param, beta, z, dist)),
                       2,
                       function(col) Ei/(Surv_L-Surv_R)*(col))
    
    dQ_alpha <- colSums(dQ_i_alpha)
    dQ_beta <- colSums(dQ_i_beta)
    dQ_u <- sum(dQ_i_u)
    dQ_v <- sum(dQ_i_v)
    
    return(c(dQ_alpha, dQ_beta, dQ_u, dQ_v) )
  }
}

optimConstransF <- function(n_covm_cure = 0, 
                            n_covm_beta = 0, 
                            J = J){
  ####
  # calculte ui and ci for constrOptim()
  # ui %*% beta >= ci
  # Input:
  ## n_covm_cure: number of covariates for logistic model.
  ## J: total number of tests less or equal to t*.
  ####
  # construct constrains ui
  ui <- matrix(0, ncol = (n_covm_cure+n_covm_beta+J), nrow = (2*(n_covm_cure+n_covm_beta) + J))
  for (j in 1:(n_covm_cure+n_covm_beta)) {
    ui[(1+(j-1)*2):((j-1)*2+2), j] <- c(1, -1)
  } # alphas and betas
  ui_Sj <- diag(1, J)
  diag(ui_Sj[-1,]) <- -1
  ui[(nrow(ui)-J+1):nrow(ui), (ncol(ui)-J+1):ncol(ui)] <- ui_Sj # S_J
  # construct conztrains ci
  # ci <- c(rep(-1e+06, 2*(n_covm_cure+n_covm_beta)), -1, rep(0, J))
  ci <- c(rep(-10, 2*(n_covm_cure+n_covm_beta)), rep(0, J))

  out <- list(ui=ui, ci=ci)
  return(out)
}

# calculate p value
pvalF <- function(param_est, param_var, alpha = 0.05){
  z <- param_est / sqrt(param_var) # std from observed fisher information matrix, with a standard error given by the inverse of the square root of the Fisher information at the true value, i.e., 1/sqrt(nI_1(θ))
  2 * pnorm(-abs(z))
}


estF <- function(L, R, x, z, delta, bs_iter = NULL, dist){
  tau <- sort(unique(c(L, R)))
  tau <- tau[!is.infinite(tau)]
  J <- length(tau)
  if(dist == "npmle"){
    params_init <- sort(runif(ncol(x)+ncol(z)+J))
    constrains <- optimConstransF(n_covm_cure = ncol(x), 
                                  n_covm_beta = ncol(z), 
                                  J = J)
    ui <- constrains$ui
    ci <- constrains$ci
    param_ests <- constrOptim(theta = params_init, 
                              f = loglikF, 
                              grad = gradF,
                              ui = ui, ci = ci, control = list(fnscale = -1),
                              L = L, R = R, 
                              x = x, z = z, delta = delta, 
                              dist = dist,
                              hessian = TRUE)
    
    
    if(!is.null(bs_iter)){
      
      # registerDoParallel(cores=detectCores(all.tests=TRUE))
      registerDoParallel(cores=8)
      
      point_estimates_bs <- foreach(i=1:bs_iter, .combine=rbind, .errorhandling="remove") %dopar% {
        # resample with replacement
        bs_id <- sample(1:length(L), length(L), replace = T)
        L_bs = L[bs_id]
        R_bs = R[bs_id]
        x_bs = x[bs_id,]
        z_bs = z[bs_id,]
        delta_bs = delta[bs_id]
        
        tau_bs <- sort(unique(c(L_bs, R_bs)))
        tau_bs <- tau[!is.infinite(tau_bs)]
        J <- length(tau_bs)
        params_init <- sort(runif(ncol(x)+ncol(z)+J))
        # optimization
        bs_ests <- constrOptim(theta = params_init, 
                               f = loglikF, 
                               grad = gradF,
                               ui = ui, ci = ci, control = list(fnscale = -1),
                               L = L_bs, R = R_bs, 
                               x = x_bs, z = z_bs, delta = delta_bs,
                               hessian = F)$par
        bs_surv <- rep(NA, length(tau))
        bs_surv[tau %in% tau_bs] <- tail(bs_ests, J)
        bs_ests <- c(head(bs_ests, ncol(x)+ncol(z)), bs_surv)
        bs_ests
      }
      
      std_bs <- c(apply(point_estimates_bs, 2, function(col) sd(col, na.rm = T)))
    }
    
    return(list(par = param_ests$par, hessian = param_ests$hessian, std_bs = std_bs))
  }
  if(dist == "weibull"){
    
    ######################## Initialization ######################## 
    params_init <- c(rnorm(ncol(x)+ncol(z)), -0.0001, 0.002)
    Ei_init <- runif(nrow(x))
    Ei_init[delta==1] <- 1
    ests <- optim(par = params_init, fn = QloglikF, gr = gradF,
                  L = L, R = R, 
                  Ei = Ei_init,
                  x = x, z = z, delta = delta, 
                  dist = dist,
                  control = list(fnscale = -1),
                  hessian = T) # fnscale=-1 to find maximum
    
    ######################## End of Initialization ######################## 
    ## Maximization ####
    qloglik <- 0
    qloglik[2] <- QloglikF(params = ests$par, Ei = Ei_init, L, R, x, z, delta, dist)
    i <- 2
    
    tol = 1e-6 # tolerance 1e-6
    while (abs(qloglik[i] - qloglik[i - 1]) >= tol) {
      # update Ei using theta^{t-1}
      Ei <- EiF(ests$par, L, R, x, z, delta, dist)
      
      ests <- optim(par = ests$par, fn = QloglikF, gr = gradF,
                    L = L, R = R, 
                    Ei = Ei_init,
                    x = x, z = z, delta = delta, 
                    dist = dist,
                    control = list(fnscale = -1),
                    hessian = T) # fnscale=-1 to find maximum
      i <- i + 1
      qloglik[i] <- QloglikF(params = ests$par, Ei = Ei, L, R, x, z, delta, dist)
      # print(c(i, est_nparms, loglik[i]))
      print(paste0(i-2, "'s Iteration ..."))
    }
    print("Done!")
    maxLoglik <- loglikF(ests$par, L, R, x, z, delta, dist)
    
    return(list(iter = i, par = ests$par, hessian = ests$hessian, maxloglik = maxLoglik))
    
  }
}

mxcureF <- function(survfun, curefun, data, bs_iter, dist){
  sdata <- data
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  temp <- c("", "formula", "data", "na.action")
  mf <- mf[match(temp, names(mf), nomatch = 0)]
  mf[[1]] <- as.name("model.frame")
  avars <- all.vars(survfun)
  survfun1 <- as.formula(paste("Surv(", avars[1], ",", avars[2], 
                               ",type='interval2')~", paste(avars[-c(1, 2)], collapse = "+"), 
                               sep = ""))
  temp.z <- terms(survfun1, data = sdata)
  temp.x <- terms(curefun, data = sdata)
  mf$formula <- temp.z
  mf <- eval(mf, envir = parent.frame())
  Y <- model.extract(mf, "response")
  # attr(temp.z, "intercept") <- 0
  # attr(temp.x, "intercept") <- 1
  Zp <- model.matrix(temp.z, sdata)[,-1]
  # colnames(Zp) <- paste("z", 1:ncol(Zp), sep = "")
  Xp <- model.matrix(temp.x, sdata)
  # colnames(Xp) <- paste("x", c(1:ncol(Xp)) - 1, sep = "")
  sdata <- data.frame(L = Y[,1],
                      R = ifelse(Y[,2]==1, Inf, Y[,2]),
                      status = ifelse(Y[,3] == 3, 1, 0))
  # list(sdata = sdata, X = Xp, Z = Zp)
  ests <- estF(L = sdata$L, R = sdata$R, x = Xp, z = Zp, delta = sdata$status, bs_iter, dist)
  # ests <- estF(L = sdata$L, R = sdata$R, x = Xp, z = Zp, delta = sdata$status, bs_iter, dist = "npmle")
  # point_ests <- head(ests$par, ncol(Xp)+ncol(Zp))
  # point_ests_std <- head(ests$std_bs, ncol(Xp)+ncol(Zp))
  # names <- c(paste0('cure_', colnames(Xp)),
  #            paste0('beta_', colnames(Zp)))
  #            # paste0('Cum_haz_', with(sdata, sort(unique(c(L, R))[!is.infinite(unique(c(L, R)))]))))
  # tb_out <- data.frame(names = names,
  #                      Estimates = point_ests,
  #                      Std = point_ests_std,
  #                      pval = pvalF(param_est = point_ests,
  #                                   param_var = point_ests_std^2)
  # )
  # cumhaz <- tail(ests$par, length(ests$par)-ncol(Xp)-ncol(Zp))
  # cumhaz_std <- tail(ests$std_bs, length(ests$par)-ncol(Xp)-ncol(Zp))
  # cumhaz_lb <- cumhaz-1.96*cumhaz_std
  # cumhaz_ub <- cumhaz+1.96*cumhaz_std
  # surv <- exp(-cumhaz)
  # surv_lb <- exp(-cumhaz_lb)
  # surv_ub <- exp(-cumhaz_ub)
  # time <- with(sdata, sort(unique(c(L, R))[!is.infinite(unique(c(L, R)))]))
  # time <- ifelse(time==0, 0.1, time)
  # time <- c(0, time)
  # ssurv <- c(1, surv)
  # surv_lb <- c(1, surv_lb)
  # surv_ub <- c(1, surv_ub)
  # surv <- data.frame(survival = ssurv, time = time, lb = surv_lb, ub = surv_ub)
  # list(coef = tb_out, surv = surv)
}


