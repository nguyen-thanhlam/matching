# MHL (Mahalanobis Distance Matching)

library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")

#--------------------
# Simulation function
#--------------------

sim <- function(n = 10000, a = 2, med0 = 200, hr.ttm = 1, hr.c = 1.5, 
                px.str = 3, scale.cens = 100, ttm.prop = 0.1, n.cova = 20) {
  
  px <- rbinom(n, 1, 0.5)
  
  cova <- matrix(0, ncol = n.cova, nrow = n)
  name <- NA
  
  for (k in 1:n.cova) {
    cova[, k] <- runif(n, -0.5, 0.5)
    name <- c(name, paste0("c", k))
  }
  colnames(cova) <- name[-1]
  
  cova.coef.ttm <- rep(c(-1, +1), length.out = n.cova)
  cova.lp.ttm <- as.vector(cova %*% cova.coef.ttm)
  
  itc.ttm <- qnorm(ttm.prop) * sqrt(1 + n.cova / 12 + px.str^2 * 0.25) - px.str * 0.5
  lp <- pnorm(itc.ttm + px.str * px + cova.lp.ttm)
  ttm <- rbinom(n, 1, lp)
  
  b0 <- med0 * log(2)^(-1 / a)
  b.ttm <- log(hr.ttm)  
  
  cova.coef.oc <- rep(c(log(hr.c), -log(hr.c)), length.out = n.cova)
  cova.lp.oc <- as.vector(cova %*% cova.coef.oc)
  t.event <- rweibull(n = n, shape = a, scale = b0 * exp((-b.ttm * ttm - cova.lp.oc) / a))
  t.cens <- rweibull(n = n, shape = a, scale = scale.cens)
  
  tobs <- pmin(t.event, t.cens)
  out <- data.frame(px, ttm, cova,
                    tobs,
                    stt = as.numeric(tobs == t.event),
                    id = 1:n)
  return(out)
}


for (algo in algo_list) {
  for (nc in n.cova.list) {
    for (scale in scale.cens.list) {
      for (ttm in ttm.prop.list) {
    
    pvalue.match = c()
    pvalue.cmatch = c()
    clog_coef = c()
    clog_pval = c()
    
    n.cova <- nc
    scale.cens <- scale
    ttm.prop <- ttm
    
    name = NA
    for (k in 1:n.cova) name = c(name,paste0("c", k))
    name = name[-1]
    
    for (z in 1:nsim) {
      set.seed(z)
      cat("mahalanobis distance, approach ", algo, nc, " covariates ",
          ", event.prop ", ifelse(scale==26.3,0.01,ifelse(scale==55.2,0.05,0.10)),
          ", ttm.prop ", ttm,
          ", iter ",z,"/",nsim,
          " \r")
      flush.console()
      
      n = ifelse(scale.cens == 26.3, 20000, 10000)
      
      data = sim(n=n, 
                 hr.ttm=hr.ttm, 
                 px.str=px.str, 
                 ttm.prop = ttm.prop, 
                 scale.cens=scale.cens,
                 n.cova = n.cova)
      #print(mean(data$ttm))
      #print(mean(data$stt))
      # -------------------------
      # No counter-matching (noc)
      # ------------------------- 
      source.noc = data[order(data$tobs),]
      j = 1
      
      sel = matrix(c(-1, -1, 0), nrow = 1) 
      
      ind = which(source.noc$stt == 1)
      
      for (i in ind) {
        if (j > n.pair) break
        if (i == n) break
        
        if (algo == 'A1') {
          elig = setdiff(source.noc$id[i:n], unique(c(sel[,1], sel[,2])))
          if (!source.noc$id[i] %in% elig | length(elig) == 1) next
          d = source.noc[source.noc$id %in% elig, ]
          d = d[order(d$tobs), ]
          d$case = c(1, rep(0, nrow(d) - 1))
        } else {
          
          elig = c(i:n)
          #elig = elig[!(source.noc$id[elig] %in% unique(c(sel[,1], sel[,2])))]
          if (length(elig) == 1) next
          
          d = source.noc[elig, ]
          d$case = c(1, rep(0, nrow(d) - 1))
        }
        
        # Mahalanobis 
        X_case <- d[1, name, drop = FALSE]
        X_ctrl <- d[-1, name, drop = FALSE]
        
        if (nrow(X_ctrl) < 2) next
        cov_ctrl <- cov(X_ctrl)
        
        # Check if covariance is singular, use generalized inverse if needed
        if (rcond(cov_ctrl) < 1e-10) {
          cov_ctrl_inv <- MASS::ginv(cov_ctrl)
          dist_vector <- mahalanobis(x = as.matrix(X_ctrl), 
                                     center = as.numeric(X_case), 
                                     cov = cov_ctrl_inv,
                                     inverted = TRUE)
        } else {
          dist_vector <- mahalanobis(x = as.matrix(X_ctrl), 
                                     center = as.numeric(X_case), 
                                     cov = cov_ctrl)
        }
        dist_list <- list(
          start = rep(1, nrow(X_ctrl)),
          end   = 2:(nrow(X_ctrl) + 1),
          d     = as.numeric(dist_vector)
        )
        
        o1 <- match(z = d$case, dist = dist_list, dat = d, ncontrol = 1)
        
        if (o1$feasible == TRUE && nrow(o1$data) >= 2) {
          ctrl.i.id = o1$data$id[2]
          sel = rbind(sel, c(source.noc$id[i], ctrl.i.id, j))
          j = j + 1
        }
      }
      
      
      if (nrow(sel) > 1) {
        selected = as.data.frame(sel[-1, , drop = FALSE]) 
        colnames(selected) = c('id.case', 'id.ctrl', 'pair')
        
        test = source.noc
        test = merge(test, selected[,c('id.case','pair')], by.x='id', by.y='id.case', all.x=TRUE)
        test = merge(test, selected[,c('id.ctrl','pair')], by.x='id', by.y='id.ctrl', all.x=TRUE)
        test$pair = ifelse(is.na(test$pair.x), 0, test$pair.x) + ifelse(is.na(test$pair.y), 0, test$pair.y)
        
        case.noc = test[test$pair.x != 0 & !is.na(test$pair.x), ]
        ctrl.noc = test[test$pair.y != 0 & !is.na(test$pair.y), ]
        
        # MH test
        pvalue.match[z] = sens.analysis.mh(cases.exposed = case.noc$ttm,
                                           referents.exposed = ctrl.noc$ttm,
                                           no.referents = 1, Gamma = 1)$lower.bound.pval
        
        # Clogit
        clog_dat <- test[test$pair != 0, ]
        clog_dat <- clog_dat[order(clog_dat$pair),]
        
        clog_mod <- clogit(stt ~ ttm + strata(pair), 
                           data = clog_dat, 
                           method = "breslow")
        clog_sum <- summary(clog_mod)
        
        clog_coef[z] <- clog_sum$coefficients[1, 1]
        clog_pval[z] <- clog_sum$coefficients[1, 5]
      }
    }
    
    tmp <- paste(algo, nc, scale, ttm, sep = ".")
    
    mhl_res[[tmp]] <- data.frame(
      algo      = algo,
      approach  = "mhl",
      n         = n,
      n.cova    = nc,
      event.prob = ifelse(scale==26.3,1,ifelse(scale==55.2,5,10)),
      ttm.prob  = ttm,
      mh_pval   = pvalue.match,
      clog_coef = clog_coef,
      clog_pval = clog_pval)
      }  
    }
  }
}

res <- do.call(rbind, mhl_res)
