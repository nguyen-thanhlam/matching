# MHL (Mahalanobis Distance Matching)

library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")

#-------------------------------------------
#                 Simulation 
#-------------------------------------------
for (strat in strategy) {
  for (algo in algo_list) {
    for (n.cova in n.cova.list) {
      for (ttm.prop in ttm.prop.list) {
        for (event.prop in event.prop.list) {
          
          pvalue.match = c()
          pvalue.cmatch = c()
          clog_coef = c()
          clog_pval = c()
          clog_ci = c()

          scale.cens <- scale.res[(scale.res$n.cova==n.cova)&(scale.res$ttm.prop==ttm.prop)&(scale.res$event.prop==event.prop),'scale.cens']
          print(scale.cens)
          name = NA
          for (k in 1:n.cova) name = c(name,paste0("c", k))
          name = name[-1]
          
          for ( z in c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1))) ) {
            set.seed(z)
            cat("Mahalanobis distance, iter ",z,"/",nsim," \r")
            flush.console()
            
            n = ifelse(event.prop == 0.01, 20000*n.pair/200, 10000*n.pair/200)
            
            data = sim(n = n, 
                       hr.ttm = hr.ttm, 
                       px.str = px.str, 
                       ttm.prop = ttm.prop, 
                       scale.cens = scale.cens,
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
              
              if (strat == "matching") {
                dist_list1 <- dist_list
              } else {
                if (var(d$px) == 0) {
                  dist_list1 <- dist_list
                } else {
                  dist_list1 <- addrevcaliper(dist = dist_list, z = d$case, 
                                              dx = d$px, rg = c(-0.5, 0.5), 
                                              stdev = TRUE, penalty = max(dist_list$d))
                }
              }
              
              o1 <- match(z = d$case, dist = dist_list1, dat = d, ncontrol = 1)
              
              if (o1$feasible == TRUE && nrow(o1$data) >= 2) {
                ctrl.i.id = o1$data$id[2]
                sel = rbind(sel, c(source.noc$id[i], ctrl.i.id, j))
                j = j + 1
              }
            }
            
            
            if (nrow(sel) > 1) {
              selected = as.data.frame(sel[-1, , drop = FALSE]) 
              colnames(selected) = c('id.case', 'id.ctrl', 'pair')
              
              case.noc = merge(source.noc, selected[, c('id.case', 'pair')], by.x='id', by.y='id.case')
              ctrl.noc = merge(source.noc, selected[, c('id.ctrl', 'pair')], by.x='id', by.y='id.ctrl')
              ctrl.noc$stt = 0
              
              
              case.noc = case.noc[order(case.noc$pair), ]
              ctrl.noc = ctrl.noc[order(ctrl.noc$pair), ]
              
              # MH test
              pvalue.match = c(pvalue.match, sens.analysis.mh(cases.exposed = case.noc$ttm,
                                                              referents.exposed = ctrl.noc$ttm,
                                                              no.referents = 1, Gamma = 1)$lower.bound.pval)
              
              # Clogit
              clog_dat <- rbind(case.noc, ctrl.noc)
              clog_dat <- clog_dat[order(clog_dat$pair),]
              
              #print(mean(clog_dat$ttm))
              #print(mean(clog_dat$px))
              
              clog_mod <- clogit(stt ~ ttm + strata(pair), 
                                 data = clog_dat, 
                                 method = "breslow")
              clog_sum <- summary(clog_mod)
              
              clog_coef <- c(clog_coef, clog_sum$coefficients[1, 1])
              clog_pval <- c(clog_pval, clog_sum$coefficients[1, 5])
              clog_ci <- rbind(clog_ci,confint(clog_mod))  
            }
          }
          
          mhl_res <- do.call(rbind, list(mhl_res,
          data.frame(
            strategy  = strat,
            algo      = algo,
            approach  = "mhl",
            n         = n,
            n.cova    = n.cova,
            event.prob = event.prop,
            ttm.prob  = ttm.prop,
            seed = c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1))),
            mh_pval   = pvalue.match,
            clog_coef = clog_coef,
            clog_pval = clog_pval,
            clog_lower = clog_ci[,1],
            clog_upper = clog_ci[,2])))
        }  
      }
    }
  }
}


