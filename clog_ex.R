# ex (EXact Matching)
library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")
#-------------------------------------------
#                 Simulation 
#-------------------------------------------
for (iparams in c(1:nrow(params))) {
    strat = params[iparams,"strat"]
    algo = params[iparams,"algo"]
    n.cova = params[iparams,"n.cova"]
    ttm.prop = params[iparams,"ttm.prop"]
    event.prop = params[iparams,"event.prop"]
    r.penalty = params[iparams,"r.penalty"]
    n.pair = params[iparams,"n.pair"]

    pvalue.match = c()
    pvalue.cmatch = c()
    clog_coef = c()
    clog_pval = c()
    clog_ci = c()
    #cnt.match.prob = c()

    scale.cens <- scale.res[(scale.res$n.cova==n.cova)&(scale.res$ttm.prop==ttm.prop)&(scale.res$event.prop==event.prop),'scale.cens']
    name = NA
    for (k in 1:n.cova) name = c(name,paste0("c", k))
    name = name[-1]
    
    for ( z in split(c(1:nsim), ceiling(seq_along(c(1:nsim))/ceiling(nsim/nworkers)))[[iphase]] ) {
      #cnt.match = 0
      set.seed(z)

      n = ifelse(event.prop == 0.01, 20000*n.pair/200, 10000*n.pair/200)
      data = sim(n = n, 
                  hr.ttm = hr.ttm, 
                  px.str = px.str, 
                  ttm.prop = ttm.prop, 
                  scale.cens = scale.cens,
                  n.cova = n.cova)
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
          d$case = c(1, rep(0, nrow(d) - 1))
        } else {          
          elig = c(i:n)
          if (length(elig) == 1) next          
          d = source.noc[elig, ]
          d$case = c(1, rep(0, nrow(d) - 1))
        }        
        # Mahalanobis 
        X_case <- d[1, name, drop = FALSE]
        X_ctrl <- d[-1, name, drop = FALSE]        
        id_ctrl <- d[-1, 'id', drop = FALSE]

        X_diff <- X_ctrl - c(X_case)
        X_zero <- which(rowSums(X_diff == 0) == ncol(X_diff))
        if (strat=="match") {
            random_index <- sample(X_zero, 1)
        } else {
            px_case = d[1, 'px', drop=FALSE]
            px_ctrl = d[-1, 'px', drop=FALSE]
            px_diff = px_ctrl - c(px_case)
            px_zero = which(rowSums(px_diff != 0) == ncol(px_diff))
            common = intersect(X_zero, px_zero)
            random_index <- sample(common, 1)
        }
        ctrl.i.id = id_ctrl[random_index,'id']
        sel = rbind(sel, c(source.noc$id[i], ctrl.i.id, j))
        j = j + 1
      }
      #cnt.match.prob = c(cnt.match.prob, cnt.match / (nrow(sel)-1))

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
        if (algo=="A1") {
            clog_mod <- clogit(stt ~ ttm + strata(pair), data = clog_dat, method = "breslow")
        } else {
            clog_mod <- clogit(stt ~ ttm + strata(pair) + cluster(id), data = clog_dat, method = "breslow")
        }             
        clog_sum <- summary(clog_mod)             
        clog_coef <- c(clog_coef, clog_sum$coefficients[1, 1])
        clog_pval <- c(clog_pval, clog_sum$coefficients[1, 5])
        clog_ci <- rbind(clog_ci,confint(clog_mod))  
      }
    }
    
    ex_res <- do.call(rbind, list(ex_res,
    data.frame(
      strategy  = strat,
      algo      = algo,
      approach  = "ex",
      n         = n,
      n.cova    = n.cova,
      event.prob = event.prop,
      ttm.prob  = ttm.prop,
      r.penalty = r.penalty,
      #cnt.match.prob = cnt.match.prob,
      seed = split(c(1:nsim), ceiling(seq_along(c(1:nsim))/ceiling(nsim/nworkers)))[[iphase]],
      mh_pval   = pvalue.match,
      clog_coef = clog_coef,
      clog_pval = clog_pval,
      clog_lower = clog_ci[,1],
      clog_upper = clog_ci[,2])))
}  


