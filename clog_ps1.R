# Propensity Score
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
    conv = 0

    scale.cens <- scale.res[(scale.res$n.cova==n.cova)&(scale.res$ttm.prop==ttm.prop)&(scale.res$event.prop==event.prop),'scale.cens']

    name = NA
    for (k in 1:n.cova) name = c(name,paste0("c", k))
    name = name[-1]

    for ( z in split(c(1:nsim), ceiling(seq_along(c(1:nsim))/ceiling(nsim/nworkers)))[[iphase]] ) {
      set.seed(z)
      # Population 2*1e4 when event prob == 1%, else 1e4
      n = ifelse(event.prop == 0.01, 20000*n.pair/200, 10000*n.pair/200)
      data = sim(n = n, 
                  hr.ttm = hr.ttm, 
                  px.str = px.str, 
                  ttm.prop = ttm.prop, 
                  scale.cens = scale.cens, 
                  n.cova = n.cova)
      #---------------------------------------
      # Model construction for PS P(D=1|C,P=0)
      #---------------------------------------
      ps.dat = data[which(data$px==1),c("ttm",name)]
      mod = glm(ttm ~ ., family = binomial(link = "probit"), data = ps.dat)
      if (!isTRUE(mod$converged)) {
        pvalue.match[z] = NA
        clog_coef[z] = NA
        clog_pval[z] = NA
        next 
      } else {
        conv = conv + 1
      }
      # PS for all individuals 
      data$ps = predict(mod, newdata = data[,name], type = 'response')
  
      source.noc = data[order(data$tobs),]
      j = 1
      sel = matrix(c(-1,-1,0),nrow=1)
      ind = which(source.noc$stt ==1)
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
        X = d[, "ps"]
        distance <- abs(log(rep(X[1],length(X)-1)/(1-rep(X[1],length(X)-1))) 
                        - log(X[2:length(X)]/(1 - X[2:length(X)])))
        dist_list <- list()
        dist_list$start = rep(1,length(X)-1)
        dist_list$end = c(2:length(X))
        dist_list$d = distance
        
        if (strat == "matching") {
          dist_list1 <- dist_list
        } else {
          if (var(d$px) == 0) {
            next
          } else {
            dist_list1 <- addrevcaliper(dist = dist_list, z = d$case, 
                                        dx = d$px, rg = c(-0.5, 0.5), 
                                        stdev = TRUE, penalty = max(dist_list$d), r.penalty=r.penalty)
          }
        }
        o1 <- DiPs::match(z = d$case, dist = dist_list1, dat = d, ncontrol=1)

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
        
        # MH test on the first 200 matched sets
        pvalue.match = c(pvalue.match, sens.analysis.mh(cases.exposed = case.noc$ttm,
                                                          referents.exposed = ctrl.noc$ttm,
                                                          no.referents = 1, Gamma = 1)$lower.bound.pval)
        
        
        # Clogit
        clog_dat <- rbind(case.noc, ctrl.noc)
        clog_dat <- clog_dat[order(clog_dat$pair),]
        clog_mod <- clogit(stt ~ ttm + strata(pair), data = clog_dat, method = "breslow")
        
        clog_sum <- summary(clog_mod)        
        clog_coef <- c(clog_coef, clog_sum$coefficients[1, 1])
        clog_pval <- c(clog_pval, clog_sum$coefficients[1, 5])
        clog_ci <- rbind(clog_ci,confint(clog_mod))  
      }
    }

            
    ps_res <- do.call(rbind, list(ps_res,
      data.frame(
        strategy  = strat,
        algo      = algo,
        approach  = "ps",
        n         = n,
        n.cova    = n.cova,
        event.prob = event.prop,
        ttm.prob  = ttm.prop,
        r.penalty = r.penalty,
        seed = split(c(1:nsim), ceiling(seq_along(c(1:nsim))/ceiling(nsim/nworkers)))[[iphase]],
        mh_pval   = pvalue.match,
        clog_coef = clog_coef,
        clog_pval = clog_pval,
        clog_lower = clog_ci[,1],
        clog_upper = clog_ci[,2])))
    }