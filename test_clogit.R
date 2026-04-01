library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")

#--------------------
# Simulation study parameters
#--------------------

nphase = 3
nsim = 3
n.pair = 200
n = 20000
px.str = 3
hr.ttm = 1

# Treatment prob = 10% and 50%
#ttm.prop.list = c(0.1, 0.5)
ttm.prop.list = c(0.1)
# Scale cens --> event prob 1%, 5%, 10%
#scale.cens.list = c(26.3, 55.2, 79.6)
scale.cens.list = c(26.3,79.6)
# Nb of cova
#n.cova.list = c(4, 10, 20)
n.cova.list = c(1)
n.hx = n
conv = 0
# List of algorithms
algo_list = c('A1', 'A2')


sim <- function(n = 10000, a = 2, med0 = 200, hr.ttm = 1, hr.c = 1.5, 
                px.str = 3, scale.cens = 100, ttm.prop = 0.1, n.cova = 20) {
  
  px <- rbinom(n, 1, 0.5)
  
  cova <- matrix(0, ncol = n.cova, nrow = n)
  name <- NA
  
  for (k in 1:n.cova) {
    cova[, k] <- rbinom(n, 1, 0.3)
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

res = c()
res_ex = c()
res_nm = c()
for (iphase in c(1:nphase)) {
for (algo in algo_list) {
    for (nc in n.cova.list) {
        for (scale in scale.cens.list) {
            for (ttm in ttm.prop.list) {
                clog_coef = c()
                clog_pval = c()
                clog_ci = c()
                clog_coef_nm = c()
                clog_pval_nm = c()
                clog_ci_nm = c()
                
                n.cova <- nc
                scale.cens <- scale
                ttm.prop <- ttm
                
                name = NA
                for (k in 1:n.cova) name = c(name,paste0("c", k))
                name = name[-1]
                
                for (z in c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1)))) {
                    set.seed(z)
                    cat("Test", algo, nc, " covariates ",
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
                    sel_nm = matrix(c(-1, -1, 0), nrow = 1) 

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
                            if (length(elig) == 1) next
                            d = source.noc[elig, ]
                            d$case = c(1, rep(0, nrow(d) - 1))
                        }
                        
                        # test
                        X_case <- d[1, name, drop = FALSE]
                        X_ctrl <- d[-1, name, drop = FALSE]
                        id_ctrl <- d[-1, 'id', drop = FALSE]
                        # exact match
                        X_diff <- X_ctrl - c(X_case)
                        zero_indices <- which(rowSums(X_diff == 0) == ncol(X_diff))
                        random_index <- sample(zero_indices, 1)

                        ctrl.i.id = id_ctrl[random_index,'id']
                        sel = rbind(sel, c(source.noc$id[i], ctrl.i.id, j))
                        
                        # no match
                        
                        sel_nm = rbind(sel_nm, c(source.noc$id[i], sample(id_ctrl[['id']], 1), j))
                        j = j + 1
                    }
                    
                    lst_clog_dat = list()
                    lst_sel = list(sel,sel_nm)
                    for (i in 1:length(lst_sel)) {
                        dat = lst_sel[[i]]
                        if (nrow(dat) > 1) {
                            selected = as.data.frame(dat[-1, , drop = FALSE]) 
                            colnames(selected) = c('id.case', 'id.ctrl', 'pair')
                            case.noc = merge(source.noc, selected[, c('id.case', 'pair')], by.x='id', by.y='id.case')
                            ctrl.noc = merge(source.noc, selected[, c('id.ctrl', 'pair')], by.x='id', by.y='id.ctrl')
                            ctrl.noc$stt = 0
                            case.noc = case.noc[order(case.noc$pair), ]
                            ctrl.noc = ctrl.noc[order(ctrl.noc$pair), ]
                
                            clog_dat <- rbind(case.noc, ctrl.noc)
                            clog_dat <- clog_dat[order(clog_dat$pair),]
                        } else {
                            clog_dat = data.frame()}
                        lst_clog_dat[[i]] = clog_dat
                    }
                        # Clogit exact match
                    clog_mod <- clogit(stt ~ ttm + strata(pair), 
                                        data = lst_clog_dat[[1]], 
                                        method = "breslow")
                    clog_sum <- summary(clog_mod)
                    
                    clog_coef <- c(clog_coef,clog_sum$coefficients[1, 1])
                    clog_pval <- c(clog_pval,clog_sum$coefficients[1, 5])
                    clog_ci <- rbind(clog_ci,confint(clog_mod))                    
                    # Clogit no match
                    clog_mod_nm <- clogit(as.formula(paste0('stt ~ ttm + ',
                                                    paste(name,collapse = " + "),
                                                    ' + strata(pair)')), 
                                        data = lst_clog_dat[[2]], 
                                        method = "breslow")
                    clog_sum_nm <- summary(clog_mod_nm)
                    
                    clog_coef_nm <- c(clog_coef_nm,clog_sum_nm$coefficients[1, 1])
                    clog_pval_nm <- c(clog_pval_nm,clog_sum_nm$coefficients[1, 5])
                    clog_ci_nm <- rbind(clog_ci_nm,confint(clog_mod_nm)['ttm',])
                }
                res_ex <- rbind(res_ex,
                data.frame(
                    algo      = algo,
                    approach  = 'exact',
                    n         = n,
                    n.cova    = nc,
                    event.prob = ifelse(scale==26.3,1,ifelse(scale==55.2,5,10)),
                    ttm.prob  = ttm,
                    seed = c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1))),
                    clog_coef = clog_coef,
                    clog_pval = clog_pval,
                    clog_lower = clog_ci[,1],
                    clog_upper = clog_ci[,2]))

                res_nm <- rbind(res_nm,
                data.frame(
                    algo      = algo,
                    approach  = 'no_match',
                    n         = n,
                    n.cova    = nc,
                    event.prob = ifelse(scale==26.3,1,ifelse(scale==55.2,5,10)),
                    ttm.prob  = ttm,
                    seed = c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1))),
                    clog_coef = clog_coef_nm,
                    clog_pval = clog_pval_nm,
                    clog_lower = clog_ci_nm[,1],
                    clog_upper = clog_ci_nm[,2]))
            }  
        }
    }
}
    res <- do.call(rbind, list(res, rbind(res_ex,res_nm)))
    View(res)
    write.csv(res, paste0("match_res_clogtest.csv"))
}


