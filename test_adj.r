library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")
library(dplyr)
#--------------------
# Simulation study parameters
#--------------------

nphase = c(8:10)
nsim = 100
n.pair.list = c(200,500)
n = 20000
px.str = 3
hr.ttm = 1

# Treatment prob = 10% and 50%
#ttm.prop.list = c(0.1, 0.5)
ttm.prop.list = c(0.1)
# Scale cens --> event prob 1%, 5%, 10%
#scale.cens.list = c(26.3, 55.2, 79.6)
scale.cens.list = c(26.3, 79.6)
# Nb of cova
#n.cova.list = c(4, 10, 20)
n.cova.list = c(4)
n.hx = n
conv = 0
# List of algorithms
algo_list = c('A2')


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

addrevcaliper<-function(dist,z,dx,rg, stdev = FALSE, penalty = 1000){
  
  stopifnot(is.vector(rg)&(length(rg)==2))
  stopifnot((rg[1]<=0)&(rg[2]>=0))
  
  
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  
  dx1<-dx[z==1]
  dx0<-dx[z==0]
  
  v <- stats::var(dx)
  sp <- sqrt(v)
  stopifnot(sp>0)
  
  if (stdev) rg <- rg *sp
  
  m=sum(z)
  dif<-dx1[dist$start]-dx0[dist$end-m]
  d <- dist$d + as.numeric(dif<rg[2])*penalty
  d <- d + as.numeric(dif>rg[1])*penalty
  
  list(d=d,start=dist$start,end=dist$end)
}

func_m_cm = function(d, name, method) {
    X_case <- d[1, name, drop = FALSE]
    X_ctrl <- d[-1, name, drop = FALSE]
    id_ctrl <- d[-1, 'id', drop = FALSE]
    if (method == "mhl") {
        d$mhl = NA
        if (nrow(X_ctrl) < 2) next
        cov_ctrl <- cov(X_ctrl)
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
    } else if (method == "ps") {
        X = d[, "ps"]
        distance <- abs(log(rep(X[1],length(X)-1)/(1-rep(X[1],length(X)-1))) 
                        - log(X[2:length(X)]/(1 - X[2:length(X)])))
        
        dist_list <- list()
        dist_list$start = rep(1,length(X)-1)
        dist_list$end = c(2:length(X))
        dist_list$d = distance
    } else {
        X = d[, "drs"]
        distance <- abs(log(rep(X[1],length(X)-1)/(1-rep(X[1],length(X)-1))) 
                        - log(X[2:length(X)]/(1 - X[2:length(X)])))
        dist_list <- list()
        dist_list$start = rep(1,length(X)-1)
        dist_list$end = c(2:length(X))
        dist_list$d = distance
    }
        
    dist_list1_m <- dist_list
    if (var(d$px) == 0) {
        dist_list1_cm <- dist_list
    } else {
        dist_list1_cm <- addrevcaliper(dist = dist_list, z = d$case, 
                                        dx = d$px, rg = c(-0.5, 0.5), 
                                        stdev = TRUE, penalty = max(dist_list$d))
    }
    o1_m <- match(z = d$case, dist = dist_list1_m, dat = d, ncontrol = 1)
    o1_cm <- match(z = d$case, dist = dist_list1_cm, dat = d, ncontrol = 1)
    if (o1_m$feasible == TRUE && nrow(o1_m$data) >= 2) {
        sel_m = c(o1_m$data$id)
    }
    if (o1_cm$feasible == TRUE && nrow(o1_cm$data) >= 2) {
        sel_cm = c(o1_cm$data$id)
    }
return(c(sel_m,sel_cm))}
safe_clogit <- function(formula, data, method) {
  warn_flag <- FALSE
  err_flag <- FALSE

  fit <- withCallingHandlers(
    tryCatch(
      {
        clogit(formula, data = data, method = method)
      },
      error = function(e) {
        err_flag <<- TRUE
        return(NULL)
      }
    ),
    warning = function(w) {
      warn_flag <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  list(
    model = fit,
    summary = if (!is.null(fit)) summary(fit) else NULL,
    converge = ifelse(warn_flag || err_flag || is.null(fit), 0, 1)
  )
}
proc_clogit = function(dat, source.dat, typ, name) {
    selected = as.data.frame(dat[-1, , drop = FALSE]) 
    colnames(selected) = c('id.case', 'id.ctrl', 'pair')
    case.noc = merge(source.dat, selected[, c('id.case', 'pair')], by.x='id', by.y='id.case')
    ctrl.noc = merge(source.dat, selected[, c('id.ctrl', 'pair')], by.x='id', by.y='id.ctrl')
    ctrl.noc$stt = 0
    case.noc = case.noc[order(case.noc$pair), ]
    ctrl.noc = ctrl.noc[order(ctrl.noc$pair), ]
    clog_dat <- rbind(case.noc, ctrl.noc)
    clog_dat <- clog_dat[order(clog_dat$pair),]
    if ((typ=='drs'|typ=='ps')) {
        clog_mod_score <- safe_clogit(as.formula(paste0('stt ~ ttm +',
                                        typ,
                                        ' + strata(pair)')),
                            data = clog_dat, 
                            method = "breslow")
        clog_sum_score <- clog_mod_score$summary
        if (!(is.null(clog_sum_score))) {
            clog_coef_score <- clog_sum_score$coefficients[1, 1]
            clog_pval_score <- clog_sum_score$coefficients[1, 5]
            clog_ci_score <- c(confint(clog_mod_score$model)['ttm',])
        } else {
            clog_coef_score <- NA
            clog_pval_score <- NA
            clog_ci_score <- c(NA, NA)
        }
        vec_score <- c(clog_coef_score, clog_pval_score, clog_ci_score, clog_mod_score$converge)
    } else (vec_score <- c(NA, NA, NA, NA, 0))
    clog_mod_cova <- safe_clogit(as.formula(paste0('stt ~ ttm + ',
                                        paste(name,collapse = " + "),
                                        ' + strata(pair)')), 
                            data = clog_dat, 
                            method = "breslow")
    clog_sum_cova <- clog_mod_cova$summary
        
    clog_coef_cova <- clog_sum_cova$coefficients[1, 1]
    clog_pval_cova <- clog_sum_cova$coefficients[1, 5]
    clog_ci_cova <- c(confint(clog_mod_cova$model)['ttm',])
    vec_cova <- c(clog_coef_cova, clog_pval_cova, clog_ci_cova, clog_mod_cova$converge)
return(c(vec_score, vec_cova))}
methods = c("mhl","ps","drs")
fin.res = c()
for (iphase in nphase) {
    for (n.pair in n.pair.list) {
        for (algo in algo_list) {
            for (n.cova in n.cova.list) {
                for (scale.cens in scale.cens.list) {
                    for (ttm.prop in ttm.prop.list) {
                        
                        name = NA
                        for (k in 1:n.cova) name = c(name,paste0("c", k))
                        name = name[-1]
                        
                        for (z in c((1+nsim*(iphase-1)):(nsim+nsim*(iphase-1)))) {
                            set.seed(z)
                            cat("Test", algo, ", iter ",z,"/",nsim," \r")
                            flush.console()
                            
                            n = ifelse(scale.cens == 26.3, 20000*n.pair/200, 10000*n.pair/200)
                            
                            data = sim(n=n, 
                                    hr.ttm=hr.ttm, 
                                    px.str=px.str, 
                                    ttm.prop = ttm.prop, 
                                    scale.cens=scale.cens,
                                    n.cova = n.cova)
                            # PS
                            ps.dat = data[which(data$px==1),c("ttm",name)]
                            mod = glm(ttm ~ ., 
                                    family = binomial(link = "probit"), 
                                    data = ps.dat)
                            data$ps = predict(mod, 
                                            newdata = data[,name],
                                            type = 'response')
                            # DRS
                            hx = sim(n = n, hr.ttm = hr.ttm, 
                                    px.str = px.str, ttm.prop = ttm.prop, 
                                    scale.cens = scale.cens, n.cova = n.cova)
                            hx = hx[,c(name,"tobs","stt")]
                            mod = coxph(Surv(tobs,stt)~., data = hx, method = 'breslow')
                            time = max(hx$tobs)
                            null = as.data.frame(matrix(0,nrow=1,ncol=n.cova))
                            names(null) = name
                            fit_surv <- survfit(mod, newdata = null)
                            try = data.frame(time = fit_surv$time, surv= fit_surv$surv)
                            s0 = try$surv[which(try$time == max(try$time))]
                            lp= as.vector(as.matrix(data)[,name] %*% coef(mod)) 
                            data$drs = s0^exp(lp)
                            # -------------------------
                            # No counter-matching (noc)
                            # ------------------------- 
                            source.noc = data[order(data$tobs),]
                            j = 1
                            
                            sel_m_cm = list(matrix(c(-1, -1, -1, -1, 0), nrow = 1),
                                            matrix(c(-1, -1, -1, -1, 0), nrow = 1),
                                            matrix(c(-1, -1, -1, -1, 0), nrow = 1))
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
                                sel = lapply(methods, function(m) func_m_cm(d = d, name = name, method = m))
                                sel_m_cm = lapply(seq(length(methods)), function(i) {rbind(sel_m_cm[[i]], c(sel[[i]],j))})
                                # no match
                                id_ctrl <- d[-1, 'id', drop = FALSE]
                                sel_nm = rbind(sel_nm, c(source.noc$id[i], sample(id_ctrl[['id']], 1), j))
                                j = j + 1
                            }
                            sel_m_cm = lapply(c(0,2), function(j) lapply(seq(length(methods)), function(i) sel_m_cm[[i]][,c(c(1,2)+j,5)]))
                            sel_m = sel_m_cm[[1]]
                            sel_cm = sel_m_cm[[2]]
                            lst_sel = list(sel_m,sel_cm,sel_nm)
                            res = lapply(1:length(lst_sel), function(i) {
                                dat = lst_sel[[i]]
                                if (is.list(dat)) {
                                    names(dat) = methods
                                    res = lapply(methods, function(m) proc_clogit(dat = dat[[m]], source = source.noc, typ = m, name = name))
                                    names(res) = methods
                                } else {
                                    res = list(proc_clogit(dat = dat, source.dat = source.noc, typ = 'nm', name = name))
                                    names(res) = 'nm'
                                }
                                return(res)
                            })
                            names(res) = c('m','cm','nm')
                            for (mth in c('m','cm','nm')) {
                                for (typ in names(res[[mth]])) {
                                    res.vec = c(res[[mth]][[typ]],mth,typ,n.pair,algo,n.cova,ifelse(scale.cens==26.3,1,ifelse(scale.cens==55.2,5,10)),ttm.prop,z)
                                    fin.res[[mth]][[typ]] <- do.call(rbind,list(fin.res[[mth]][[typ]], res.vec))
                                    fin.res[[mth]][[typ]] = as.data.frame(fin.res[[mth]][[typ]])
                                    colnames(fin.res[[mth]][[typ]]) = c('coef.sc','pval.sc','lw.sc','up.sc','conv.sc','coef.cv','pval.cv','lw.cv','up.cv','conv.cv','method','type','n.pair','algo','n.cova','event.prob','ttm.prop','seed')
                                }
                            }
                        }
                        fin.out = as.data.frame(do.call(rbind,unlist(fin.res, recursive = FALSE)))
                    }  
                }
            }
        }
    }
    View(fin.out)
    write.csv(fin.out, paste0("match_res_clogadj_phase810.csv"))
}


