library(parallel)

# Simulation parameters
nsim = 100
n = 20000
px.str = 3
hr.ttm = 1.5
strategy = c("matching", "counter-matching")
algo_list = c('A1', 'A2')
n.pair.list = c(200)
ttm.prop.list = c(0.1, 0.5)
event.prop.list = c(0.01, 0.05, 0.10)
n.cova.list = c(2)
#r.penalty.list = c(0.5, 1, 5, 10)
r.penalty.list = c(1)
params <- expand.grid(strat = strategy,
                    algo = algo_list,
                    n.cova = n.cova.list,
                    ttm.prop = ttm.prop.list,
                    event.prop = event.prop.list,
                    r.penalty = r.penalty.list,
                    n.pair = n.pair.list,
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
params$r.penalty[params$strat=='matching'] = NA
params = params[!(params$n.cova>2 & params$n.pair==200), ]

params = unique(params)

sim <- function(n = 1e6, a = 2, med0 = 200, hr.ttm = 1, hr.c = 1.5, 
                px.str = 3, scale.cens = 100, ttm.prop = 0.1, n.cova = 20) {  
    px <- rbinom(n, 1, 0.5)

    cova <- matrix(0, ncol = n.cova, nrow = n)
    name <- NA
    ####################### test binary ##########################################
    for (k in 1:n.cova) { 
        cova[, k] <- rbinom(n, 1, 0.3)
        #cova[, k] <- runif(n, -0.5, 0.5)
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
scanning.scale = function(event.prop.list) {
    res = c()
    for (n.cova in n.cova.list) {
        for (ttm.prop in ttm.prop.list) {
            for (event.prop in event.prop.list) {
                set.seed(1)
                h = 5e-3
                scale.cens = 50
                data = sim(n = 1e6, 
                            hr.ttm = hr.ttm, 
                            px.str = px.str, 
                            ttm.prop = ttm.prop, 
                            scale.cens = scale.cens, 
                            n.cova = n.cova)
                check = event.prop - mean(data$stt)
                n = 1
                if (abs(check)<h) {new.scale = scale.cens}
                while (abs(check)>h) {
                    prev.check = check
                    a = ifelse(check <0, -1, 1)
                    new.scale = scale.cens + scale.cens*n*0.1*a
                    cat(event.prop,new.scale,n,"\r")
                    flush.console()
                    set.seed(1)
                    data = sim(n = 1e6, 
                            hr.ttm = hr.ttm, 
                            px.str = px.str, 
                            ttm.prop = ttm.prop, 
                            scale.cens = new.scale, 
                            n.cova = n.cova)
                    check = event.prop - mean(data$stt)
                    if (check*prev.check <0) {a = a*(-1); n = n - exp(check)} else {n = n + exp(check)}
                    }
                set.seed(1)
                data = sim(n = 1e6, 
                          hr.ttm = hr.ttm, 
                          px.str = px.str, 
                          ttm.prop = ttm.prop, 
                          scale.cens = new.scale, 
                          n.cova = n.cova)
                res = rbind(res, c(n.cova,ttm.prop,event.prop,mean(data$stt),new.scale))            
            }              
        }
    }
    res = data.frame(res)
    colnames(res) = c('n.cova','ttm.prop','event.prop','pop.event.prop','scale.cens')
    return(res)
}

addrevcaliper<-function(dist,z,dx,rg, stdev = FALSE, penalty = max(dist$d), r.penalty = 10) {
    stopifnot(is.vector(rg)&(length(rg)==2))
    stopifnot((rg[1]<=0)&(rg[2]>=0))  
    stopifnot(is.vector(z))
    stopifnot(all((z==1)|(z==0)))
    penalty = r.penalty*penalty 
    dx1<-dx[z==1]
    dx0<-dx[z==0]  
    v <- stats::var(dx)
    sp <- sqrt(v)
    stopifnot(sp>0)
    if (stdev) rg <- rg *sp  
    m=sum(z)
    dif<-dx1[dist$start]-dx0[dist$end-m]
    d <- dist$d + as.numeric((dif>rg[1])&(dif<rg[2]))*penalty
    list(d=d,start=dist$start,end=dist$end)
}
scale.res = scanning.scale(event.prop.list = event.prop.list)
nworkers <- detectCores()
cl <- makeCluster(nworkers)
tasks <- c(1:nworkers)

worker <- function(task) {
    env <- new.env(parent = globalenv())
    env$params <- params
    env$nsim <- nsim
    env$n <- n
    env$px.str <- px.str
    env$hr.ttm <- hr.ttm
    env$scale.res <- scale.res
    env$sim <- sim
    env$addrevcaliper <- addrevcaliper
    env$nworkers <- nworkers
    env$iphase <- task
    #env$ex_res <- c()
    env$mhl_res <- c()
    env$drs_res <- c()
    env$ps_res <- c()

    for (script in c("clog_mhl1.R","clog_drs1.R", "clog_ps1.R")) { #"clog_ex.R",
        source(paste0("/work/ttkle/matching/",script), local = env) #"/work/ttkle/matching/",
    }
    list(#ex = env$ex_res,
        mhl = env$mhl_res,
        ps  = env$ps_res,
        drs = env$drs_res)
}

clusterExport(cl, c("params", "nsim", "n", "px.str", "hr.ttm",
                    "sim", "addrevcaliper", "scale.res", "tasks", "worker", "nworkers"))

clusterEvalQ(cl, {
    library(remotes)
    library(DiPs)
    library(survival)
    library(SensitivityCaseControl)
    library(stats)
    NULL
})
results <- parLapplyLB(cl, tasks, worker)
stopCluster(cl)

nms <- unique(unlist(lapply(results, names)))
res <- do.call(rbind, lapply(nms, function(nm) {do.call(rbind, lapply(results, `[[`, nm))}))
write.csv(res, paste0("/work/ttkle/matching/res_pw_checkvar.csv")) 