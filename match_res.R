# ----------------
# SIMULATION nsim
# ----------------

library("remotes")
library("DiPs")
library("survival")
library("SensitivityCaseControl")

#--------------------
# Simulation study parameters
#--------------------
strategy = c("matching", "counter-matching")

pair = c()
pair.cm = c()
cm.fail = c()

# Simulation parameters
nphase = c(1:10)
nsim = 1
n.pair = 200
n = 20000
px.str = 3
hr.ttm = 1.5

# Treatment prob = 10% and 50%
ttm.prop.list = c(0.1, 0.5)

# Scale cens --> event prob 1%, 5%, 10%
event.prop.list = c(0.01, 0.05, 0.10)

# Nb of cova
n.cova.list = c(4, 10, 20)

n.hx = n
conv = 0

# List of algorithms
algo_list = c('A1','A2')

#--------------------
# Simulation function
#--------------------

sim <- function(n = 1e6, a = 2, med0 = 200, hr.ttm = 1, hr.c = 1.5, 
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
    return(res)}
# -------------------------
# Counter-matching function
# -------------------------

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

# -------------------
# Combine all results
# -------------------
names <- c("clog_mhl1", "clog_ps1", "clog_drs1")

#d <- lapply(names, function(x) {
#  source(paste0(x, ".r"), local = TRUE)
#  return(res) 
#})
scale.res = scanning.scale(event.prop.list = event.prop.list)

res = c()
mhl_res <- c()
drs_res <- c()
ps_res <- c()
for (iphase in nphase) {
  for (x in names) {
    source(paste0(x, ".r"), local = TRUE)
  }
  res <- rbind(mhl_res, ps_res, drs_res)
  write.csv(res, paste0("match_resA1_pwtest.csv"))
}

#names(res) <- names
#res_match <- do.call(rbind, res)
#
#write.csv(res, "match_res1_test.csv")