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
mhl_res <- c()
drs_res <- c()
ps_res <- c()

pair = c()
pair.cm = c()
cm.fail = c()

# Simulation parameters
nsim = 100
n.pair = 200
px.str = 3
hr.ttm = 1

# Treatment prob = 10% and 50%
ttm.prop.list = c(0.1, 0.5)

# Scale cens --> event prob 1%, 5%, 10%
scale.cens.list = c(26.3, 55.2, 79.6)

# Nb of cova
n.cova.list = c(4, 10, 20)
conv = 0

# List of algorithms
algo_list = c('A1', 'A2')

# -------------------
# Combine all results
# -------------------
names <- c("clog_mhl1", "clog_ps1", "clog_drs1")
#directory <- "C:/R/matchsim1/Code/matching/"

d <- lapply(names, function(x) {
  source(paste0(x, ".r"), local = TRUE)
  return(res) 
})

names(d) <- names
res_match <- do.call(rbind, d)

write.csv(res_match, "match_res1.csv")
