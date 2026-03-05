#--------------------
# Simulation study parameters
#--------------------
drs_res <- c()
pair = c()
pair.cm = c()
cm.fail = c()

# Simulation parameters
nsim = 1
n.pair = 200
n = 20000
px.str = 3
hr.ttm = 1

# Treatment prob = 10% and 50%
ttm.prop.list = c(0.1, 0.5)

# Scale cens --> event prob 1%, 5%, 10%
scale.cens.list = c(26.3, 55.2, 79.6)

# Nb of cova
n.cova.list = c(4, 10, 20)
n.hx = n
conv = 0

# List of algorithms
algo_list = c('A1', 'A2')