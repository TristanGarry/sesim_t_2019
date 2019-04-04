### SETTING UP METACOMMUNITY ###

### Connectivity

# number of patches
numCom <- 5
# absent
d.a <- matrix(c(0,0,0,0,0), nrow=5, ncol=5)
diag(d.a) <- 0
# linear
d.l <- matrix(c(0,1.0,0,0,0,1.0), nrow=5, ncol=5)
# circular
d.c <- d.l
d.c[5,1] <- d.c[1,5] <- 1.0
# global
d.g <- matrix(c(1.0,1.0,1.0,1.0,1.0), nrow=5, ncol=5)
diag(d.g) <- 0
# connectivity treatment
connec <- c("d.a", "d.l", "d.c", "d.g")

# dispersal rate
rate <- 0.0001

# species-specific dispersal ability
disp <- c(0.11, 0.11, 0.11, 0.64, 0.62, 0.53, 0.90, 0.62, 0.44, 0.33, 0.71, 1.00)

# dispersal kernel
k <- 1

### Patch quality

# growth rates
r.source <- c(8.4,  16.8, 20.3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
r.sink   <- c(0.84, 1.68, 2.03, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# carrying capacity
K.source <- c(5000000000,5000000000,5000000000, 321000, 1216000, 1411000, 44000, 6823000, 95000, 280000, 190000, 407000)
K.sink   <- c(500000000, 500000000, 500000000,  1600,   5400,    62000,   100,   41000,   2000,  28000,  31000,  2000)

# number of source patches
ss.prop <- c(0,1,2,3,4,5)

### Species

# food web
source("./fw-T.R")
# number of species
nSp <- nrow(FW)

# cell mass used to calculate productivity
mass <- c(4.3e-08, 1.96e-07, 3.76e-06, 1.52e-08, 4.77e-09, 6.9e-08, 9.68e-08, 8.05e-08, 2.27e-07)

# number of environmental variables
nEnv <- 1
# environmental niche amplitude
eA <- 1
# strength of environmental effect
eff.source <- 10 ; eff.sink <- 1

# environmental fluctuation period
eP <- 0



### Simulation parameters

# number of replicates
reps <- 10

# simulation length
Tmax <- 30

# sampling frequency
sampfreq <- 1
sampleV <- seq(0, Tmax, by=sampfreq)



### ANALYSIS PARAMETERS ### 
species_w_bacteria = c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12')
species = c('V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12')
bacteria = c('V1', 'V2', 'V3')
non_species = c('time', 'dispersal', 'kernel', 'replicate', 'patch', 'connectivity', 'quality')



