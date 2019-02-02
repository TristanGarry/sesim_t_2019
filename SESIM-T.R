###############################################
#           SESIM - microcosm model           #
# for patch quality x connectivity experiment #
###############################################

library(vegan)
library(igraph)
library(reshape)
library(ggplot2)


### Setting up parameters
source("./sesim_parameters.R")

### Species
# food web
source("./fw.R")
# number of species
nSp <- nrow(FW)

### Simulation parameters

# number of replicates
reps <- 10
xseed <- as.integer(runif(reps)*100000)

# simulation length
Tmax <- 30

# sampling frequency
sampfreq <- 1
sampleV <- seq(0, Tmax, by=sampfreq)
# empty array for saving sampled data
X_save <- array(data=NA, dim=c(numCom, nSp, Tmax/sampfreq))

# data frames for results
results <- data.frame(treatment = rep(1:(dim(experiment)[1]), each=numCom*reps),
                      r = rep(1:reps, each=numCom),
                      patch = rep(1:numCom))

### Experimental design

# experiment factorial design
experiment <- data.frame(connectivity = rep(NA, each=length(connec)*length(ss.prop)), 
                         patchquality = rep(NA))


# Model

treatment_id = 1
loop_id = 1

for(conn in seq_along(connec)){
  for(ss in ss.prop){
    repeat{
      for(r in 1:reps){
        
        set.seed(xseed[r])
        
        cc <- get(connec[conn])
        sampling = 1
        
        # growth rate
        C <- matrix(c(rep(r.source, times=ss), rep(r.sink, times=(5-ss))), nrow=nSp, ncol=5)
        # carrying capacity
        K <- matrix(c(rep(K.source, times=ss), rep(K.sink, times=(5-ss))), nrow=nSp, ncol=5)
        
        # matrix of species interaction
        BB  <- data.matrix(as.data.frame(FW))
        matriz  <- BB
        novamatriz <- matriz
        for(j in 1:ncol(matriz)){
          for(i in 1:nrow(matriz)){
            novamatriz[i,j] <- runif(1, min=min(0, matriz[i,j]), max=max(0, matriz[i,j]))
          }
        }
        BB <- novamatriz/100
        
        # dispersal ability
        d_exp <- exp(-k*cc) - diag(nrow(cc))
        dispersal_m <- apply(d_exp, 1, function(x) x/sum(x))
        dispersal_m[is.na(dispersal_m)] <- 0
        
        # species environmental optimum
        Env_Opt <- matrix(1, nSp, nEnv)
        # environmental condition in each patch per species
        Env <- matrix(rep(1), nrow=nSp, ncol=numCom)
        # strength of environment effect
        enveff <- c(rep(eff.source, times=ss), rep(eff.sink, times=(5-ss)))
        # environmental match
        enviro  <- t(apply( 1-(abs(Env-Env_Opt[,1])), 1, function(x) abs(x * enveff)))
        envPrey <- t(apply(   (abs(Env-Env_Opt[,1])), 1, function(x) abs(x * enveff)))
        envPred <- t(apply( 1-(abs(Env-Env_Opt[,1])), 1, function(x) abs(x * enveff)))
        #envPrey[envPrey == 0] <- 0.1 ; envPred[envPred == 0] <- 0.1 ; enviro[enviro == 0] <- 0.1
        
        # species initial abundances
        source <- c(25000000,25000000,25000000,20,100,20,100,100,20,20,20,100)
        sink   <- c(2500000, 2500000, 2500000, 20,100,20,100,100,20,20,20,100)
        X <- matrix(c(rep(source, times=ss), rep(sink, times=(5-ss))), nrow=nSp, ncol=5)
        Xd <- t(X)
        
        for(l in 1:(Tmax)){
          
          # update matrix of species interaction
          if(l != 1){
            BB  <- data.matrix(as.data.frame(FW))
            matriz  <- BB
            novamatriz <- matriz
            for(j in 1:ncol(matriz)){
              for(i in 1:nrow(matriz)){
                novamatriz[i,j] <- runif(1, min=min(0, matriz[i,j]), max=max(0, matriz[i,j]))
              }
            }
            BB <- novamatriz/100
          }
          
          # patch specific growth rate
          diff <- C/enviro
          diff[is.na(diff)] <- 0
          diff[is.infinite(diff)] <- 0
          diff[diff <= 0] <- 0
          # growth
          growth <- ((X*diff)*((K-X)/K)) 
          growth[is.na(growth)] <- 0
          for(patch in 1:5){
            if(all(X[1:3 ,patch] == 0)) { growth[4:10,patch] <- 0  }
            if(all(X[7:10,patch] == 0)) { growth[11  ,patch] <- 0  }
            if(all(X[7:11,patch] == 0)) { growth[12  ,patch] <- 0  }
          }
          #growth[growth < 0] <- 0
          
          ## species interactions
          # BB + and -
          BBpos <- BBneg <- BB
          BBpos[BBpos < 0] <- 0 ; BBneg[BBneg > 0] <- 0
          # prey and predator abundances 
          Xprey <- Xpredator <- X
          Xprey[12,] <- 0 ; Xpredator[1:3,] <- 0
          # interactions
          impact  <- BBneg %*% ((Xprey*(envPrey)) + (Xpredator*envPred))
          #impact[((Xprey*envPrey) & (Xpredator*envPred)) <= 0] <- 0 
          benefit <- BBpos %*% ((Xprey*(envPrey)) + (Xpredator*envPred))
          #benefit[((Xprey*envPrey) & (Xpredator*envPred)) <= 0] <- 0
          
          # full interaction
          #interactions <- growth + (BB%*%X)
          interactions <- growth + (benefit+impact)
          # migrants
          #Migrants   <- apply(X, 2, function(x) x * (rate/disp))
          mMigrants <- matrix(NA, nSp, numCom)
          for(j in 1:numCom){
            mMigrants[,j] <- X[,j]*(1/(4/(length(dispersal_m[j,][dispersal_m[j,] > 0 & dispersal_m[j,] != Inf]))))
          }
          Migrants <- (apply(mMigrants, 2, function(x) x * (rate/disp)))
          Migrants[is.na(Migrants)] <- 0
          Migrants[is.infinite(Migrants)] <- 0
          # immigrants
          Immigrants <- matrix(NA, nSp, numCom)
          for(i in 1:nSp){
            for(j in 1:numCom){
              Immigrants[i,j] <- sum(Migrants[i,]*dispersal_m[j,])
            }
          }
          Immigrants[is.na(Immigrants)] <- 0
          Immigrants[is.infinite(Immigrants)] <- 0
          
          # Lotka-Volterra model
          Xt <- X + interactions + Immigrants - Migrants
          
          Xhold <- Xt
          Xhold[(Xhold < 0)] <- 0
          Xhold[is.na(Xhold)] <- 0
          Xhold[is.infinite(Xhold)] <- 0
          X <- Xhold
          Xd <- t(X)
          if(all(is.nan(X))) break
          
          if(l==sampleV[l/sampfreq+1] && l<=Tmax){
            X_save[,,sampling] <- Xd
            sampling <- sampling + 1
          } # end sampling
          
        } # end all time steps
        
        # saving sampled data for one replicate and one treatment level
        treat <- data.frame(time=(rep(seq_len(dim(X_save)[3]),each=dim(X_save)[1]))*sampfreq, 
                            dispersal = rep(rate), 
                            kernel = rep(k), 
                            replicate = rep(r), 
                            patch = rep(1:5),
                            connectivity = rep(conn),
                            quality = rep(ss))
        X_saved <- as.data.frame(apply(X_save, 2, cbind))
        assign(paste("rep_sampled_", loop_id, "_", r, sep=""), cbind(treat, X_saved))
        
        # adding treatment level to experiment 
        experiment$connectivity [treatment_id] <- connec[conn]
        experiment$patchquality [treatment_id] <- ss
        
        # adding results
        results$connectivity [(((loop_id-1)*5)+1):(loop_id*5)] <- connec[conn]
        results$patchquality [(((loop_id-1)*5)+1):(loop_id*5)] <- ss
        # diversity
        results$div_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="simpson")
        results$div_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="simpson")
        # beta diversity
        results$div_b [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="simpson")-mean(vegan::diversity(Xd[,-1:-3], index="simpson"))
        # richness
        results$rich_l [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber(Xd[,-1:-3])
        results$rich_r [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber((apply(Xd[,-1:-3], 2, sum)))
        # evenness
        results$eve_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="simpson")/log(specnumber(Xd[,-1:-3])) 
        results$eve_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="simpson")/log(specnumber((apply(Xd[,-1:-3], 2, sum)))) 
        # productivity
        productivity <- Xd[,-1:-3] * mass
        results$prod_l [(((loop_id-1)*5)+1):(loop_id*5)] <- apply(productivity, 1, sum)
        results$prod_r [(((loop_id-1)*5)+1):(loop_id*5)] <- sum(productivity)
        
        loop_id <- loop_id + 1
        
      } # finished all replicates for one treatment level
      
      # regional species abundance time series
      visual <- melt(X_save)
      #mainplotitle <- paste("abund_replicate", r, treatment_id, ".png", sep="")	
      #png(filename=mainplotitle, width=29, height=21, units="cm", res=300)
      #print(
      #  ggplot(visual, aes(x = X3, y = log(value), colour=factor(X2))) +
      #  stat_summary(fun.y=sum, na.rm=T, geom="line", size=1.2) + 
      #  labs(x = "Time", y = "Abundance") +
      #  theme_bw()
      #  )
      #dev.off()
      
      id <- treatment_id
      treatment_id <- id+1
      
      break
      
    } # end repeat
  } 
} # end

write.csv(results, file="cp-results.csv")

names_samples <- grep("rep_sampled_", x=ls(), value=T)
rep_samples <- do.call(rbind, mget(names_samples))
str(rep_samples)
write.csv(rep_samples, file="cp-timeseries.csv")
