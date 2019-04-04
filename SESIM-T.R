###############################################
#           SESIM - microcosm model           #
# for patch quality x connectivity experiment #
###############################################

library(vegan)
library(igraph)
library(reshape)
library(ggplot2)

source('./sesim_parameters.R')

### SET UP DATA STORAGE ###

# empty array for saving sampled data
X_save <- array(data=NA, dim=c(numCom, nSp, Tmax/sampfreq))

# experiment factorial design
experiment <- data.frame(connectivity = rep(NA, each=length(connec)*length(ss.prop)), 
                         patchquality = rep(NA))

# data frames for results
results <- data.frame(treatment = rep(1:(dim(experiment)[1]), each=numCom*reps),
                      r = rep(1:reps, each=numCom),
                      patch = rep(1:numCom))

# Model

treatment_id = 1
loop_id = 1

for(conn in seq_along(connec)){
  for(ss in ss.prop){
    for(r in 1:reps){
        
      cc <- get(connec[conn])
      sampling = 1
        
      # growth rate
      growth <- C <- matrix(c(rep(r.source, times=ss), rep(r.sink, times=(5-ss))), nrow=nSp, ncol=5)
      # carrying capacity
      K <- matrix(c(rep(K.source, times=ss), rep(K.sink, times=(5-ss))), nrow=nSp, ncol=5)
        
      # matrix of species interaction
      BB  <- data.matrix(as.data.frame(FW))
      matriz  <- novamatriz <- BB
      for(j in 1:ncol(matriz)){
        for(i in 1:nrow(matriz)){
          novamatriz[i,j] <- runif(1, min=min(0, matriz[i,j]), max=max(0, matriz[i,j]))
        }
      }
      BB <- novamatriz/100
        
      # dispersal ability
      dispersal_m <- exp(-k*cc) - diag(nrow(cc))
      dispersal_m[dispersal_m == 1] <- 0
      nlinks <- rowSums(dispersal_m != 0)
      maxlinks <- 4
      
      # species environmental optimum
      Env_Opt <- matrix(runif(nSp, 0, 1), nSp, nEnv)
      #Env_Opt <- matrix(1, nSp, nEnv)
      # environmental condition in each patch per species
      Env <- matrix(rep(1), nrow=nSp, ncol=numCom)
      # strength of environment effect
      enveff <- c(rep(eff.source, times=ss), rep(eff.sink, times=(5-ss)))
      # environmental match
      enviro <- t(apply(1-abs(Env_Opt[,1] - Env), 1, function(x) x/enveff))

      # species initial abundances
      source <- c(25000000,25000000,25000000,100,100,100,100,100,100,100,100,100)
      sink   <- c(2500000, 2500000, 2500000, 100,100,100,100,100,100,100,100,100)
      X <- matrix(c(rep(source, times=ss), rep(sink, times=(5-ss))), nrow=nSp, ncol=5)
      Xd <- t(X)
        
      for(l in 1:(Tmax)){
        
        # growth
        growth <- C
        
        # species interactions
        interactions <- BB %*% (X*enviro)
        
        # migrants
        migrants <- matrix(0, nSp, numCom)
        for(i in 1:numCom){
          migrants[,i] <- (X[,i]*(rate*disp))*(1/(maxlinks/(length(dispersal_m[,i][dispersal_m[,i]>0]))))
        }
        
        immigrants <- matrix(0, nSp, numCom)
        for(i in 1:nSp){
          for(j in 1:numCom){
            immigrants[i,j] <- sum(migrants[i,][dispersal_m[j,]!=0]/nlinks[dispersal_m[j,]!=0])
          }
        }
        
        # lotka-volterra equation
        Xt <- X + X * (growth + interactions) + immigrants - migrants
        
        # extinctions
        Xt[Xt < 0] <- 0 ; Xt[is.na(Xt)] <- 0 ; Xt[is.infinite(Xt)] <- 0
        X <- Xt
        Xd <- t(X)
        
        # sample
        if(l==sampleV[l/sampfreq+1] && l<=Tmax){
          X_save[,,sampling] <- Xd
          sampling <- sampling+1
        }
        
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
      results$div_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="shannon")
      results$div_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="shannon")
      # beta diversity
      results$div_b [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="shannon")-mean(vegan::diversity(Xd[,-1:-3], index="simpson"))
      # richness
      results$rich_l [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber(Xd[,-1:-3])
      results$rich_r [(((loop_id-1)*5)+1):(loop_id*5)] <- specnumber((apply(Xd[,-1:-3], 2, sum)))
      # evenness
      results$eve_l [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity(Xd[,-1:-3], index="shannon")/log(specnumber(Xd[,-1:-3])) 
      results$eve_r [(((loop_id-1)*5)+1):(loop_id*5)] <- vegan::diversity((apply(Xd[,-1:-3], 2, sum)), index="shannon")/log(specnumber((apply(Xd[,-1:-3], 2, sum)))) 
      # productivity
      productivity <- log(Xd[,-1:-3]+10) * mass
      results$prod_l [(((loop_id-1)*5)+1):(loop_id*5)] <- apply(productivity, 1, sum)
      results$prod_r [(((loop_id-1)*5)+1):(loop_id*5)] <- sum(productivity)
      
      loop_id <- loop_id + 1
      
    } # finished all replicates for one treatment level
      
    # regional species abundance time series
    visual <- melt(X_save)
    #mainplotitle <- paste("abund_replicate", r, treatment_id, ".png", sep="")	
    #png(filename=mainplotitle, width=29, height=21, units="cm", res=300)
    #print(
    #  ggplot(visual, aes(x = Var3, y = log(value), colour=factor(Var2))) +
    #  stat_summary(fun.y=sum, na.rm=T, geom="line", size=1.2) + 
    #  labs(x = "Time", y = "Abundance") +
    #  facet_wrap(~ Var1) +
    #  theme_bw()
    #  )
    #dev.off()
    
    id <- treatment_id
    treatment_id <- id+1
    
  } 
} # end

write.csv(results, file="cp-results.csv")

names_samples <- grep("rep_sampled_", x=ls(), value=T)
rep_samples <- do.call(rbind, mget(names_samples))
str(rep_samples)
write.csv(rep_samples, file="cp-timeseries.csv")
rm(list=names_samples)

# # assemble results
# results.data <- data.frame(connectivity=results$connectivity,
#                            patchquality=results$patchquality,
#                            scale=rep(c("local","regional"), each=(dim(results)[1])),
#                            metric=rep(c("1diversity", "2richness", "3evenness", "4productivity"), each=(dim(results)[1])*2),
#                            value=c(results$div_l, results$div_r, results$rich_l, results$rich_r,
#                                    results$eve_l, results$eve_r, results$prod_l, results$prod_r))
# 
# # results plot
# ggplot(results.data[results.data$scale == "local",], aes(x=(patchquality), y=value)) +
#   labs(x="Patch quality", y="Metrics") +
#   geom_ribbon(stat="summary", fun.ymin=function(x) {mean(x)-(sd(x)/sqrt(length(x)))},
#               fun.ymax=function(x) {mean(x)+(sd(x)/sqrt(length(x)))}, alpha=0.3, aes(fill=connectivity)) +
#   stat_summary(fun.y=mean, na.rm=T, geom="line", size=0.5, show.legend=T, aes(colour=connectivity)) +
#   stat_summary(fun.y=mean, na.rm=T, geom="point", size=3, show.legend=T, aes(colour=connectivity)) +
#   #scale_fill_manual(values=colors) +
#   #scale_colour_manual(values=colors) +
#   facet_wrap(~metric, nrow=2, scales="free", labeller=as_labeller(c("1diversity"="a)","2richness"="b)","3evenness"="c)","4productivity"="d)"))) +
#   theme_gray() + theme(panel.background=element_rect(fill="white"), panel.border=element_rect(linetype="solid", fill=NA),
#         axis.ticks=element_line(size=0.3), axis.ticks.length=unit(0.2, "cm"),
#         axis.title.y=element_text(size=rel(1.2)), axis.title.x=element_text(size=rel(1.2)),
#         axis.text.y=element_text(size=rel(1.1)), axis.text.x=element_text(size=rel(1.1)),
#         strip.text=element_text(size=rel(1)), strip.background=element_rect(fill="white"),
#         strip.text.x=element_text(angle = 0, hjust = 0))

