  require(vegan) # this script requires the vegan package
ShanDiversity_Sim <- function(x, nsim=1000, by.d=1) { # create script and set default values for input
  maxn <- max(rowSums(x))*1.05 # the number of sample sizes considered is set as the maximum observed sample size + 5%
  step.d <- seq(1,maxn,by=by.d) # this step simply creates a sequesnce from 1 to maxn by the interval set using the by.d variable
  
  prob.d <- colSums(x) # this step produces the pool from which random samples will be drawn by simply summing the values across all units by type (column)
  divlist <- list() # create an empty list for output
  
  for (i in step.d) {divlist[[length(divlist)+1]] <- vegan::diversity(rmultinom(nsim,i,prob.d), MARGIN=2)} 
  mean.d <- rapply(divlist,mean) # this rapply command recursively calculates the mean richness for each sample size in divlist
  sd.d <- rapply(divlist,sd) # this rapply command recursively calculates the standard deviation of richness for each sample size in divlist
  Resid <- data.frame(matrix(NA, nrow = nrow(x), ncol = 1))
  rownames(Resid) <- rownames(x)
  colnames(Resid) <- c("n")
  Resid$n <- rowSums(x)
  Resid$ShanDiv.obs <- vegan::diversity(x,MARGIN=1)
  Resid$ShanDiv.sim.mean <- NA
  Resid$ShanDiv.sim.sd <- NA
  
  for (j in 1:nrow(x)) { 
    Resid$ShanDiv.sim.mean[j] <- mean.d[Resid[j,1]] # Resid$n[j]
    Resid$ShanDiv.sim.sd[j] <- sd.d[Resid[j,1]] # Resid$n[j]
  }
  Resid$ShanDiv.resid <- ((Resid$ShanDiv.obs) - (Resid$ShanDiv.sim.mean))
  Resid$ShanDiv.resid.sim.z <- ((Resid$ShanDiv.obs) - (Resid$ShanDiv.sim.mean))/Resid$ShanDiv.sim.sd
  Resid$ShanDiv.resid.sim.p <- pnorm(Resid$ShanDiv.obs, mean = Resid$ShanDiv.sim.mean, sd = Resid$ShanDiv.sim.sd, lower.tail = TRUE)
  Resid$ShanDiv.resid.obs.z <- ((Resid$ShanDiv.obs) - (mean(Resid$ShanDiv.obs, na.rm=T)))/sd(Resid$ShanDiv.obs, na.rm=T)
  nm <- deparse(substitute(x))
  
  assign(paste0("ShanDivSim.",nm,".Resid"), Resid, envir = .GlobalEnv)
  
  CI <- data.frame(step.d,mean.d,sd.d)
  colnames(CI) <- c("step.d","mean.d","sd.d")
  CI$upper0.66 <- CI$mean.d+(CI$sd.d*qnorm((1-0.66)/2))
  CI$lower0.66 <- CI$mean.d-(CI$sd.d*qnorm((1-0.66)/2))
  CI$upper0.95 <- CI$mean.d+(CI$sd.d*qnorm((1-0.95)/2))
  CI$lower0.95 <- CI$mean.d-(CI$sd.d*qnorm((1-0.95)/2))
  
  assign(paste0("ShanDivSim.",nm,".RarefCurv"), CI, envir = .GlobalEnv)
  return(Resid)
  #return(CI)
  
}