require(vegan) # this script requires the vegan package

FishDiv_Sim <- function(x, nsim=1000, by.d=1) { # create script and set default values for input
  
  maxn <- max(rowSums(x))*1.05 # the number of sample sizes considered is set as the maximum observed sample size + 5%
  step.d <- seq(1,maxn,by=by.d) # this step simply creates a sequesnce from 1 to maxn by the interval set using the by.d variable
  prob.d <- colSums(x) # this step produces the pool from which random samples will be drawn by simply summing the values across all units by type (column)
  divlist <- list() # create an empty list for output
  
  for (i in step.d) {divlist[[length(divlist)+1]] <- microbiome::alpha((rmultinom(nsim,i,prob.d)), index = "diversity_fisher")} 
  mean.d <- rapply(divlist,mean) # this rapply command recursively calculates the mean richness for each sample size in divlist
  sd.d <- rapply(divlist,sd) # this rapply command recursively calculates the standard deviation of richness for each sample size in divlist
  Resid <- data.frame(matrix(NA, nrow = nrow(x), ncol = 1))
  rownames(Resid) <- rownames(x)
  colnames(Resid) <- c("n")
  Resid$n <- rowSums(x)
  xxx <- as.data.frame(t(as.matrix(x)))
  Resid$FishDiv.obs <- microbiome::alpha(xxx, index = "diversity_fisher")
  Resid$FishDiv.sim.mean <- NA
  Resid$FishDiv.sim.sd <- NA
  
  for (j in 1:nrow(x)) { 
    Resid$FishDiv.sim.mean[j] <- ifelse(Resid[j,1] == 0, NA, mean.d[Resid[j,1]]) # Resid$n[j]
    Resid$FishDiv.sim.sd[j] <- ifelse(Resid[j,1] == 0, NA, sd.d[Resid[j,1]]) # Resid$n[j]
  }
  Resid$FishDiv.resid <- ((Resid$FishDiv.obs) - (Resid$FishDiv.sim.mean))
  Resid$FishDiv.resid.sim.z <- ((Resid$FishDiv.obs) - (Resid$FishDiv.sim.mean))/Resid$FishDiv.sim.sd
  Resid$FishDiv.resid.sim.p <- pnorm(Resid$FishDiv.obs, mean = Resid$FishDiv.sim.mean, sd = Resid$FishDiv.sim.sd, lower.tail = TRUE)
  meanR = mean(Resid$FishDiv.obs, na.rm=T)
  sdR = sd(Resid$FishDiv.obs, na.rm=T)
  Resid$FishDiv.resid.obs.z <- ((Resid$FishDiv.obs) - meanR)/sdR
  Resid$SurvReg <- Attr$SurvReg
  Resid$Label <- Attr$NameShort
  nm <- deparse(substitute(x))
  assign(paste0("FishDivSim.",nm,".Resid"), Resid, envir = .GlobalEnv)
  
  CI <- data.frame(step.d,mean.d,sd.d)
  colnames(CI) <- c("step.d","mean.d","sd.d")
  CI$upper0.66 <- CI$mean.d+(CI$sd.d*qnorm((1-0.66)/2))
  CI$lower0.66 <- CI$mean.d-(CI$sd.d*qnorm((1-0.66)/2))
  CI$upper0.95 <- CI$mean.d+(CI$sd.d*qnorm((1-0.95)/2))
  CI$lower0.95 <- CI$mean.d-(CI$sd.d*qnorm((1-0.95)/2))
  
  assign(paste0("FishDiv.",nm,".RarefCurv"), CI, envir = .GlobalEnv)
  return(Resid)
  #return(CI)
  
}