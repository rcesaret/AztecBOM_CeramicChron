Rich_Sim_Curv <- function(x, nsim=1000, by.d=1) { # create script and set default values for input
  
  maxn <- 2133 # the number of sample sizes considered is set as the maximum observed sample size + 5%
  step.d <- seq(1,2133,by=by.d) # this step simply creates a sequesnce from 1 to maxn by the interval set using the by.d variable
  prob.d <- colSums(x) # this step produces the pool from which random samples will be drawn by simply summing the values across all units by type (column)
  divlist <- list() # create an empty list for output
  
  for (i in step.d) {divlist[[length(divlist)+1]] <- specnumber(rmultinom(nsim,i,prob.d),MARGIN=2)} 
  mean.d <- rapply(divlist,mean) # this rapply command recursively calculates the mean richness for each sample size in divlist
  sd.d <- rapply(divlist,sd) # this rapply command recursively calculates the standard deviation of richness for each sample size in divlist
  
  nm <- deparse(substitute(x))
  
  CI <- data.frame(step.d,mean.d,sd.d)
  colnames(CI) <- c("step.d","mean.d","sd.d")
  CI$upper0.66 <- CI$mean.d+(CI$sd.d*qnorm((1-0.66)/2))
  CI$lower0.66 <- CI$mean.d-(CI$sd.d*qnorm((1-0.66)/2))
  CI$upper0.95 <- CI$mean.d+(CI$sd.d*qnorm((1-0.95)/2))
  CI$lower0.95 <- CI$mean.d-(CI$sd.d*qnorm((1-0.95)/2))
  
  assign(paste0("RichSimCurv.",nm,".RarefCurv"), CI, envir = .GlobalEnv)
  
}