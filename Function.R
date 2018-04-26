Q1.similation <- function(Time,N){
  
  alpha <- c(0.0100,0.0025,0.0100)
  gamma <- c(0.0025,0.0100,0.0025)
  sigma_e <- c(0.005,0.005,0.001)
  sigma_v <- c(0.001,0.001,0.005)
  y <- z <- matrix(0,Time,1)
  bhat <- b.se <- t_ratio <- p_value <- matrix(NA,N,1)
  bhat.mean <- b.se_mean <- t_ratio_mean <-reject_freq <- matrix(NA,1,3)
  
  for(i in 1:3){
    for (j in 1:N) {
      
      for(t in 1:(Time-1)){
        y[t+1] <- alpha[i] + y[t] + rnorm(1,0,sigma_e[i])
        z[t+1] <- gamma[i] + z[t] + rnorm(1,0,sigma_v[i])
      }
      
      data1 <- data.frame(cbind(c(1:Time),y,z))
      colnames(data1) <- c("Time", "y", "z")
      
      reg <- summary(lm(y~z+Time, data = data1))
      
      bhat[j] <- reg$coefficients[2,1]
      b.se[j] <- reg$coefficients[2,2]
      t_ratio[j] <- reg$coefficients[2,3]
      p_value[j] <- reg$coefficients[2,4]
    }
    
    bhat.mean[i] <- abs(mean(bhat))
    b.se_mean[i] <- mean(b.se)
    t_ratio_mean[i] <- mean(t_ratio)
    
    reject_freq[i] <- length(p_value[which(p_value<0.05)])/N
    
    cat("Completed experiment", i,"\n")
  }
  
  return(list(Time=Time, bhat.mean=bhat.mean,b.se_mean=b.se_mean, t_ratio_mean=t_ratio_mean, reject_freq=reject_freq))
}

