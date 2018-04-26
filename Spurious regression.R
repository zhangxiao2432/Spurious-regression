rm(list = ls())

library("reshape")
library("car")
library("grid")
library("gridExtra")
library("xtable")
library("plm")
library("ggplot2")

source("Function.R")

set.seed(1030)

N <- 10000
Experiment <- seq(10,80,10)

table1 <- c()

for(i in 1:length(Experiment)){

  N0 <- 1
  #N0 <- 2
  #N0 <- Experiment[i]
  T0 <- Experiment[i]
  
  alpha <- matrix(NA,N,1)
  beta <- beta.ESE <- beta.t_value <- beta <- beta.p_value <- matrix(NA,N,1) 
  R_square <- DW <- matrix(NA,N,1)
  
  # Simulation
  for(k in 1:N){
    
    # DGP
    u <- matrix(rnorm(N0*(T0+1000)),N0,T0+1000)
    e <- matrix(rnorm(N0*(T0+1000),0,2),N0,T0+1000)
    x <- y <- matrix(0,N0,T0+1000)
    
    
    for(i in 1:N0){
      for(t in 1:(T0+999)){
        
        x[i,t+1] <- x[i,t] + e[i,t+1]
        y[i,t+1] <- y[i,t] + u[i,t+1]
        
        
      }
    }
    
    x <- x[,1001:(1000+T0)]
    y <- y[,1001:(1000+T0)]
    
    if(N0==1){
      
      panel <- as.data.frame(cbind(c(1:T0),y,x))
      colnames(panel) <- c("T","Y","X")
      rownames(panel) <- c()
      
      reg <- summary(lm(Y ~ X, data=panel))
      
      beta[k] <- reg$coefficients[2,1]
      beta.ESE[k] <- reg$coefficients[2,2]
      beta.t_value[k] <- reg$coefficients[2,3]
      beta.p_value[k] <- reg$coefficients[2,4]
      
      R_square[k] <- reg$r.squared
      DW.temp <- durbinWatsonTest(lm(Y ~ X, data=panel))
      DW[k] <- DW.temp$dw
      
    } else {
      
      # x
      x <- as.data.frame(cbind(c(1:nrow(x)),x)) 
      x.long <- reshape(x, varying = list(names(x)[2:ncol(x)]),
                        idvar = "V1", timevar = "T", times = c(1:T0), direction = "long")
      
      # y
      y <- as.data.frame(cbind(c(1:nrow(y)),y)) 
      y.long <- reshape(y, varying = list(names(y)[2:ncol(y)]),
                        idvar = "V1", timevar = "T", times = c(1:T0), direction = "long")
      
      panel <- cbind(y.long,x.long)
      panel <- panel[,c(1,2,3,6)]
      colnames(panel) <- c("N","T","Y","X")
      rownames(panel) <- c()
      
      reg <- summary(plm(Y ~ X, data=panel))
      
      beta[k] <- reg$coefficients[1]
      beta.ESE[k] <- reg$coefficients[2]
      beta.t_value[k] <- reg$coefficients[3]
      beta.p_value[k] <- reg$coefficients[4]
      
      R_square[k] <- reg$r.squared[1]
      DW.temp <- pdwtest(reg)
      DW[k] <- DW.temp$statistic
    }
 
    remove("u","e")
    
    if(k %% 5000 == 0){
      cat("Completed: ", k,"/10000 \n")
    }
    
  }
  
  # mean
  alpha.mean <- mean(alpha)
  beta.mean <- mean(beta)
  beta.ESE.mean <- mean(beta.ESE)
  beta.t_value.mean <- mean(beta.t_value)
  beta.p_value.mean <- mean(beta.p_value)
  R_square.mean <- mean(R_square)
  DW.mean <- mean(DW)
  beta.SSD <- sd(beta)
  beta.t_value.SSD <- sd(beta.t_value)
  
  # Table1
  table1.temp <- c(N0,beta.mean,beta.SSD,beta.ESE.mean,beta.t_value.mean,
                              beta.t_value.SSD,beta.p_value.mean,R_square.mean,DW.mean)
  
  table1 <- rbind(table1, table1.temp)
  
  cat("Completed T0 = ", T0," \n")
}

colnames(table1) <- c("N(T)","beta","SSD","ESE","t_beta","SSD_t","P_beta","R_square","DW")
row.names(table1) <- c()

# Create PDF file
pdf("Table 1.pdf", width = 12, height = 12, onefile = TRUE)
grid.table(format(round(table1,digits = 4), nsmall = 4))
dev.off()   
  
# Create Latex code
xtable(table1)

#================================================================================
#Q1
Biases <- Est.std.err <- t_ratio <- reject_freq <- dt <- c()

for (t in 40:100) {
  
  Time <- t
  cat("Time=", t,"/100 \n")
  data.temp <- Q1.similation(Time,1000)
  
  temp1 <- matrix(data.temp$bhat.mean,1,3)
  temp2 <- matrix(data.temp$b.se_mean,1,3)
  temp3 <- matrix(data.temp$t_ratio_mean,1,3)
  temp4 <- matrix(data.temp$reject_freq,1,3)
  
  dt <- rbind(dt,Time)
  Biases <- rbind(Biases, temp1)
  Est.std.err <- rbind(Est.std.err, temp2)
  t_ratio <- rbind(t_ratio, temp3)
  reject_freq <- rbind(reject_freq, temp4)
  
  remove(temp1,temp2,temp3,temp4)
  
}

h1 <- ggplot(NULL,aes(x=dt)) + ggtitle("Biases")+
  geom_line(aes(y=Biases[,1], col="b1"),size=1,linetype="solid")+ 
  geom_line(aes(y=Biases[,2], col="b2"),size=1,linetype="longdash")+ 
  geom_line(aes(y=Biases[,3], col="b3"),size=1,linetype="twodash")+
  labs(x="", y="")+
  scale_color_manual(values=c("b1"="#00AFBB", "b2"="#f8766d","b3"="#E7B800"))+
  theme_bw()+theme(legend.position = c(0.8,0.8),plot.title = element_text(hjust=0.5))


h2 <- ggplot(NULL,aes(x=dt)) + ggtitle("Estimated standard errors")+
  geom_line(aes(y=Est.std.err[,1], col="b1"),size=1,linetype="solid")+ 
  geom_line(aes(y=Est.std.err[,2], col="b2"),size=1,linetype="longdash")+ 
  geom_line(aes(y=Est.std.err[,3], col="b3"),size=1,linetype="twodash")+
  labs(x="", y="")+
  scale_color_manual(values=c("b1"="#00AFBB", "b2"="#f8766d","b3"="#E7B800"))+
  theme_bw()+theme(legend.position = c(0.8,0.4),plot.title = element_text(hjust=0.5))


h3 <- ggplot(NULL,aes(x=dt)) + ggtitle("Mean t-ratios")+
  geom_line(aes(y=t_ratio[,1], col="b1"),size=1,linetype="solid")+ 
  geom_line(aes(y=t_ratio[,2], col="b2"),size=1,linetype="longdash")+ 
  geom_line(aes(y=t_ratio[,3], col="b3"),size=1,linetype="twodash")+
  labs(x="", y="")+
  scale_color_manual(values=c("b1"="#00AFBB", "b2"="#f8766d","b3"="#E7B800"))+
  theme_bw()+theme(legend.position = c(0.8,0.8),plot.title = element_text(hjust=0.5))


h4 <- ggplot(NULL,aes(x=dt)) + ggtitle("t rejection frequencies")+
  geom_line(aes(y=reject_freq[,1], col="b1"),size=1,linetype="solid")+ 
  geom_line(aes(y=reject_freq[,2], col="b2"),size=1,linetype="longdash")+ 
  geom_line(aes(y=reject_freq[,3], col="b3"),size=1,linetype="twodash")+
  labs(x="", y="")+
  scale_color_manual(values=c("b1"="#00AFBB", "b2"="#f8766d","b3"="#E7B800"))+
  theme_bw()+theme(legend.position = c(0.2,0.8),plot.title = element_text(hjust=0.5))


pdf("Figure 1.pdf", width=12, height=8)
grid.arrange(h1,h2,h3,h4, ncol=2)
dev.off()



