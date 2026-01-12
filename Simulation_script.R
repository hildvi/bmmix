################################################################################
## Script for making table 1 in Gangsei & Vinje (2026)
## Authors: Hilde Vinje and Lars Erik Gangsei

## Clean workspace
rm(list = ls())  

## Load the beta3ext and bmmix packages
devtools::install_git(url = 'https://github.com/hildvi/bmmix')
devtools::install_git(url = 'https://github.com/hildvi/beta3ext')


library(beta3ext)
library(bmmix)

## Set dimensions and parameters in model

nrep <- 1000 # Number of simulations

### Dimentions
w <- 4
p <- 3
n <- 15

### Set parameters
sigma2 <- 3
sigma2u <- 0.5

delta <- w*sigma2u/(sigma2 + w*sigma2u)
beta <- c(1,2,-0.5)#as.matrix(rnorm(p))

### Put parameters in one vector
par_vec <- c(delta,sigma2,sigma2u,as.vector(beta))
Par_mat <- matrix(par_vec,nrep,length(par_vec),byrow = TRUE)

## Set up model and environment for evaluation
bmlmerform <- formula(y ~ X + (1|id))
ID <- matrix(rep(1:n,w),n,w,byrow=FALSE)

eval <- matrix(0,0,7, dimnames = 
                 list(NULL,
              c('delta','sigma2','sigma2u','(intercept)','x1','x2','method')))

## Start simulation 
### Set seeds to replicate results.
set.seed(15)
new_seeds <- sample(1:(10*nrep),nrep)

for(i in 1:nrep)
{
  set.seed(new_seeds[i]) #New seed 
  
### Simulate new design matrices and random effects
  X <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1))
  U <- matrix(rep(rnorm(n,0,sqrt(sigma2u)),w),n,w,byrow=FALSE)
  E <- matrix(rnorm(n*w,0,sqrt(sigma2)),n,w)
  Y <- matrix(X%*%beta,n,w,byrow=FALSE) + U + E

### Put data in to data frame
  mixdata <- data.frame(y = as.vector(Y),
                        X = I(kronecker(as.matrix(rep(1,w)),X[,-1])),
                        id = as.vector(ID))
  
### Get the result from Bayesian inference/ empirical Bayes
  bmlmer_res <- bmlmer(mod_form = bmlmerform,data = mixdata,
                       nsim = 5000,beta0=NULL,nu1 = 2,nu2 = 0.01,nu3 = 0.01,
                       simreturn = FALSE,empirical_bayes = TRUE,nu_max = c(2.0001,10,10))
  
  print(bmlmer_res$eb_nu)
  eval <- rbind(eval,cbind(cbind(bmlmer_res$Quantiles[c(1,3),],1))) #Store results

### Get results from standard frequentist approach
  lme4_res <- confint(lme4::lmer(formula = bmlmerform,data = mixdata))
  lme4_res[1:2,] <- lme4_res[1:2,]^2 
  if(rownames(lme4_res)[1]==".sig01"){lme4_res[1:2,] <- lme4_res[c(2,1),]}
  
  #eval[1,] <- eval[1,] + ((bmlmer_res$Quantiles[1,] < par_vec)&
  #                          (bmlmer_res$Quantiles[3,] > par_vec))
  
  
  lme4_res <- t(rbind(NA,lme4_res,2))
  eval <- rbind(eval,lme4_res)# store results
  
  #eval[2,] <- eval[2,] + ((lme4_res[,1] < par_vec)&
  #                          (lme4_res[,2] > par_vec))
  
  
 # print(bmlmer_res$Quantiles)
 # print(lme4_res )
  #print(c(i,nrep))
  if(i %in% seq(100,nrep,by = 100))
  {
    print('Bayes')
    print(rbind(
      colSums((eval[seq(1,dim(eval)[1],by = 4),-7]<Par_mat[1:i,])&
                (eval[seq(2,dim(eval)[1],by = 4),-7]>Par_mat[1:i,]))/i,
      apply(eval[seq(2,dim(eval)[1],by = 4),-7]-eval[seq(1,dim(eval)[1],by = 4),-7],2,mean)))
    
    print('LME4')
    print(rbind(
      colSums((eval[seq(3,dim(eval)[1],by = 4),-7]<Par_mat[1:i,])&
                (eval[seq(4,dim(eval)[1],by = 4),-7]>Par_mat[1:i,]))/i,
      apply(eval[seq(4,dim(eval)[1],by = 4),-7]-eval[seq(3,dim(eval)[1],by = 4),-7],2,mean)))
  }
}

plot(eval[seq(1,dim(eval)[1],by = 4),1],ylim = c(0,1),pch = '_')
points(eval[seq(2,dim(eval)[1],by = 4),1],ylim = c(0,1),pch = '_')
abline(h = delta,col = 'red')
plot(eval[seq(1,dim(eval)[1],by = 4),2],ylim = c(0,5),pch = '_')
points(eval[seq(2,dim(eval)[1],by = 4),2],ylim = c(0,1),pch = '_')
abline(h = sigma2,col = 'red')
plot(eval[seq(1,dim(eval)[1],by = 4),3],ylim = c(0,5),pch = '_')
points(eval[seq(2,dim(eval)[1],by = 4),3],ylim = c(0,1),pch = '_')
abline(h = sigma2u,col = 'red')

Par_mat <- matrix(par_vec,nrep,6,byrow=TRUE)

rbind(
colSums((eval[seq(1,dim(eval)[1],by = 4),-7]<Par_mat)&
          (eval[seq(2,dim(eval)[1],by = 4),-7]>Par_mat))/nrep,
apply(eval[seq(2,dim(eval)[1],by = 4),-7]-eval[seq(1,dim(eval)[1],by = 4),-7],2,mean))

rbind(
colSums((eval[seq(3,dim(eval)[1],by = 4),-7]<Par_mat)&
          (eval[seq(4,dim(eval)[1],by = 4),-7]>Par_mat))/nrep,
apply(eval[seq(4,dim(eval)[1],by = 4),-7]-eval[seq(3,dim(eval)[1],by = 4),-7],2,mean))



