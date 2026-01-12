################################################################################
## Script for making table 1 in Gangsei & Vinje (2026)
## Authors: Hilde Vinje and Lars Erik Gangsei

## Clean workspace
rm(list = ls())  

## Load the beta3ext and bmmix packages
#devtools::install_git(url = 'https://github.com/hildvi/bmmix')
#devtools::install_git(url = 'https://github.com/hildvi/beta3ext')


require(beta3ext)
require(bmmix)

## Set dimensions and parameters in model

nrep <- 1000 # Number of simulations

### Dimentions
w <- 4
p <- 3
n <- 20

### Set parameters
sigma2 <- 4
sigma2u <- 0.5

delta <- w*sigma2u/(sigma2 + w*sigma2u)
beta <- c(1,2,-0.5)#as.matrix(rnorm(p))

### Put parameters in one vector
par_vec <- c(delta,sigma2,sigma2u,as.vector(beta))
Par_mat <- matrix(par_vec,nrep,length(par_vec),byrow = TRUE)

## Set up model and environment for evaluation
bmlmerform <- formula(y ~ X + (1|id))
ID <- matrix(rep(1:n,w),n,w,byrow=FALSE)

bayes_eval <- freq_eval <- matrix(NA,0,6, dimnames = 
                 list(NULL,
              c('delta','sigma2','sigma2u','(intercept)','x1','x2')))



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
                       simreturn = FALSE,empirical_bayes = TRUE,
                       nu_min = c(2,0,0),nu_max = c(2.0001,10,10))
  
  #print(bmlmer_res$eb_nu)
  bayes_eval <- rbind(bayes_eval,
                      rbind(bmlmer_res$Quantiles[c(1,3),],
                            bmlmer_res$Mean)) #Store results

### Get results from standard frequentist approach
  freq_mod <- lme4::lmer(formula = bmlmerform,data = mixdata)
  lme4_res <- confint(freq_mod)
  lme4_res[1:2,] <- lme4_res[1:2,]^2 
  
  if(rownames(lme4_res)[1]==".sig01"){lme4_res[1:2,] <- lme4_res[c(2,1),]}
  
  lme4_res <- t(rbind(NA,cbind(lme4_res,
                    c(as.data.frame(lme4::VarCorr(freq_mod))$vcov[c(2,1)],
                      summary(freq_mod)$coef[,1]))))
  
  #eval[1,] <- eval[1,] + ((bmlmer_res$Quantiles[1,] < par_vec)&
  #                          (bmlmer_res$Quantiles[3,] > par_vec))
  
  
  #lme4_res <- t(rbind(NA,lme4_res,2))
  freq_eval <- rbind(freq_eval,lme4_res)# store results
}

## Table for manuscript
Tab1 <- data.frame(
      Real_value = par_vec,
      Bayes_OK_per = colSums((bayes_eval[seq(1,3*nrep,by = 3),-7]<Par_mat)&
          (bayes_eval[seq(2,3*nrep,by = 3),-7]>Par_mat))/nrep,
      Bayes_CI_widh = apply(bayes_eval[seq(2,3*nrep,by = 3),-7]-
              bayes_eval[seq(1,3*nrep,by = 3),-7],2,mean),
      Bayes_MSE = apply(bayes_eval[seq(3,3*nrep,by = 3),-7]-
                          Par_mat,2,var),
      Bayes_bias = apply(bayes_eval[seq(3,3*nrep,by = 3),-7]-
                           Par_mat,2,mean),
      Freq_OK_per =colSums((freq_eval[seq(1,3*nrep,by = 3),-7]<Par_mat)&
                (freq_eval[seq(2,3*nrep,by = 3),-7]>Par_mat))/nrep,
      Bayes_CI_widh = apply(freq_eval[seq(2,3*nrep,by = 3),-7]-
              freq_eval[seq(1,3*nrep,by = 3),-7],2,mean),
      Freq_MSE = apply(freq_eval[seq(3,3*nrep,by = 3),-7]-
                         Par_mat,2,var),
      Freq_bias = apply(freq_eval[seq(3,3*nrep,by = 3),-7]-
                           Par_mat,2,mean))


### Write table in Latex format
xtable::xtable(Tab1)



