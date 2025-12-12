## Simulate data
library(beta3ext)
set.seed(15)
nrep <- 1000
w <- 4
p <- 3
n <- 15
new_seeds <- sample(1:(10*nrep),nrep)
sigma2 <- 4
sigma2u <- 0.25
delta <- w*sigma2u/(sigma2 + w*sigma2u)
beta <- as.matrix(rnorm(p))
par_vec <- c(delta,sigma2,sigma2u,as.vector(beta))
Par_mat <- matrix(par_vec,nrep,length(par_vec),byrow = TRUE)

bmlmerform <- formula(y ~ X + (1|id))
ID <- matrix(rep(1:n,w),n,w,byrow=FALSE)

eval <- matrix(0,0,7, dimnames = 
                 list(NULL,
                      c('delta','sigma2','sigma2u','(intercept)','x1','x2','method')))

for(i in 1:nrep)
{
  set.seed(new_seeds[i])
  X <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1))
  U <- matrix(rep(rnorm(n,0,sqrt(sigma2u)),w),n,w,byrow=FALSE)
  E <- matrix(rnorm(n*w,0,sqrt(sigma2)),n,w)
  Y <- matrix(X%*%beta,n,w,byrow=FALSE) + U + E
  
  mixdata <- data.frame(y = as.vector(Y),
                        X = I(kronecker(as.matrix(rep(1,w)),X[,-1])),
                        id = as.vector(ID))
  
  bmlmer_res <- bmlmer(formula = bmlmerform,data = mixdata,
                       nsim = 5000,beta0=NULL,nu1 = 2,nu2 = 0.01,nu3 = 0.01,
                       simreturn = FALSE,empirical_bayes = TRUE,nu_max = c(2.0001,10,10))
  
  print(bmlmer_res$eb_nu)
  
  lme4_res <- confint(lme4::lmer(formula = bmlmerform,data = mixdata))
  lme4_res[1:2,] <- lme4_res[1:2,]^2 
  if(rownames(lme4_res)[1]==".sig01"){lme4_res[1:2,] <- lme4_res[c(2,1),]}
  
  #eval[1,] <- eval[1,] + ((bmlmer_res$Quantiles[1,] < par_vec)&
  #                          (bmlmer_res$Quantiles[3,] > par_vec))
  eval <- rbind(eval,cbind(cbind(bmlmer_res$Quantiles[c(1,3),],1)))
  
  lme4_res <- t(rbind(NA,lme4_res,2))
  eval <- rbind(eval,lme4_res)
  
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



