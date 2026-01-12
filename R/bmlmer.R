#' @title Bayesian Analyses of Simple and Balanced Linear Mixed-Effects Models. 
#' @description Bayesian analysis of simple linear mixed-effects model (LMM) applied 
#' to balanced data, via a closed form posterior distribution.
#' @details See Gangsei & Vinje ... for details.  
#' @usage bmlmer(formula, data, nu1,mu1,nu2,mu2,beta0,nu0,Upsilon0)
#' @param a bla bla
#' @return bla bla 
#' 
#' @note bla bla
#' 
#' @references Gangsei, L.E.& Vinje, H. 2025. A closed form solution for 
#' Bayesian analysis of a simple linear mixed model.
#' 
#'
#' @author Hilde Vinje\cr
#' hilde.vinje@@nmbu.no\cr
#' 
#' Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' 
#' @import lme4
#' @import mvtnorm
#' 
#' @examples
#' \dontrun{
#' 
#' ## Simulate data
#' set.seed(11)
#' nrep <- 10000
#' w <- 4
#' p <- 3
#' n <- 15
#' new_seeds <- sample(1:(10*nrep),nrep)
#' sigma2 <- 1.5
#' sigma2u <- 0.75
#' beta <- as.matrix(rnorm(p))
#' par_vec <- c(sigma2,sigma2u,as.vector(beta))
#' bmlmerform <- formula(y ~ X + (1|id))
#' ID <- matrix(rep(1:n,w),n,w,byrow=FALSE)
#' 
#' eval <- matrix(0,2,5, dimnames = 
#'           list(c('bmlmer','lme4'),
#'           c('sigma2','sigma2u','(intercept)','x1','x2')))
#' 
#' for(i in 1:nrep)
#' {
#' set.seed(new_seeds[i])
#' X <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1))
#' U <- matrix(rep(rnorm(n,0,sqrt(sigma2u)),w),n,w,byrow=FALSE)
#' E <- matrix(rnorm(n*w,0,sqrt(sigma2)),n,w)
#' Y <- matrix(X%*%beta,n,w,byrow=FALSE) + U + E
#' 
#' mixdata <- data.frame(y = as.vector(Y),
#' X = I(kronecker(as.matrix(rep(1,w)),X[,-1])),
#' id = as.vector(ID))
#' 
#' bmlmer_res <- bmlmer(formula = bmlmerform,data = mixdata,
#'                       nsim = 5000,nu1=0.01,nu2=0.01,beta0=NULL,nu3=0.01)
#'
#' lme4_res <- confint(lme4::lmer(formula = bmlmerform,data = mixdata))
#' 
#' eval[1,] <- eval[1,] + ((bmlmer_res$Quantiles[1,] < par_vec)&
#'                         (bmlmer_res$Quantiles[3,] > par_vec))
#'                         
#' eval[2,] <- eval[2,] + ((lme4_res[,1] < par_vec)&
#'                         (lme4_res[,2] > par_vec))
#' }
#' 
#' print(eval/nrep)
#' 
#' }
#'
#'
#'
#' @export
# Probability density distribution (pdf)
bmlmer <- function(mod_form,data,nu1=1,mu1=0.5,nu2=1,mu2=1,
                   beta0=NULL,nu3=1,Upsilon0=NULL,nsim = 10^4,
                   alpha = 0.05,simreturn = FALSE,
                   empirical_bayes = FALSE,nu_start = NULL,
                   nu_max = NULL,nu_min = NULL)
{
  # Split in random terms and fixed terms
  fixed_form <- reformulas::nobars(mod_form)
  random_terms <- reformulas::findbars(mod_form)
  
  grouping_factors <- sapply(random_terms, function(term) {
    as.character(term[[3]])})
  
  # Get all in matrix form
  id <- data[,grouping_factors]
  
  nn <- length(unique(id))
  ww <- dim(data)[1]/nn
  
  if(any(table(id)!=ww)){stop('The design has to be balanced')}
  
  # The response in one matrix
  YY <- matrix(data[order(id),as.character(fixed_form[[2]])],nn,ww,
               byrow = TRUE)
  
  if(any(is.na(YY))){stop('Missing values in response not allowed')}
  
  # Average per grouping factor
  yy_bar <- rowMeans(YY)
  
  # Get beta_hat (response variable estimates)
  fixed_mod <- lm(fixed_form,data = data)
  
  # The predictors in matrix form
  XX <- model.matrix(fixed_mod)
  XX <- XX[order(id)[seq(1,nn*ww,by = ww)],]
  
  Mn <- t(XX)%*%XX/nn
  
  if(is.null(beta0)){beta0 <- matrix(0,dim(XX)[2],1)}
  
  # Get different quadratic sums
  QQ1 <- sum((YY-matrix(yy_bar,nn,ww,byrow = FALSE))^2)
  
  beta_hat <- as.matrix(coef(fixed_mod))
  QQ2 <- ww*sum((yy_bar - XX%*%beta_hat)^2)
  
  # Use empirical Bayes to set nu-parameters 
  if(empirical_bayes == TRUE)
  {
    eb_res <- bm_evidence_max(mu1=mu1,mu2=mu2,beta0=beta0,
                              beta_hat=beta_hat,n=nn,w=ww,Q1=QQ1,Q2=QQ2,
                              Mn=Mn,nu_start = nu_start,nu_max = nu_max,
                              nu_min = nu_min)
    #nu1 <- eb_res$estimate[1]
    #nu2 <- eb_res$estimate[2]
    #nu3 <- eb_res$estimate[3]
    nu1 <- eb_res$par[1]
    nu2 <- eb_res$par[2]
    nu3 <- eb_res$par[3]
  }
  
  if(is.null(Upsilon0)){Upsilon0 <- nu3*Mn}
  QQ3 <- t(beta_hat-beta0)%*%(nn*Mn)%*%solve(nn*Mn + Upsilon0)%*%Upsilon0%*%(beta_hat-beta0)
  
  # Get posterior hyperparameters for extended beta distribution of 3th kind
  kappa1 <- as.numeric((QQ1 + QQ2 + 2*nu2/mu2 + ww*QQ3)/2)
  kappa2 <- as.numeric((QQ2 + ww*QQ3)/2)
  #print(kappa2/kappa1)
  phi1_post <- nn*ww/2 + nu2
  phi2_post <- nu1*mu1
  phi3_post <- nn/2 + nu1 
  
  #print(c(phi1_post,phi2_post,phi3_post,kappa1,kappa2))
  
  # Sample from the posterior for delta
  delta_post <- rbeta3ext(n = nsim,shape1 = phi1_post,shape2 = phi2_post,
                          shape3 = phi3_post,delta = kappa2/kappa1)
  
  # Hyperparameters for 1/sigma2
  shape1_post <- phi1_post
  rate1_post <- kappa1-delta_post*kappa2
  
  # Draw samples for 1/sigma2, calculate sigma2 and sigmau2
  sigma2inv_post <- rgamma(nsim,shape = shape1_post,rate = rate1_post)
  sigma2_post <- 1/sigma2inv_post 
  sigma2u_post <- delta_post*sigma2_post/((1-delta_post)*ww)
  
  
  # Draw samples for beta
  beta_weights <- sqrt(1/(ww*(1-delta_post)*sigma2inv_post))
  
  Error0_post <- (mvtnorm::rmvnorm(n = nsim, mean = rep(0,dim(XX)[2]),
                                  sigma = solve(nn*Mn + Upsilon0))*
                    matrix(beta_weights,nsim,dim(XX)[2],byrow = FALSE))
  
  beta_tilde <- solve(nn*Mn + Upsilon0)%*%((nn*Mn)%*%beta_hat+Upsilon0%*%beta0)
  
  Beta_post <- matrix(beta_tilde,nsim,dim(XX)[2],byrow = TRUE) + Error0_post
  
  # Send the result in an appropriate way
  ResMat <- cbind(delta_post,sigma2_post,sigma2u_post,Beta_post)
  colnames(ResMat) <- c('delta','sigma2','sigma2u',rownames(beta_hat))
  
  # Find and return approximate posteriors
  pdf_s <- lapply(as.data.frame(ResMat),density)
  
  posteriors <- lapply(pdf_s,function(x){approxfun(x$y,x$x, rule = 2)})
  
  res <- list(Mean = colMeans(ResMat),
              Quantiles = apply(ResMat,2,quantile,
                                         probs = c(alpha/2,0.5,1-alpha/2)),
              Covariance = cov(ResMat),
              Correlation = cor(ResMat))
  
  if(simreturn == TRUE){res$simulations <- data.frame(sigma2 = sigma2_post,
                                                sigma2u = sigma2u_post,
                                                Beta = I(Beta_post))}
  
  if(empirical_bayes == TRUE){res$eb_nu <- eb_res$par}
  
  return(res)
}