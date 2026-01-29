#' @title Bayesian Analyses of Simple and Balanced Linear Mixed-Effects Models. 
#' @description Function to set prior sample sizes based on maxiizing model evidence.
#' @details See Gangsei & Vinje (2026) for details.  
#' @usage bm_evidence_max(mu1,mu2,beta0,n,w,Q1,Q2,Mn)
#' @param mu1,mu2,beta0. Prior hyper parameters for expected values.
#'        mu1 (0<mu1<1) and mu2 are prior expectations
#'        for the beta and gamma distributions respectively. beta0 (length p) is prior
#'        expected value for fixed regression parameters. 
#' @param n,w integers, number of independent samples and number of replicates 
#'        per sample respectively.
#' @param Q1,Q2 positive scalars. Quadratic sms from data.
#' @param Mn positive definite matrix (p x p). Mn = X'X/n
#' @param nu_start,nu_max,nu_min optional. All vectors of length 3 setting start, 
#'        minimum and maximum numbers for the three prior sizes (nu1, nu2, and nu3)
#'        if empirical Bayes is applied to set hyperparameters. 
#'        
#' @return The function returns a list of the same class as returned by function 
#'        optim(). 
#' 
#' @note Uses function optim
#' 
#' @references Gangsei, L.E.& Vinje, H. 2026. A closed form solution for 
#' Bayesian analysis of a simple linear mixed model.
#' 
#'
#' @author Hilde Vinje\cr
#' hilde.vinje@@nmbu.no\cr
#' 
#' Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' 
#' @import gsl
#' 
#' @examples
#' \dontrun{
#' 
#' }
#'
#'
#'
#' @export
# Probability density distribution (pdf)
bm_evidence_max <- function(mu1=0.5,mu2=1,beta0=NULL,beta_hat,n,w,Q1,Q2,Mn,
                            nu_start = NULL,nu_max = NULL,nu_min = NULL)
{
  # Set beta0 till 0 if not defined before
  if(is.null(beta0)){beta0 <- as.matrix(rep(0,dim(Mn)[1]))}
  if(is.null(nu_start)){nu_start <- c(2,1,1)}
  if(is.null(nu_max)){nu_max <- rep(n/2,3)}
  if(is.null(nu_min)){nu_min <- c(2,0,0)}
  
  # Function to minimze the logarithm of log likelihood.
  ev_min <- function(nu,mu1,mu2,beta0,beta_hat,n,w,Q1,Q2,Mn,
                     nu_min = rep(0,3),nu_max=rep(Inf,3))
  {
    if(any(nu<nu_min)&&any(nu>nu_max)){res <- -10^50
    }else{
     beta_hat <- as.matrix(beta_hat)
     beta0 <- as.matrix(beta0)
     Mn <- as.matrix(Mn)
     
     Q3 <- (nu[3]/(n+nu[3]))*t(beta_hat-beta0)%*%(n*Mn)%*%(beta_hat-beta0)
     kappa1 <- (Q1+Q2+2*nu[2]/mu2+w*Q3)/2
     kappa2 <- (Q2+w*Q3)/2
     phi1_post <- n*w/2 + nu[2]
     phi2_post <- nu[1]*mu1
     phi3_post <- n/2 + nu[1]*(1-mu1)
     p <- dim(Mn)[1]
     
     res <- (-(n*w/2)*log(2*pi)
              +log(gsl::hyperg_2F1(phi2_post,phi1_post,phi3_post+phi2_post,kappa2/kappa1))
              -(phi1_post)*log(kappa1)
              +(p/2)*(log(nu[3])-log(n+nu[3]))
              +nu[2]*(log(nu[2])-log(mu2))
              +lgamma(phi1_post)
              -lgamma(nu[2])
              -lgamma(phi3_post+phi2_post)
              +lgamma(phi3_post)
              +lgamma(nu[1])
              -lgamma(nu[1]*(1-mu1)))
     if(is.na(res)||abs(res)==Inf){res <- -10^50}
    }
     return(-res)
  }
  
  # Use the non-linar optimizer nlm() to solve for nu
  #res <- nlm(f = ev_min,p = nu_start,mu1=mu1,mu2=mu2,beta0=beta0,
  #           beta_hat = beta_hat,n=n,w=w,Q1=Q1,Q2=Q2,Mn = Mn,nu_min = nu_min,
  #           nu_max = nu_max)
  res <- optim(par = nu_start,fn = ev_min,gr = NULL,
               mu1=mu1,mu2=mu2,beta0=beta0,beta_hat = beta_hat,n=n,w=w,Q1=Q1,
               Q2=Q2,Mn = Mn,method = "L-BFGS-B",nu_min = nu_min,nu_max = nu_max,
               lower = nu_min,upper = nu_max)
  #print(res)
  return(res)
}