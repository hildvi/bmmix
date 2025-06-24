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
#' +47 959 01 450 
#' 
#' Lars Erik Gangsei\cr
#' lars.erik.gangsei@@vetinst.no\cr
#' +47 950 61 231
#' 
#' @import lme4
#' 
#' @examples
#' \dontrun{
#' 
#' ## Simulate data
#' set.seed(11)
#' sigma2 <- 1.5
#' sigma2u <- 0.75
#' w <- 4
#' p <- 3
#' n <- 15
#' X <- cbind(1,matrix(rnorm(n*(p-1)),n,p-1))
#' beta <- as.matrix(rnorm(p))
#' U <- matrix(rep(rnorm(n,0,sqrt(sigma2u)),w),n,w,byrow=FALSE)
#' E <- matrix(rnorm(n*w),n,w)
#' Y <- matrix(X%*%beta,n,w,byrow=FALSE) + U + E
#' ID <- matrix(rep(1:n,w),n,w,byrow=FALSE)
#' 
#' mixdata <- data.frame(y = as.vector(Y),
#' X = I(kronecker(as.matrix(rep(1,w)),X[,-1])),
#' id = as.vector(ID))
#' 
#' bmlmerform <- formula(y ~ X + (1|id))
#' 
#' }
#'
#'
#'
#' @export
# Probability density distribution (pdf)
bmlmer <- function(formula, data,nu1=1,mu1=1,nu2=1,mu2=1,
                   beta0=NULL,nu0=1,Upsilon0=NULL)
{
  # Split in random terms and fixed terms
  fixed_form <- lme4::nobars(formula)
  random_terms <- lme4::findbars(formula)
  
  grouping_factors <- sapply(random_terms, function(term) {
    as.character(term[[3]])
  })
  
  # Get all in matrix form
  id <- data[,grouping_factors]
  nn <- length(unique(id))
  ww <- dim(data)[1]/nn
  YY <- matrix(data[order(id),as.character(fixed_form[[2]])],nn,ww,
               byrow = TRUE)
  
  fixed_form
  
  
  res <- 1
}