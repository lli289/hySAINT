#' Performing mutation\cr
#'
#' This function gives mutant from parents.
#'
#' @param myParent A numeric matrix with dimension \code{numElite} by \code{r1 + r2}.
#' @param EVAoutput The output from function \code{EVA}.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of
#'        screened main effects.
#' @param r2 At most how many interaction effects do you want to include in your model?
#' @param initial.temp Initial temperature. Default is 1000.
#' @param cooling.rate A numeric value represents the speed at which the
#'        temperature decreases. Default is 0.95.
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#'        unknown. Users can estimate sigma using various methods such as the MAD estimator:
#'        \code{sigma = mad(residuals(lm(y ~ X)))}.
#' @param varind A numeric vector that specifies the indices of variables to be extracted from \code{X}.
#' @param lambda A numeric value defined by users. The number needs to satisfy the condition:
#'        \eqn{\lambda\geq 5.1/log(2)}. Default is 10.
#'
#' @return Mutant. A numeric matrix with dimensions \code{numElite} by \code{r1+r2}.
#' @export
#'
#' @seealso \code{\link{EVA}}, \code{\link{Initial}}.
#'
#' @examples
#' set.seed(0)
#' X <- matrix(rnorm(100*10,1,0.1), 100, 10)
#' epl <- rnorm(100,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
#' EVAoutput <- EVA(X, y, r1 = 5, sigma = 0.01)
#' myParent <- Initial(X = X, y = y, EVAoutput, r1 = 5, r2 = 2)
#' Mutation(myParent, EVAoutput, r1 = 5, r2 = 2, X = X, y = y, sigma = 0.1)


Mutation <- function(myParent, EVAoutput, r1, r2, initial.temp = 1000, cooling.rate = 0.95, X, y, heredity = "Strong", sigma, varind = NULL, lambda = 10){

  interaction.ind <- make_pairs(dim(X)[2])

  mutate.or.not <- stats::rbinom(1,1,prob = 0.7)
  if (mutate.or.not == 1){
    mutant <- Saint(myParent, EVAoutput, r1, r2, initial.temp, cooling.rate, X, y,heredity, sigma, varind, interaction.ind, lambda)
  }else{
    mutant <- myParent
  }
  return(mutant)
}
