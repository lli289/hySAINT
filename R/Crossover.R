#' Performing crossover\cr
#'
#' This function gives offspring from parents. It performs crossover at a fixed
#' probability of 0.6.
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param myParent A numeric matrix with dimension \code{numElite} by \code{r1 + r2}.
#' @param EVAoutput The output from function \code{EVA}.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of screened main effects.
#' @param r2 At most how many interaction effects do you want to include in your model?
#' @param numElite Number of elite parents. Default is 40.
#'
#' @return Offspring. If crossover occurred, it returns a numeric matrix with dimensions
#' \code{choose(numElite,2)} by \code{r1+r2}. Otherwise, \code{numElite} by \code{r1 + r2}.
#'
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
#' Offsprings <- Crossover(X, myParent, EVAoutput, r1 = 5, r2 = 2)


Crossover <- function(X, myParent, EVAoutput, heredity = "Strong", r1, r2, numElite = 40){

  # DETERMINE WHETHER TO PERFORM CROSSOVER
  cross.or.not <- stats::rbinom(1,1,prob = 0.6)
  if (cross.or.not == 1){
    offsprings <- sample_prob(X, myParent, EVAoutput, heredity, r1, r2)
  }else{
    offsprings <- myParent
  }
  return(offsprings)
}
