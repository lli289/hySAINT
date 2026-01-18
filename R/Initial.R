#' Creating initial parents\cr
#'
#' This function creates the initial population for the genetic algorithm
#' by sampling from pre-screened main effects and interactions.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param EVAoutput The output from function \code{EVA}
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param r1 At most how many main effects do you want to include in your model?.
#'        For high-dimensional data, \code{r1} cannot be larger than the number of
#'        screened main effects.
#' @param r2 At most how many interaction effects do you want to include in your model?
#' @param numElite Number of elite parents. Default is 40.
#'
#' @return Initial parents. A numeric matrix with dimensions \code{numElite} by \code{r1+r2}.
#' @export
#'
#' @seealso \code{\\link{EVA}}
#'
#' @examples
#' set.seed(0)
#' X <- matrix(rnorm(100*10,1,0.1), 100, 10)
#' epl <- rnorm(100,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,3]+X[,1]*X[,2]+X[,1]*X[,3]+epl
#' EVAoutput <- EVA(X, y, r1 = 5, sigma = 0.01)
#' myParent <- Initial(X = X, y = y, EVAoutput, r1 = 5, r2 = 2)

Initial <- function(X, y, EVAoutput, heredity = "Strong", r1, r2, numElite = 40){

  p <- ncol(X)
  interaction.ind <- make_pairs(p)

  # Extract M_act and I_act from EVAoutput with proper sizing
  M_act <- EVAoutput$ranked.mainpool[1:min(r1*2, length(EVAoutput$ranked.mainpool))]

  # Extract interaction pairs from ranked.intermat
  if (nrow(EVAoutput$ranked.intermat) > 0) {
    I_act <- EVAoutput$ranked.intermat[1:min(r2*2, nrow(EVAoutput$ranked.intermat)), 1:2, drop = FALSE]
  } else {
    I_act <- matrix(ncol = 2, nrow = 0)
  }

  # Use init_population_matrix for efficient population initialization
  initial_parents <- init_population_matrix(
    M_act = M_act,
    I_act = I_act,
    interaction.ind = interaction.ind,
    p = p,
    N = numElite,
    r1 = c(1L, as.integer(r1)),
    r2 = c(0L, as.integer(r2)),
    heredity = heredity,
    seed = NULL
  )

  return(initial_parents)
}
