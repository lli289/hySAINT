#' ABC Evaluation\cr
#'
#' Gives ABC score for each fitted model. For a model I, the ABC is defined as
#' \deqn{ABC(I)=\sum\limits_{i=1}^n\bigg(Y_i-\hat{Y}_i^{I}\bigg)^2+2r_I\sigma^2+\lambda\sigma^2C_I.}
#' When comparing ABC of fitted models to the same dataset, the smaller the ABC, the better fit.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects.
#' @param y Response variable. A \code{n}-dimensional vector.
#' @param heredity Whether to enforce Strong, Weak, or No heredity. Default is "Strong".
#' @param sigma The standard deviation of the noise term. In practice, sigma is usually
#'        unknown. Users can estimate sigma using various methods such as the MAD estimator:
#'        \code{sigma = mad(residuals(lm(y ~ X)))}, or from external packages.
#'        See examples for details.
#' @param varind A numeric vector that specifies the indices of variables to be extracted from \code{X}.
#'        Default is "No".
#' @param lambda A numeric value defined by users. The number needs to satisfy the condition:
#'        \eqn{\lambda\geq 5.1/log(2)}. Default is 10.
#'
#' @return A numeric value is returned. It represents the ABC score of the fitted model.
#' @export
#'
#' @references
#' Ye, C. and Yang, Y., 2019. \emph{High-dimensional adaptive minimax sparse estimation with interactions.}
#'
#' @examples # When sigma is known
#' set.seed(0)
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y <- 1+X[,1]+X[,2]+X[,1]*X[,2] + epl
#'  ABC(X, y, sigma = 0.01, varind = c(1,2,5))
#'
#' @examples # When sigma is not known
#' full <- Extract(X, varind = c(1:(dim(X)[2]+choose(dim(X)[2],2))))
#' sigma_est <- mad(residuals(lm(y ~ full)))  # Estimate sigma using MAD

ABC <- function(X, y, heredity = "Strong", sigma, varind = NULL, lambda = 10){
  # Filter out zeros and NAs before checking duplicates
  valid_varind <- varind[!is.na(varind) & varind != 0]
  if (length(valid_varind) > 0 && any(duplicated(valid_varind))){
    # Remove duplicates instead of stopping
    varind <- varind[!duplicated(varind) | varind == 0]
  }

  p <- dim(X)[2]
  interaction.ind <- make_pairs(p)
  pi0 <- 1e-10
  pi1 <- pi2 <- pi3 <- (1-pi0)/3
  n <- dim(X)[1]
  sigma <- sigma
  data_extract <- Extract(X, varind, interaction.ind)
  data <- data_extract
  
  # Check if data has valid dimensions
  if (ncol(data) == 0 || nrow(data) == 0) {
    # Return a very large penalty if no valid data
    return(1e10)
  }
  
  r.I <- Matrix::rankMatrix(data)[1]

  # Extract k1 (main effects) and k2 (interactions)
  k2 <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(data))) > p)
  k1 <- sum(as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(data))) <= p)
  
  # Use lm.fit for better numerical stability and efficiency
  if (r.I < ncol(data)){
    # Handle rank-deficient case with orthogonalization
    data_orth <- pracma::orth(data)
    fit <- lm.fit(data_orth, y)
    yhat <- fitted(fit)
  } else {
    # Full rank case - direct lm.fit
    fit <- lm.fit(data, y)
    yhat <- fitted(fit)
  }
  SSE <- sum((y - yhat)^2)
  # Calculate model complexity C_I based on heredity condition
  # Formula: ABC(I) = SSE + 2*r_I*sigma^2 + lambda*sigma^2*C_I
  if (heredity == "Strong"){
    # Strong heredity: interaction requires both parents
    C.I.strong <- -log(pi1) + log(min(p, n)) + log(min(mychoose(k1), n)) + 
                  log(choose(p, k1)) + log(choose(mychoose(k1), k2))
    ABC <- SSE + 2*r.I*sigma^2 + lambda*sigma^2*C.I.strong
  }
  else if (heredity == "Weak"){
    # Weak heredity: interaction requires at least one parent
    K <- k1*p - choose(k1, 2) - k1
    C.I.weak <- -log(pi2) + log(min(p, n)) + log(min(K, n)) + 
                log(choose(p, k1)) + log(choose(K, k2))
    ABC <- SSE + 2*r.I*sigma^2 + lambda*sigma^2*C.I.weak
  }
  else if (heredity == "No"){
    # No heredity: interaction independent of parents
    C.I.no <- -log(pi3) + log(min(p, n)) + log(min(choose(p, 2), n)) + 
              log(choose(p, k1)) + log(choose(choose(p, 2), k2))
    ABC <- SSE + 2*r.I*sigma^2 + lambda*sigma^2*C.I.no
  }
  return(ABC)
}
