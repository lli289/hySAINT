#' Extracting columns and generating required interaction effects from data\cr
#'
#' This function simplifies the data preparation process by enabling users to
#' extract specific columns from their dataset \code{X}, and automatically
#' generating any necessary interaction effects based on \code{varind}.
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#'        \code{n} observations by \code{p} main effects. Note that the interaction
#'        effects should not be included in \code{X} because this function
#'        automatically generates the corresponding interaction effects if needed.
#' @param varind A numeric vector that specifies the indices of variables to be
#'        extracted from \code{X}. Duplicated values are not allowed. See Example
#'        for details.
#' @param interaction.ind A two-column numeric matrix. Each row represents a unique
#'        interaction pair, with the columns indicating the index numbers of the variables
#'        involved in each interaction. If NULL (default), will be generated automatically.
#'
#' @return A numeric matrix is returned.
#' @export
#'
#' @examples # Generate data
#' set.seed(0)
#' X <- matrix(rnorm(20), ncol = 4)
#' y <- X[, 2] + rnorm(5)
#'
#' @examples # Extract X1 and X1X2 from X1, ..., X4
#' Extract(X, varind = c(1,5))
#'
#' @examples # Extract X5 from X1, ..., X4
#' Extract(X, varind = 5)
#'
#' @examples # Extract using duplicated values
#' try(Extract(X, varind = c(1,1))) # this will not run
#'
Extract <- function(X, varind, interaction.ind = NULL){
  if (is.null(interaction.ind)){
    interaction.ind <- make_pairs(ncol(X))
  }
  # Filter out zeros and NAs before checking duplicates
  valid_varind <- varind[!is.na(varind) & varind != 0]
  if (length(valid_varind) > 0 && any(duplicated(valid_varind))){
    # Remove duplicates instead of stopping
    varind <- varind[!duplicated(varind) | varind == 0]
  }
  ncoln <- length(varind)
  nrown <- nrow(X)
  p <- ncol(X)
  mainind <- varind[which(varind%in%1:p)]
  
  # Handle main effects
  if (length(mainind) > 0) {
    if (length(mainind) == 1) {
      mainvars.matrix <- matrix(X[, mainind], ncol=1)
    } else {
      mainvars.matrix <- X[, mainind, drop=FALSE]
    }
  } else {
    mainvars.matrix <- NULL
  }
  
  # Handle interaction effects
  interind <- varind[which(varind>p)]
  
  if (length(interind) > 0) {
    a <- interaction.ind[interind-p, , drop=FALSE]
    intervars.matrix <- matrix(0, nrow = nrown, ncol = length(interind))
    
    for (i in 1:length(interind)) {
      if (nrow(a) == 1) {
        intervars.matrix[,i] <- X[,a[1,1]]*X[,a[1,2]]
      } else {
        intervars.matrix[,i] <- X[, a[i,1]]*X[,a[i,2]]
      }
    }
  } else {
    intervars.matrix <- NULL
  }
  
  # Combine matrices
  if (!is.null(mainvars.matrix) && !is.null(intervars.matrix)) {
    data_extract <- cbind(mainvars.matrix, intervars.matrix)
  } else if (!is.null(mainvars.matrix)) {
    data_extract <- mainvars.matrix
  } else if (!is.null(intervars.matrix)) {
    data_extract <- intervars.matrix
  } else {
    data_extract <- matrix(0, nrow = nrown, ncol = 0)
  }
  
  # Ensure data_extract is a matrix
  data_extract <- as.matrix(data_extract)
  
  # Set column names safely
  expected_cols <- length(mainind) + length(interind)
  if (ncol(data_extract) == expected_cols && expected_cols > 0) {
    col_names <- character(expected_cols)
    if (length(mainind) > 0) {
      col_names[1:length(mainind)] <- paste0("X", mainind)
    }
    if (length(interind) > 0) {
      col_names[(length(mainind)+1):expected_cols] <- paste0("X", interind)
    }
    colnames(data_extract) <- col_names
  }
  return(data_extract)
}

# BELOW ARE THE HELPER FUNCTIONS THAT DO NOT NEED TO BE DOCUMENTED

# OPTIMIZED: Vectorized interaction pair generation - much faster than loops
make_pairs <- function(p){
  if (p < 2) return(matrix(ncol = 2, nrow = 0))
  # Vectorized approach - creates all pairs efficiently
  i <- rep.int(seq_len(p - 1L), times = rev(seq_len(p - 1L)))
  j <- sequence(rev(seq_len(p - 1L))) + rep.int(seq_len(p - 1L), times = rev(seq_len(p - 1L)))
  cbind(i, j)
}

# MODIFIED CHOOSE FUNCTION
mychoose <- function(k1){
  if (k1 == 1){
    return(1)
  }else{
    return(choose(k1,2))
  }
}

# MODIFIED DIM FUNCTION
mydim <- function(x){
  if (is.matrix(x)) dim(x)
  else return(1)
}

# MODIFIED SORT FUNCTION
sort_zeros <- function(vec){
  non_zeros <- vec[vec != 0]
  zeros <- vec[vec == 0]
  if (length(zeros) > 0){
    return(c(sort(non_zeros), zeros))
  }else{
    return(c(sort(non_zeros)))
  }
}

# MODIFIED SAMPLE FUNCTION
mysample <- function(x, size, replace = FALSE, prob = NULL){
  if (length(x) == 1) return(x)
  if (length(x) > 1) return(sample(x, size, replace, prob))
}

# Helper function for partial correlation (used in pre_screen)
partial_cor <- function(W, y, X1, X2) {
  # Compute partial correlation of W and y given X1, X2
  # This is cor(resid(W ~ X1 + X2), resid(y ~ X1 + X2))
  fit_W <- lm.fit(cbind(1, X1, X2), W)
  fit_y <- lm.fit(cbind(1, X1, X2), y)
  cor(fit_W$residuals, fit_y$residuals)
}

################################################################################
##  HELPER: PRE-SCREENING  (Step 0) -------------------------------------------
################################################################################
pre_screen <- function(X, y, m = 500, L = 1000,
                       pairs = make_pairs(ncol(X)),
                       verbose = TRUE) {

  p   <- ncol(X);           stopifnot(m <= p)
  cor_main <- abs(cor(X, y))
  M0  <- order(cor_main, decreasing = TRUE)[seq_len(m)]

  sel_pairs <- pairs[pairs[,1] %in% M0 & pairs[,2] %in% M0, , drop = FALSE]

  # Vectorize partial correlation computation where possible
  pcor <- numeric(nrow(sel_pairs))
  for (idx in seq_len(nrow(sel_pairs))) {
    j <- sel_pairs[idx, 1];  k <- sel_pairs[idx, 2]
    W <- X[, j] * X[, k]
    pcor[idx] <- abs(partial_cor(W, y, X[, j], X[, k]))
  }
  I0_idx <- order(pcor, decreasing = TRUE)[seq_len(min(L, length(pcor)))]
  I0     <- sel_pairs[I0_idx, , drop = FALSE]

  # Vectorized residual calculation
  Z <- matrix(NA_real_, nrow = nrow(X), ncol = nrow(I0))
  for (l in seq_len(nrow(I0))) {
    j <- I0[l, 1];  k <- I0[l, 2]
    W <- X[, j] * X[, k]
    Z[, l] <- lm.fit(cbind(1, X[, j], X[, k]), W)$resid
  }
  colnames(Z) <- paste0("Z_", I0[,1], "_", I0[,2])

  S <- cbind(X[, M0, drop = FALSE], Z)
  if (verbose)
    message("Pre-screened ", ncol(S), " variables (", length(M0), " mains, ",
            ncol(Z), " interactions)")

  list(S      = S,
       M0     = M0,
       I0_mat = I0,
       pairs  = pairs)
}

################################################################################
##  HELPER: ELASTIC-NET  (Step 1) ---------------------------------------------
################################################################################
enet_step <- function(S_list, y, heredity = c("Strong", "Weak", "No"),
                      alpha = 0.5, lambda = NULL, seed = 1L) {

  heredity <- match.arg(heredity)
  S   <- S_list$S
  M0  <- S_list$M0
  I0  <- S_list$I0_mat          # L × 2
  pairs<- S_list$pairs

  set.seed(seed)
  fit <- glmnet::cv.glmnet(S, y, family = "gaussian",
                           alpha  = alpha,
                           lambda = lambda,
                           standardize = TRUE)
  coefs <- coef(fit, s = "lambda.min")[, 1]
  name_vec <- names(coefs)
  active <- which(abs(coefs) > 0 & name_vec != "(Intercept)")

  ## Map back to mains vs. interactions
  p_main <- length(M0)
  idx0 <- active - 1L

  act_main_idx <- idx0[idx0 <= p_main]
  act_int_idx  <- idx0[idx0 >  p_main] - p_main

  M_act <- M0[act_main_idx]
  I_act <- I0[act_int_idx, , drop = FALSE]

  if (heredity != "No" && nrow(I_act)) {
    parents <- unique(as.vector(I_act))
    if (heredity == "Strong") {
      M_act <- union(M_act, parents)
    } else {
      missing <- setdiff(parents, M_act)
      if (length(missing)) {
        for (row in seq_len(nrow(I_act))) {
          for (j in I_act[row, ]) {
            if (!(j %in% M_act)) { M_act <- c(M_act, j); break }
          }
        }
      }
    }
  }

  list(M_act = sort(unique(M_act)),
       I_act = I_act,
       coef   = coefs)
}

################################################################################
##  HELPER: INITIAL GA POPULATION  (Step 2) -----------------------------------
################################################################################
init_population_matrix <- function(M_act,
                                   I_act,
                                   interaction.ind,
                                   p,
                                   N        = 100,
                                   r1       = c(2L, 10L),
                                   r2       = c(0L, 10L),
                                   heredity = c("Strong", "Weak", "No"),
                                   seed     = NULL) {
  heredity <- match.arg(heredity)
  if (!is.null(seed)) set.seed(seed)

  full_lab <- apply(interaction.ind, 1, paste, collapse = "-")
  act_lab  <- apply(I_act,         1, paste, collapse = "-")
  full_row_idx <- match(act_lab, full_lab)
  code_act <- p + full_row_idx

  max_m    <- r1[2]
  max_i    <- r2[2]
  max_vars <- max_m + max_i
  pop_mat  <- matrix(0L, nrow = N, ncol = max_vars)

  for (i in seq_len(N)) {
    size_m <- sample(seq(r1[1], min(max_m, length(M_act))), 1)
    size_i <- sample(seq(r2[1], min(max_i, nrow(I_act))), 1)

    M_sel <- sort(sample(M_act, size_m))
    I_idx <- if (size_i > 0) sort(sample(seq_len(nrow(I_act)), size_i)) else integer(0)

    if (length(I_idx) > 0 && heredity == "Strong") {
      pairs <- I_act[I_idx, , drop = FALSE]
      ok    <- apply(pairs, 1, function(par) all(par %in% M_sel))
      I_idx <- I_idx[ok]
    } else if (length(I_idx) > 0 && heredity == "Weak") {
      pairs <- I_act[I_idx, , drop = FALSE]
      ok    <- apply(pairs, 1, function(par) any(par %in% M_sel))
      I_idx <- I_idx[ok]
    }

    inter_codes <- code_act[I_idx]

    var_ind <- sort(c(M_sel, inter_codes))
    pop_mat[i, seq_along(var_ind)] <- var_ind
  }

  pop_mat
}

# FUNCTION USED TO REMOVE DUPLICATED ROWS
remove.dup.rows <- function(Matrix, r1, r2){
  duplicates <- duplicated(Matrix) | duplicated(Matrix, fromLast = TRUE)
  any_duplicates <- any(duplicates)
  if (any_duplicates == TRUE){
    row_dup <- as.numeric(which(duplicates))
    rm_ind <- row_dup[seq(1, length(row_dup), by = 2)]
    no_dup_new <- Matrix[-rm_ind,]
    Matrix <- no_dup_new[order(unlist(no_dup_new[,((r1+r2)+1)]), na.last = TRUE),]
  }else{
    Matrix <- Matrix
  }
  return(Matrix)
}


# FUNCTION USED TO MATCH IDX TO VARIABLE NAMES
predictor_match <- function(candidate.model, p, r1,r2,interaction.ind){
  model.match <- list()
  if (mydim(candidate.model)[1]==1){
    ee <- as.numeric(candidate.model[1:(r1+r2)])
    ee1 <- as.numeric(ee[which(ee%in% 1:p)])
    ee2 <- as.numeric(ee[which(ee > p)])
    if (length(ee2) == 0){
      model.match[[1]] <-c(paste0("X", ee1))
    }
    if (length(ee1) == 0){
      model.match[[1]] <- c( paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
    }
    if (!length(ee2) ==0 & !length(ee1) ==0){
      model.match[[1]] <-c(paste0("X", ee1), paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
    }
  }else{
    for (i in 1:nrow(candidate.model)) {
      ee <- candidate.model[i, 1:(r1+r2)]
      ee1 <- as.numeric(ee[which(ee%in% 1:p)])
      ee2 <- as.numeric(ee[which(ee > p)])
      if (length(ee2) == 0){
        model.match[[i]] <-c(paste0("X", ee1))
      }
      if (length(ee1) == 0){
        model.match[[i]] <- c( paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
      }
      if (!length(ee2) ==0 & !length(ee1) ==0){
        model.match[[i]] <-c(paste0("X", ee1), paste0("X", interaction.ind[ee2-p,1], "X", interaction.ind[ee2-p,2]))
      }
    }
  }
  return(model.match = model.match)
}

# FUNCTION USED TO MATCH MAIN EFFECT
mainRank_match <- function(a,b){

  ranks <- numeric(length(a))
  for (i in 1:length(a)) {
    ranks[i] <- which(b == a[i])
  }
  return(list(mainVec = b,
              mainSIS = a,
              ranks = ranks))
}

# OPTIMIZED: FUNCTION USED DURING CROSSOVER - Vectorized and more efficient
sample_prob <- function(X, myParent, EVAoutput, heredity = "Strong", r1, r2){
  p <- dim(X)[2]
  InterRank.ind.mat <- EVAoutput$ranked.intermat
  max_model_size <- dim(myParent)[2]
  n_rows <- nrow(myParent)
  
  # Compute row probabilities
  top_20_percent <- ceiling(0.2*n_rows)
  row_probabilities <- rep(0.1, n_rows)
  row_probabilities[1:top_20_percent] <- 0.9
  
  # Sample even number of rows
  if (n_rows == 3 | n_rows == 2){
    num_rows_to_sample <- 2
  }else{
    even_numbers <- seq(2, (n_rows-2), by = 2)
    num_rows_to_sample <- mysample(even_numbers, 1)
  }
  
  sampled_rows <- myParent[sample(n_rows, size = num_rows_to_sample, 
                                   replace = TRUE, prob = row_probabilities), , drop = FALSE]
  CPMatrix <- sampled_rows
  n_pairs <- choose(dim(CPMatrix)[1], 2)
  
  # Always create a proper matrix
  if (n_pairs == 0) n_pairs <- 1
  CMatrix <- matrix(0, nrow = n_pairs, ncol = max_model_size)
  
  # OPTIMIZED: Vectorize crossover operations where possible
  pair_count <- 0
  for (i in 1:(dim(CPMatrix)[1]-1)) {
    for (j in ((i+1):(dim(CPMatrix)[1]))) {
      pair_count <- pair_count + 1
      
      # Get non-zero variables from both parents
      parent1_vars <- CPMatrix[i,][CPMatrix[i,] != 0]
      parent2_vars <- CPMatrix[j,][CPMatrix[j,] != 0]
      crossind <- c(parent1_vars, parent2_vars)
      
      # Count occurrences
      counts <- table(crossind)
      unique_numbers <- as.numeric(names(counts[counts == 1]))
      common_numbers <- as.numeric(names(counts[counts > 1]))
      
      # Separate main and interaction effects
      bind <- c(common_numbers, unique_numbers)
      bind_main <- bind[bind <= p]
      bind_inter <- bind[bind > p]
      
      # Identify variables from both parents
      bothmain <- bind_main[bind_main %in% common_numbers]
      bothinter <- bind_inter[bind_inter %in% common_numbers]
      
      # Set selection probabilities - prioritize common variables
      main_selection_rates <- rep(0.3, length(bind_main))
      if (length(bothmain) > 0) {
        main_selection_rates[1:length(bothmain)] <- 0.7
      }
      
      inter_selection_rates <- rep(0.3, length(bind_inter))
      if (length(bothinter) > 0) {
        inter_selection_rates[1:length(bothinter)] <- 0.7
      }
      
      # Sample main effects
      n_main_sample <- min(mysample(c(round(r1/2,0):r1), 1), length(bind_main))
      if (length(bind_main) > 0 && n_main_sample > 0) {
        sample_bind_main <- sort(mysample(bind_main, n_main_sample,
                                          replace = FALSE, prob = main_selection_rates))
        sample_bind_main <- sort(unique(sample_bind_main))
      } else {
        sample_bind_main <- numeric(0)
      }
      
      # Sample interaction effects
      n_inter_sample <- min(mysample(c(1:r2), 1), length(bind_inter))
      if (length(bind_inter) > 0 && n_inter_sample > 0) {
        sample_bind_inter <- sort(mysample(bind_inter, n_inter_sample,
                                           replace = FALSE, prob = inter_selection_rates))
      } else {
        sample_bind_inter <- numeric(0)
      }
      
      sampled_ind <- sort(c(sample_bind_main, sample_bind_inter))
      main_ind <- sampled_ind[sampled_ind <= p]
      
      # Apply heredity constraints - optimized
      inter_idx <- sampled_ind[sampled_ind > p]
      if (length(inter_idx) > 0) {
        maps <- InterRank.ind.mat[InterRank.ind.mat[,3] %in% inter_idx, , drop = FALSE]
        
        if (nrow(maps) > 0) {
          if (heredity == "Strong"){
            # Both parents must be in main_ind
            if (nrow(maps) == 1) {
              valid_mask <- maps[1] %in% main_ind & maps[2] %in% main_ind
            } else {
              valid_mask <- maps[,1] %in% main_ind & maps[,2] %in% main_ind
            }
            bbb <- if (nrow(maps) == 1) {
              if (valid_mask) maps[3] else numeric(0)
            } else {
              maps[valid_mask, 3, drop = TRUE]
            }
          } else if (heredity == "Weak"){
            # At least one parent must be in main_ind
            if (nrow(maps) == 1) {
              valid_mask <- maps[1] %in% main_ind | maps[2] %in% main_ind
            } else {
              valid_mask <- maps[,1] %in% main_ind | maps[,2] %in% main_ind
            }
            bbb <- if (nrow(maps) == 1) {
              if (valid_mask) maps[3] else numeric(0)
            } else {
              maps[valid_mask, 3, drop = TRUE]
            }
          } else {
            # No heredity constraint
            bbb <- inter_idx
          }
        } else {
          bbb <- numeric(0)
        }
      } else {
        bbb <- numeric(0)
      }
      
      # Fill offspring
      fill <- c(main_ind, bbb)
      # Remove any NA values
      fill <- fill[!is.na(fill)]
      # Safely check all conditions
      can_fill <- (!is.null(fill) && length(fill) > 0 && 
                   !is.null(pair_count) && !is.na(pair_count) && 
                   pair_count > 0 && pair_count <= nrow(CMatrix))
      
      if (isTRUE(can_fill)) {
        # Ensure we don't exceed matrix dimensions
        n_cols_to_fill <- min(length(fill), ncol(CMatrix))
        CMatrix[pair_count, 1:n_cols_to_fill] <- fill[1:n_cols_to_fill]
      }
    }
  }
  
  CMatrix[is.na(CMatrix)] <- 0
  return(CMatrix)
}

# FUNCTION USED DURING INITIAL
sampleGpool <- function(X, EVAoutput, heredity, r1, r2){
  p <- dim(X)[2]
  n <- dim(X)[1]
  if (p < n ){
    initial.main.numbs <- EVAoutput$ranked.mainpool[1:r1]
  }else{
    initial.main.numbs <- EVAoutput$ranked.mainpool
  }
  initial.inter.pool <- EVAoutput$ranked.intermat

  initial.main.select <- initial.main.numbs[1:r1]
  if (heredity == "Strong"){
    mapinter.initial.main.select <- initial.inter.pool[initial.inter.pool[, 1] %in% initial.main.select & initial.inter.pool[, 2] %in% initial.main.select, ]
  }else if (heredity == "Weak"){
    mapinter.initial.main.select <- initial.inter.pool[initial.inter.pool[, 1] %in% initial.main.select | initial.inter.pool[, 2] %in% initial.main.select, ]
  }else if (heredity == "No"){
    mapinter.initial.main.select <- initial.inter.pool
  }
  initial.inter.numbs <- as.numeric(mapinter.initial.main.select[,3])
  top_1_percent <- ceiling(0.2*length(initial.inter.numbs))
  # top_1_percent <- ceiling(0.01*length(initial.inter.numbs))
  inter_filler_probabilities <- rep(0.1, length(initial.inter.numbs))
  inter_filler_probabilities[1:top_1_percent] <- 0.9
  initial.inter.select <- mysample(initial.inter.numbs, 1, replace = FALSE, prob = inter_filler_probabilities)

  filler.initial <- c(initial.main.select, initial.inter.select)
  return(filler.initial)
}

# FUNCTION USED TO GENERATE INTERACTION POOL FOR DIFFERENT HEREDITY CONDITIONS
inter_pool_generate <- function (a, p, heredity = "Strong", interaction.ind){
  if (heredity == "No"){
    interpooltemp <- interaction.ind
  }
  else if (heredity == "Weak"){
    # Fix matrix dimension bug: safer approach with proper empty handling
    subset1_mask <- interaction.ind[,1] %in% a
    subset2_mask <- interaction.ind[,2] %in% a

    if (!any(subset1_mask) && !any(subset2_mask)) {
      # No interactions found
      interpooltemp <- matrix(ncol = 2, nrow = 0)
      colnames(interpooltemp) <- c("V1", "V2")
    } else {
      # Get subsets with proper handling
      if (any(subset1_mask)) {
        subset1 <- interaction.ind[subset1_mask, , drop=FALSE]
        order1 <- order(stats::na.exclude(match(interaction.ind[subset1_mask, 1], a)))
        subset1 <- subset1[order1, , drop=FALSE]
      } else {
        subset1 <- matrix(ncol = 2, nrow = 0)
      }

      if (any(subset2_mask)) {
        subset2 <- interaction.ind[subset2_mask, , drop=FALSE]
        order2 <- order(stats::na.exclude(match(interaction.ind[subset2_mask, 2], a)))
        subset2 <- subset2[order2, , drop=FALSE]
      } else {
        subset2 <- matrix(ncol = 2, nrow = 0)
      }

      # Combine and remove duplicates
      if (nrow(subset1) == 0) {
        interpooltemp <- subset2[!duplicated(subset2), , drop=FALSE]
      } else if (nrow(subset2) == 0) {
        interpooltemp <- subset1[!duplicated(subset1), , drop=FALSE]
      } else {
        df <- rbind(subset1, subset2)
        interpooltemp <- df[!duplicated(df), , drop=FALSE]
      }
    }
  }
  else if (heredity == "Strong"){
    # Handle boundary case: need at least 2 main effects to create pairs
    if (length(a) < 2) {
      # Return empty interaction pool if insufficient main effects
      interpooltemp <- matrix(ncol = 2, nrow = 0)
      colnames(interpooltemp) <- c("V1", "V2")
    } else {
      interpooltemp <- t(utils::combn(sort(a),2))
    }
  }
  
  # Handle empty interaction pool case
  if (nrow(interpooltemp) == 0) {
    # Return empty result matrix with proper structure
    intercandidates.mat <- matrix(ncol = 3, nrow = 0)
    colnames(intercandidates.mat) <- c("a", "b", "idx")
    return(intercandidates.mat)
  }
  
  intercandidates.ind <- match(do.call(paste, as.data.frame(interpooltemp)), do.call(paste, as.data.frame(interaction.ind)))+ p
  intercandidates.mat <- as.matrix(cbind(interpooltemp,intercandidates.ind))
  colnames(intercandidates.mat) <- c("a", "b", "idx")
  return(intercandidates.mat)
}

# FUNCTION USED TO EVALUATE INTERACTION
interaction_evaluation <- function(a, b, X, y, heredity = "Strong", sigma = NULL, varind = NULL, interaction.ind = NULL, lambda = 10){
  p <- dim(X)[2]
  if (heredity == "Strong" | heredity == "Weak"){
    interscore <- lapply(1:length(b), function(i) {
      ABC(X, y, heredity, sigma, varind = c(a, b[i]), lambda)})
  }else if(heredity == "No"){
    if (length(b) <= choose(1000,2)){
      interscore <- lapply(1:length(b), function(i) {
        ABC(X, y, heredity, sigma, varind = c(a, b[i]), lambda)})
    }else{
      interscore <- lapply(1:length(b), function(i){
        abs(as.numeric(energy::dcor(Extract(X, varind = c(b[i]), interaction.ind), y)))
      })
    }
  }
  return(interscore)
}

# FUNCTION USED TO EVALUATE MY ENTIRE FAMILY
myfamily.eva <- function(Myfamily, r1, r2, X , y, heredity = "Strong",
                         sigma = NULL, varind = NULL, interaction.ind = NULL , lambda = 10){
  G.current <- lapply(1:nrow(Myfamily), function(i) {
    ABC(X , y, heredity, sigma, varind = c(as.numeric(Myfamily[i,])[as.numeric(Myfamily[i,]) != 0]), lambda)})
  G.current.matrix <- as.matrix(cbind(Myfamily, G.current))
  colnames(G.current.matrix) <- NULL
  if (mydim(G.current.matrix)[1] == 1){
    G.current.matrix.ordered <- G.current.matrix
  }else{
    G.current.matrix.ordered <-  G.current.matrix[order(unlist(G.current.matrix[,((r1+r2)+1)]), na.last = TRUE),]
  }
  return(G.current.matrix.ordered)
}

# FUNCION USED TO PERFORM SIMULATED ANNEALING
Saint <- function(myParent, EVAoutput, r1, r2, initial.temp = 1000, cooling.rate = 0.95, X, y,
                  heredity = "Strong", sigma = NULL, varind = NULL, interaction.ind = NULL, lambda = 10){
  p <- dim(X)[2]

  # CALCULATE CURRENT ABC FOR MYPARENT
  E.best <- lapply(1:nrow(myParent), function(i){
    ABC(X , y, heredity, sigma, varind = c(as.numeric(myParent[i,])[as.numeric(myParent[i,]) != 0]), lambda)})

  mutant <- matrix(0, nrow = mydim(myParent)[1], ncol = (r1+r2))
  for (i in 1:nrow(myParent)) {
    if (i == 1){current.temp <- initial.temp}

    # DETERMINE CURRENT MAIN
    current.main.numbs <- as.numeric(myParent[i,][which(myParent[i,]%in% 1:p)])
    current.main.size <- length(current.main.numbs)
    current.inter.numbs <- as.numeric(myParent[i,][which(myParent[i,]> p)])

    # DETERMINE WHETHER MAIN EFFECTS IS ADDED OR SUBTRACTED
    # IF 0 THEN DELETE. IF 1 THEN ADD
    subtract.or.add <- mysample(0:1,1)

    # IF 1 IS SELECTED. AND THE MAIN EFFECT IS LESS THAN R1
    # ALLOW ADD
    if (subtract.or.add == 1 & current.main.size < r1){
      not.in.current.main.numbs <- as.numeric(setdiff(EVAoutput$ranked.mainpool, current.main.numbs))

      if (length(not.in.current.main.numbs) == 0){
        new.main.numbs <- current.main.numbs
        new.main.size <- length(new.main.numbs)
      }else{
        new.sample.main.numbers <- mysample(not.in.current.main.numbs, 1, replace = FALSE)
        new.main.numbs <- c(current.main.numbs, new.sample.main.numbers)
        new.main.size <- length(new.main.numbs)
      }

      # IF 1 IS SELECTED. AND THE MAIN EFFECT IS EQUAL TO R1
      # FORCE DELETE
    }else if (subtract.or.add == 1 & current.main.size == r1){
      remove.numbs <- mysample(current.main.numbs, 1, replace = FALSE)
      new.main.numbs <- as.numeric(setdiff(current.main.numbs, remove.numbs))
      new.main.size <- length(new.main.numbs)

      # IF 0 IS SELECED. AND THE MAIN EFFECT WILL BE GREATER THAN ZERO
      # ALLOW DELETE
    }else if (subtract.or.add == 0 & ((current.main.size - 1) >= 1)){
      remove.numbs <- mysample(current.main.numbs, 1, replace = FALSE)
      new.main.numbs <- as.numeric(setdiff(current.main.numbs, remove.numbs))
      new.main.size <- length(new.main.numbs)

      # IF 0 IS SELECTED. AND THE MAIN EFFECT WILL BE LOWER TO ZERO
      # ALLOW DELETE FOR NO HEREDITY
      # FORCE ADD FOR STRONG AND WEAK HEREDITY
    } else if (subtract.or.add == 0 & ((current.main.size - 1) < 1)){
      if (heredity == "No"){
        new.main.numbs <- NULL
      }else{
        not.in.current.main.numbs <- as.numeric(setdiff(EVAoutput$ranked.mainpool, current.main.numbs))
        new.sample.main.numbers <- mysample(not.in.current.main.numbs, 1, replace = FALSE)
        new.main.numbs <- c(current.main.numbs, new.sample.main.numbers)
        new.main.size <- length(new.main.numbs)
      }
    }

    # CALCULATE CURRENT ABC FOR NEW MAIN.
    # THIS IS THE ABC SCORE AFTER NEW MAIN NUMBS ADDED
    E.current.m <- ABC(X, y, heredity, sigma, varind = c(new.main.numbs, current.inter.numbs), lambda)

    # IF E CURRENT M < E BEST THEN ACCEPT CURRENT MAIN
    # IF E CURRENT M IS NOT BETTER. THEN DO SA
    if (as.numeric(E.current.m) < as.numeric(E.best[[i]])){
      best.current.main.numbs <- new.main.numbs
    }else{
      prob.value <- exp((as.numeric(E.best[[i]]) - as.numeric(E.current.m))/current.temp)
      rand.value <- runif(1)
      if (prob.value <= rand.value) {best.current.main.numbs <- new.main.numbs}else{
        best.current.main.numbs <- current.main.numbs
      }
    }

    # GENERATE NEW INTERACTION POOL. AND CHECK WHETHER IF HEREDITY CONDITION IS
    # SATISFIED FOR CURRENT INTERACTION
    new.interaction.pool <- inter_pool_generate(best.current.main.numbs, p, heredity, interaction.ind)
    new.interaction.pool.idx <- as.numeric(new.interaction.pool[,3])

    if (heredity == "Strong" | heredity == "Weak"){
      current.inter.numbs.h <- intersect(current.inter.numbs, as.numeric(new.interaction.pool[,3]))
      if (length(current.inter.numbs.h) > 0){
        current.inter.numbs <- current.inter.numbs.h
      }else{
        current.inter.numbs <- NULL
      }
    }else{
      current.inter.numbs <- current.inter.numbs
    }
    current.inter.size <- length(current.inter.numbs)

    # DETERMINE WHETHER INTERACTION EFFECTS IS ADDED OR SUBTRACTED
    # IF 0 THEN DELETE. IF 1 THEN ADD
    inter.subtract.or.add <- mysample(0:1,1)

    # IF CURRENT INTER NUMBER IS NULL
    if (is.null(current.inter.numbs)){
      new.inter.numbs <- mysample(new.interaction.pool.idx, 1, replace = FALSE)

      # IF CURRENT INTER NUMBER IS NOT NULL. DO SAMPLE
    }else{

      # IF 1 IS SELECTED. AND THE INTERACTION EFFECT IS LESS THAN R2
      # ALLOW ADD
      if (inter.subtract.or.add == 1 & current.inter.size < r2){
        not.in.current.inter.numbs <- setdiff(new.interaction.pool.idx, current.inter.numbs)
        if (length(not.in.current.inter.numbs) == 0){
          new.inter.numbs <- current.inter.numbs
          new.inter.size <- length(new.inter.numbs)
        }else{
          new.sample.inter.numbers <- mysample(not.in.current.inter.numbs, 1, replace = FALSE)
          new.inter.numbs <- c(current.inter.numbs, new.sample.inter.numbers)
          new.inter.size <- length(new.inter.numbs)
        }

        # IF 1 IS SELECTED. AND THE INTERACTION EFFECT IS EQUAL TO R2
        # FORCE DELETE
      }else if (inter.subtract.or.add == 1 & current.inter.size == r2){
        remove.numbs <- mysample(current.inter.numbs, 1, replace = FALSE)
        new.inter.numbs <- setdiff(current.inter.numbs, remove.numbs)
        new.inter.size <- length(new.inter.numbs)

        # IF 0 IS SELECED. AND THE INTERACTION EFFECT WILL BE GREATER THAN ZERO
        # ALLOW DELETE
      }else if (inter.subtract.or.add == 0 & ((current.inter.size - 1) >= 1)){
        remove.numbs <- mysample(current.inter.numbs, 1, replace = FALSE)
        new.inter.numbs <- setdiff(current.inter.numbs, remove.numbs)
        new.inter.size <- length(new.inter.numbs)

        # IF 0 IS SELECTED. AND THE INTERACTION EFFECT WILL BE LOWER TO ZERO
        # FORCE ADD
      } else if (inter.subtract.or.add == 0 & ((current.inter.size - 1) < 1)){
        not.in.current.inter.numbs <- setdiff(new.interaction.pool.idx, current.inter.numbs)
        if (length(not.in.current.inter.numbs) == 0){
          new.inter.numbs <- current.inter.numbs
          new.inter.size <- length(new.inter.numbs)
        }else{
          new.sample.inter.numbers <- mysample(not.in.current.inter.numbs, 1, replace = FALSE)
          new.inter.numbs <- c(current.inter.numbs, new.sample.inter.numbers)
          new.inter.size <- length(new.inter.numbs)
        }
      }
    }

    # CALCULATE CURRENT ABC FOR NEW MAIN AND NEW INTER.
    # THIS IS THE ABC SCORE AFTER NEW MAIN NUMB AND NEW INTER NUMB ADDED
    E.current <- ABC(X, y, heredity, sigma, varind = c(best.current.main.numbs, new.inter.numbs), lambda)

    # IF E CURRENT < E.BEST THEN ACCEPT CURRENT INTER
    # IF E CURRENT IS NOT BETTER. THEN SA
    if (as.numeric(E.current) < as.numeric(E.best[[i]])){
      best.current.inter.numbs <- new.inter.numbs
    }else{
      prob.value <- exp((as.numeric(E.best[[i]]) - as.numeric(E.current))/current.temp)
      rand.value <- runif(1)
      if (prob.value <= rand.value) {best.current.inter.numbs <- new.inter.numbs}else{
        best.current.inter.numbs <- current.inter.numbs
      }
    }

    filler <- c(best.current.main.numbs, best.current.inter.numbs)
    mutant[i, 1:length(filler)] <- filler
    current.temp <- current.temp*cooling.rate
  }
  return(mutant)
}


