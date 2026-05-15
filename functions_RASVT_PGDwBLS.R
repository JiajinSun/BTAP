### file of functions for the algorithms and simulations
### the file contains 3 parts: 
    ### main algorithm functions, utility functions for main algorithms, and other functions used in simulation

### Working directory requirement for scripts using this file:
### Set working directory to BTAP before running scripts that source this file.
### Example: getwd() should return the BTAP folder path.


### Main algorithm functions:
### (i) RASVT (range adaptive singular value thresholding): initialization 
### (ii) PGDwBLS (projected gradient descent with backtracking line search): 
    ### fixed-step / line-search GD routines with optional high-precision (Rmpfr) checks.

### Function: RASVT
### Purpose:
###   Construct the RA-SVT/USVT-style initializer (Z0, alpha0) from observed matrix A.
### Inputs:
###   A      : n x n numeric matrix, observed adjacency matrix.
###   k      : integer scalar, latent dimension.
###   family : character in {'poisson','bernoulli','gaussian'}.
###   gamma_n: numeric scalar in (0,1), trimming quantile level.
### Output:
###   List with fields:
###     Z0    : n x k numeric matrix, initial estimate of latent positions.
###     alpha0: numeric vector (length n), initial estimate of degree/intercept parameters.
RASVT <- function(A,k,family, gamma_n=NULL){
  # get_link_functions() returns g and g^{-1} for the chosen distribution family (g is mean function)
  link_funs <- get_link_functions(family)
  g_inverse <- link_funs$g_inverse
  
  ### set parameters
  n <- nrow(A)
  if(is.null(gamma_n)){
    gamma_n <- 0.1 * n^(-1/(k+4))
  }

  if(family == "bernoulli" || family == "poisson")
    tau <- sqrt(log(n) *sum(A)/n)
  if(family == "gaussian")
    tau <- sqrt(log(n) * n) 
  svdA <- svd(A)
  P_hat <- svdA$u %*% diag(svdA$d * I(svdA$d > tau)) %*% t(svdA$v)
  
  q_low <- as.numeric(quantile(P_hat, gamma_n, names = FALSE, type = 8))
  q_up <- as.numeric(quantile(P_hat, 1 - gamma_n, names = FALSE, type = 8))
  eps <- 1e-8
  
  ### project into support of g_inverse
  if(family == "bernoulli"){
    Lower <- max(q_low, eps)
    Upper <- min(q_up, 1 - eps)
  } else if(family == "poisson"){
    Lower <- max(q_low, eps)
    Upper <- max(q_up, Lower + eps)
  } else if(family == "gaussian"){
    Lower <- q_low
    Upper <- q_up
  } else {
    stop("Unsupported family in RASVT.")
  }
  if(!(Lower < Upper))
    stop("Invalid projection interval: Lower must be strictly smaller than Upper.")
  
  P_hat_proj <- pmin(pmax(P_hat, Lower), Upper)
  Theta_hat <- g_inverse(P_hat_proj)
  
  alpha0 <- solve(n * diag(rep(1,n)) + rep(1,n) %*% t(rep(1,n))) %*% Theta_hat %*% rep(1,n)
  # J <- diag(rep(1,n)) - rep(1,n) %*% t(rep(1,n)) / n
  R <- Theta_hat - (alpha0 %*% t(rep(1,n)) + rep(1,n) %*% t(alpha0))
  R <- (R + t(R))/2
  eg <- eigen(R, symmetric = TRUE)
  vals <- eg$values
  vecs <- eg$vectors
  vals_pos <- vals * (vals > 0)
  # R <- vecs %*% diag(vals_pos, n, n) %*% t(vecs)
  zvals <- pmax(vals_pos[1:k], 0)
  Z0 <- vecs %*% rbind(diag(sqrt(zvals), nrow = k, ncol = k), matrix(0, n-k, k))
  
  return(list(Z0 = Z0, alpha0 = alpha0))
}


### Function: RASVT_PGDwBLS
### Purpose:
###   Run one model fit (fixed step or line-search GD) from RASVT initialization and
###   return convergence/runtime/backtracking diagnostics and final score statistics.
### Inputs:
###   A             : n x n numeric matrix, observed adjacency matrix.
###   k             : integer scalar, latent dimension.
###   family        : character in {'poisson','bernoulli','gaussian'}.
###   eta_init      : numeric scalar, initial step size for fixed or line-search GD.
###   method        : character in {'fixed','linesearch'}.
###   eps_grad      : numeric scalar, stopping threshold on max abs Y-gradient entry.
###   max_gd_steps  : positive integer, maximum GD steps.
###   beta_ls       : numeric scalar in (0,1), backtracking shrink factor.
###   c_ls          : numeric scalar, line-search constant C_ls.
###   Rprime_line   : nonnegative integer, max backtracking reductions per GD step.
### Output:
###   If fitting succeeds, returns a list with:
###     error            : FALSE
###     Z_rasvt_raw      : n x k matrix, raw RASVT initializer before centering
###     Z_rasvt          : n x k matrix, centered RASVT initializer
###     alpha_rasvt      : length-n vector, RASVT alpha initializer
###     Z                : n x k matrix, final GD estimate of latent positions
###     alpha            : length-n vector, final GD estimate of alpha
###     num_of_steps     : integer, number of GD iterations executed
###     score_Z_max_abs  : scalar, max(abs(grad_Z)) at final iterate
###     score_alpha_max_abs: scalar, max(abs(grad_alpha)) at final iterate
###     score_Y_max_abs  : scalar, max over concatenated Y-gradient entries
###     score_Y_max_which: length-2 integer vector, index of max(abs(grad_Y))
###     score_Y_abs_mean : scalar, mean(abs(grad_Y)) at final iterate
###     converged        : 1 if stop criterion reached, else 0
###     exploded         : 1 if non-finite score encountered, else 0
###     bt_mean          : mean backtracking reductions per GD step
###     bt_max           : max backtracking reductions in any GD step
###     bt_frac0         : fraction of GD steps with zero backtracking
###     bt_frac_hit_cap  : fraction of GD steps hitting max_backtrack
###     ls_fail_frac     : fraction of GD steps where line-search failed
###     elapsed          : runtime in seconds
###   If fitting fails, returns:
###     list(error = TRUE, message = <error text>, elapsed = <seconds>)
RASVT_PGDwBLS <- function(A, k, family, eta_init, method,
                          eps_grad=1e-2, max_gd_steps=2e3, 
                          beta_ls=0.5, c_ls=1, Rprime_line){
  max_backtrack <- if(method == "linesearch") Rprime_line else 0L
  n <- nrow(A)
  one <- rep(1, n)
  J <- diag(n) - one %*% t(one) / n
  # get_link_functions() supplies the mean-link g used to form P = g(Theta).
  g <- get_link_functions(family)$g
  
  t0 <- proc.time()[["elapsed"]]
  result <- tryCatch({
    # RASVT() computes spectral initialization (Z0, alpha0) before GD starts.
    init <- RASVT(A, k, family = family)
    Z0_rasvt_raw <- init$Z0
    alpha0_rasvt <- as.vector(init$alpha0)
    Z0 <- J %*% Z0_rasvt_raw
    alpha0 <- alpha0_rasvt
    
    Theta0 <- alpha0 %*% t(one) + one %*% t(alpha0) + Z0 %*% t(Z0)
    P <- g(Theta0)
    M <- (A - P) - diag(diag(A - P))
    grad0.Z <- M %*% Z0
    grad0.a <- M %*% one
    
    bt_history <- integer(0)          # backtracking count per GD step
    line_ok_history <- logical(0)     # whether line-search accepted a step per GD step
    t <- 1L                           # GD iteration number
    score_z_now <- suppressWarnings(max(abs(grad0.Z)))
    score_alpha_now <- suppressWarnings(max(abs(grad0.a)))
    score_Y_now <- suppressWarnings(max(score_z_now, score_alpha_now))
    
    while((t < max_gd_steps) && is.finite(score_Y_now) && (score_Y_now > eps_grad)){
      t <- t + 1L
      
      d.Z <- J %*% grad0.Z
      d.a <- grad0.a
      eta_ls <- eta_init
      
      Theta0 <- alpha0 %*% t(one) + one %*% t(alpha0) + Z0 %*% t(Z0)
      
      if(max_backtrack > 0){
        # negloglik_term() returns entrywise negative log-likelihood at current Theta0.
        NLL0 <- negloglik_term(Theta0, A, family)
        diag(NLL0) <- 0
        # Line-search uses two expansions:
        # (a) Global update along d(Y) = (d.Z, d.a):
        #     Theta_full(eta) = Theta0 + eta * B_full + eta * t(B_full) + eta^2 * C_full,
        #     B_full = d.a %*% t(1_n) + d.Z %*% t(Z0),   C_full = d.Z %*% t(d.Z).
        #     Global check:
        #       [L(Theta_full(eta)) - L(Theta0)] - c_ls * n * eta^2 * d_inner_global <= 0.
        # (b) Individual row-wise update (for all rows i):
        #     Theta_row(eta) = Theta0 + eta * B_row,
        #     where row i of Theta_row(eta) updates Theta0_i for evaluating L_i.
        #       B_row = grad0.a %*% t(1_n) + grad0.Z %*% t(Z0).
        #     Individual checks:
        #       [L_i(Theta_row(eta)) - L_i(Theta0)] - c_ls * n * eta^2 * d_inner_i[i] <= 0,
        #       for each i = 1,...,n.
        B_full <- d.a %*% t(one) + d.Z %*% t(Z0)
        C_full <- d.Z %*% t(d.Z)
        B_row <- grad0.a %*% t(one) + grad0.Z %*% t(Z0)
        d_inner_global <- sum((-grad0.Z) * d.Z) + sum((-grad0.a) * d.a)
        d_inner_i <- -rowSums(grad0.Z^2) - grad0.a^2   ## n-dim vector
        
        # check_linesearch_conditions() evaluates the global and individual checks above.
        line_search <- check_linesearch_conditions(
          eta_ls, Theta0, NLL0, B_full, C_full, B_row, A, family, c_ls, n,
          d_inner_global, d_inner_i
        )
        line_search_ok <- line_search$ok
        bt <- 0L
        while((!line_search_ok) && (bt < max_backtrack)){
          eta_ls <- beta_ls * eta_ls    ## backtrack
          # Re-check acceptance rules at shrunken step size beta_ls * eta_ls.
          line_search <- check_linesearch_conditions(
            eta_ls, Theta0, NLL0, B_full, C_full, B_row, A, family, c_ls, n,
            d_inner_global, d_inner_i
          )
          line_search_ok <- line_search$ok
          bt <- bt + 1L
        }
        bt_history <- c(bt_history, bt)
        line_ok_history <- c(line_ok_history, line_search_ok)
      } else {        ### if max_backtrack = 0, i.e., fixed step GD
        bt_history <- c(bt_history, 0L)
        line_ok_history <- c(line_ok_history, TRUE)
      }
      
      ### GD
      Z0 <- Z0 + eta_ls * d.Z
      alpha0 <- alpha0 + eta_ls * d.a
      
      ### calculate score to use in the next step
      Theta0 <- alpha0 %*% t(one) + one %*% t(alpha0) + Z0 %*% t(Z0)
      P <- g(Theta0)
      M <- (A - P) - diag(diag(A - P))
      grad0.Z <- M %*% Z0
      grad0.a <- M %*% one
      score_z_now <- suppressWarnings(max(abs(grad0.Z)))
      score_alpha_now <- suppressWarnings(max(abs(grad0.a)))
      score_Y_now <- suppressWarnings(max(score_z_now, score_alpha_now))
    }

    ### record last 1 step score max
    score_Z_max_abs <- score_z_now
    score_alpha_max_abs <- score_alpha_now
    score_Y_max_abs <- score_Y_now
    if (is.finite(score_Y_now)) {
      abs_gradY <- cbind(abs(grad0.Z), abs(grad0.a))
      score_Y_max_abs_ind <- which.max(abs_gradY)
      score_Y_max_which <- c(
        score_Y_max_abs_ind - n * (ceiling(score_Y_max_abs_ind / n) - 1),
        ceiling(score_Y_max_abs_ind / n)
      )
      score_Y_abs_mean <- mean(abs_gradY)
    } else {
      score_Y_max_which <- c(NA_integer_, NA_integer_)
      score_Y_abs_mean <- NA_real_
    }
    
    list(
      error = FALSE,
      Z_rasvt_raw = Z0_rasvt_raw,
      Z_rasvt = J %*% Z0_rasvt_raw,
      alpha_rasvt = alpha0_rasvt,
      Z = Z0,
      alpha = as.vector(alpha0),
      num_of_steps = t,
      score_Z_max_abs = score_Z_max_abs,
      score_alpha_max_abs = score_alpha_max_abs,
      score_Y_max_abs = score_Y_max_abs,
      score_Y_max_which = score_Y_max_which,
      score_Y_abs_mean = score_Y_abs_mean,
      converged = as.integer(is.finite(score_Y_now) && (score_Y_now <= eps_grad)),
      exploded = as.integer(!is.finite(score_Y_now)),
      bt_mean = if(length(bt_history) > 0) mean(bt_history) else NA_real_,
      bt_max = if(length(bt_history) > 0) max(bt_history) else NA_real_,
      bt_frac0 = if(length(bt_history) > 0) mean(bt_history == 0) else NA_real_,
      bt_frac_hit_cap = if(length(bt_history) > 0) mean(bt_history >= max_backtrack) else NA_real_,
      ls_fail_frac = if(length(line_ok_history) > 0) mean(!line_ok_history) else NA_real_
    )
  }, error = function(e){
    list(error = TRUE, message = conditionMessage(e))
  })
  
  elapsed <- proc.time()[["elapsed"]] - t0
  if(isTRUE(result$error)){
    return(list(
      error = TRUE,
      message = result$message,
      elapsed = elapsed
    ))
  }
  result$elapsed <- elapsed
  result
}






### Utility function set for the two main algorithm functions, including functions to
### get link function, calculate negative log-likelihood in ordinary precision and high precision for different distributions,
### and check line search conditions in ordinary precision and high precision

### Function: get_link_functions
### Purpose:
###   Return canonical mean-link g(·) and inverse-link g^{-1}(·) for the selected family.
### Inputs:
###   family: character in {'poisson','bernoulli','gaussian'}.
### Output:
###   List with fields:
###     g        : function mapping Theta -> mean parameter.
###     g_inverse: function mapping mean parameter -> Theta.
get_link_functions <- function(family){
  if(family == 'poisson')
    return(list(g = function(x) exp(x), g_inverse = function(x) log(x)))
  if(family == 'bernoulli')
    return(list(g = function(x) 1/(1+exp(-x)), g_inverse = function(x) log(x/(1-x))))
  if(family == 'gaussian')
    return(list(g = function(x) x, g_inverse = function(x) x))
  stop("Unsupported family in get_link_functions. Use 'poisson', 'bernoulli', or 'gaussian'.")
}


### Function: negloglik_term
### Purpose:
###   Return elementwise negative log-likelihood contributions for one exponential-family choice.
### Inputs:
###   theta : numeric scalar/vector/matrix, evaluated parameter(s) Theta_ij.
###   obs   : numeric scalar/vector/matrix, same shape as theta, observed A_ij.
###   family: character in {'poisson','bernoulli','gaussian'}.
### Output:
###   Numeric object with same shape as theta/obs; each entry is -ell(theta_ij; obs_ij).
negloglik_term <- function(theta, obs, family){
  if(family == 'bernoulli'){
    # log1p(x) computes log(1+x) stably; this is a stable softplus implementation.
    softplus <- ifelse(theta > 0, theta + log1p(exp(-theta)), log1p(exp(theta)))
    return(softplus - obs * theta)
  }
  if(family == 'poisson')
    return(exp(theta) - obs * theta)
  if(family == 'gaussian')
    return(0.5 * (obs - theta)^2)
  stop("Unsupported family in negloglik_term. Use 'poisson', 'bernoulli', or 'gaussian'.")
}

### Function: nll_term_mpfr_vec
### Purpose:
###   High-precision (mpfr) version of elementwise negative log-likelihood; 
###   this function is used in the mpfr version of checking global/individual linesearch conditions
### Inputs:
###   theta_mpfr: mpfr vector (length m), evaluated parameter(s) values. m would be n^2 (for global) or n (for individual) in our code
###   obs       : numeric vector (length m), observed A_{ij} values aligned with theta_mpfr.
###   family    : character in {'poisson','bernoulli','gaussian'}.
###   precBits  : integer scalar, mpfr precision in bits.
### Output:
###   mpfr vector (length m), elementwise negative log-likelihood terms.
nll_term_mpfr_vec <- function(theta_mpfr, obs, family, precBits = 192){
  obs_mpfr <- Rmpfr::mpfr(obs, precBits = precBits)
  if(family == "bernoulli"){
    softplus <- Rmpfr::mpfr(rep(0, length(theta_mpfr)), precBits = precBits)
    pos <- as.numeric(theta_mpfr) > 0
    if(any(pos)){
      softplus[pos] <- theta_mpfr[pos] + log1p(exp(-theta_mpfr[pos]))
    }
    if(any(!pos)){
      softplus[!pos] <- log1p(exp(theta_mpfr[!pos]))
    }
    return(softplus - obs_mpfr * theta_mpfr)
  }
  if(family == "poisson"){
    return(exp(theta_mpfr) - obs_mpfr * theta_mpfr)
  }
  if(family == "gaussian"){
    half <- Rmpfr::mpfr(0.5, precBits = precBits)
    return(half * (obs_mpfr - theta_mpfr)^2)
  }
  stop("Unsupported family in nll_term_mpfr_vec.")
}


### Function: check_linesearch_conditions
### Purpose:
###   Evaluate the proposed global and rowwise individual conditions for one
###   candidate step eta_ls, with optional mpfr recomputation near numerical boundaries.
    # Line-search uses two expansions:
    # (a) Global update along d(Y) = (d.Z, d.a):
    #     Theta_full(eta) = Theta0 + eta * B_full + eta * t(B_full) + eta^2 * C_full,
    #     B_full = d.a %*% t(1_n) + d.Z %*% t(Z0),   C_full = d.Z %*% t(d.Z).
    #     Global check:
    #       [L(Theta_full(eta)) - L(Theta0)] - c_ls * n * eta^2 * d_inner_global <= 0.
    # (b) Individual row-wise update (for all rows i):
    #     Theta_row(eta) = Theta0 + eta * B_row,
    #     where row i of Theta_row(eta) updates Theta0_i for evaluating L_i.
    #       B_row = grad0.a %*% t(1_n) + grad0.Z %*% t(Z0).
    #     Individual checks:
    #       [L_i(Theta_row(eta)) - L_i(Theta0)] - c_ls * n * eta^2 * d_inner_i[i] <= 0,
    #       for each i = 1,...,n.
### Inputs:
###   eta_ls        : numeric scalar, candidate step size.
###   Theta0        : n x n numeric matrix, current Theta matrix.
###   NLL0          : n x n numeric matrix, current elementwise NLL matrix.
###   B_full        : n x n numeric matrix, first-order full-update term.
###   C_full        : n x n numeric matrix, second-order full-update term.
###   B_row         : n x n numeric matrix, first-order rowwise-update term.
###   A             : n x n numeric matrix, observed adjacency matrix.
###   family        : character in {'poisson','bernoulli','gaussian'}.
###   c_ls          : numeric scalar, line-search constant C_ls.
###   n             : integer scalar, number of nodes.
###   d_inner_global: numeric scalar, full-update inner-product term.
###   d_inner_i     : numeric vector (length n), rowwise inner-product terms.
### Output:
###   List with fields:
###     ok           : logical scalar, TRUE iff both checks pass.
###     global_ok    : logical scalar, global-condition status.
###     individual_ok: logical scalar, rowwise individual-condition status.
check_linesearch_conditions <- function(eta_ls, Theta0, NLL0,
                                        B_full, C_full, B_row, A,
                                        family, c_ls = 1, n,
                                        d_inner_global, d_inner_i){
  trigger_eps <- n * 1e-16   # near-zero margin threshold; if |LHS| < trigger_eps, recompute that check in mpfr precision
  # Global condition: full update Y + eta d(Y)
  Theta_full <- Theta0 + eta_ls * B_full + eta_ls * t(B_full) + (eta_ls^2) * C_full
  # negloglik_term() computes the full-matrix NLL at Theta_full.
  NLL_full <- negloglik_term(Theta_full, A, family)
  diag(NLL_full) <- 0
  deltaL <- sum((NLL_full - NLL0)[upper.tri(NLL0)])
  lhs_global <- deltaL - c_ls * n * eta_ls^2 * d_inner_global  # LHS of the global condition
  global_ok <- (lhs_global <= 0)  # global condition satisfied or not
  if(abs(lhs_global) < trigger_eps){
    # recompute_cond_global_mpfr_stablediff() re-evaluates near-zero global margin in mpfr.
    rec_global <- recompute_cond_global_mpfr_stablediff(
      eta_ls = eta_ls, Theta0 = Theta0, NLL0 = NLL0, B_full = B_full, C_full = C_full,
      A = A, family = family, c_ls = c_ls, n = n, d_inner_global = d_inner_global,
      precBits = 192
    )
    if(isTRUE(rec_global$available)){
      lhs_global <- rec_global$lhs_global  # LHS of the global condition (after mpfr recompute)
      global_ok <- rec_global$global_ok  # global condition satisfied or not (after mpfr recompute)
    }
  }
  if(!global_ok)
    return(list(ok = FALSE, global_ok = FALSE, individual_ok = FALSE))
  
  # Individual condition: rowwise updates Y + eta d_i(Y), all i evaluated in vectorized form
  Theta_row <- Theta0 + eta_ls * B_row
  # negloglik_term() computes rowwise-updated NLL matrix at Theta_row.
  NLL_row <- negloglik_term(Theta_row, A, family)
  diag(NLL_row) <- 0
  deltaLi <- rowSums(NLL_row - NLL0)
  lhs_individual_vec <- deltaLi - c_ls * n * eta_ls^2 * d_inner_i  # LHS of the individual condition
  individual_ok <- (max(lhs_individual_vec) <= 0)  # individual condition satisfied or not
  
  # Targeted refinement only for near-zero rowwise margins:
  # directly recompute these rows with mpfr
  near_i <- which(abs(lhs_individual_vec) < trigger_eps)
  if(length(near_i) > 0){
    for(i in near_i){
      # recompute_cond_individual_mpfr_stablediff() rechecks row-i near-zero margin in mpfr.
      lhs_i_mpfr <- recompute_cond_individual_mpfr_stablediff(
        i = i, eta_ls = eta_ls, Theta0 = Theta0, B_row = B_row,
        A = A, family = family, c_ls = c_ls, n = n, d_inner_i = d_inner_i,
        precBits = 192
      )
      if(is.finite(lhs_i_mpfr)){
        lhs_individual_vec[i] <- lhs_i_mpfr  # LHS of the individual condition for row i (after mpfr recompute)
      }
    }
    individual_ok <- (max(lhs_individual_vec) <= 0) # individual condition satisfied or not (after mpfr recompute)
  }
  return(list(
    ok = (global_ok && individual_ok),
    global_ok = global_ok,
    individual_ok = individual_ok
  ))
}


### Function: recompute_cond_global_mpfr_stablediff
### Purpose:
###   Recompute the global line-search condition left-hand side in mpfr arithmetic:
###     [L(Y+eta d) - L(Y)] - c_ls * n * eta^2 * <grad, d>.
###   Theta expansion used by this function:
###     Theta_full(eta) = Theta0 + eta * B_full + eta * t(B_full) + eta^2 * C_full,
###   with
###     B_full = d_a %*% t(1_n) + d_Z %*% t(Z0),
###     C_full = d_Z %*% t(d_Z).
### Inputs:
###   eta_ls        : numeric scalar, candidate step size eta.
###   Theta0        : n x n numeric matrix
###   NLL0          : n x n numeric matrix, current elementwise NLL matrix at Theta0.
###   B_full        : n x n numeric matrix, first-order full-update expansion term.
###   C_full        : n x n numeric matrix, second-order full-update term.
###   A             : n x n numeric matrix, observed adjacency matrix.
###   family        : character in {'poisson','bernoulli','gaussian'}.
###   c_ls          : numeric scalar, line-search constant C_ls.
###   n             : integer scalar, number of nodes.
###   d_inner_global: numeric scalar, inner product term < -grad, d > for full update.
###   precBits      : integer scalar, mpfr precision in bits.
### Output:
###   List with fields:
###     available: logical scalar, TRUE if Rmpfr is available.
###     global_ok : logical scalar, TRUE iff recomputed global condition is satisfied.
###     lhs_global: numeric scalar, recomputed global-condition LHS value.
recompute_cond_global_mpfr_stablediff <- function(eta_ls, Theta0, NLL0, B_full, C_full, A, family, c_ls, n,
                                                  d_inner_global, precBits = 192){
  if(!requireNamespace("Rmpfr", quietly = TRUE)){
    return(list(available = FALSE, global_ok = FALSE, lhs_global = NA_real_))
  }
  eta_mpfr <- Rmpfr::mpfr(eta_ls, precBits = precBits)
  Theta0_mpfr <- Rmpfr::mpfr(Theta0, precBits = precBits)
  B_full_mpfr <- Rmpfr::mpfr(B_full, precBits = precBits)
  C_full_mpfr <- Rmpfr::mpfr(C_full, precBits = precBits)
  
  Theta_full_mpfr <- Theta0_mpfr + eta_mpfr * B_full_mpfr + eta_mpfr * t(B_full_mpfr) + (eta_mpfr^2) * C_full_mpfr
  # nll_term_mpfr_vec() evaluates per-edge NLL in high precision (Rmpfr).
  NLL_full_mpfr <- matrix(
    nll_term_mpfr_vec(as.vector(Theta_full_mpfr), as.vector(A), family, precBits = precBits),
    nrow = nrow(Theta0), ncol = ncol(Theta0)
  )
  NLL0_mpfr <- Rmpfr::mpfr(NLL0, precBits = precBits)
  diag(NLL_full_mpfr) <- Rmpfr::mpfr(0, precBits = precBits)
  diag(NLL0_mpfr) <- Rmpfr::mpfr(0, precBits = precBits)
  
  deltaL_mpfr <- sum((NLL_full_mpfr - NLL0_mpfr)[upper.tri(NLL0)])
  rhs_global_mpfr <- Rmpfr::mpfr(c_ls * n * eta_ls^2 * d_inner_global, precBits = precBits)
  lhs_global <- as.numeric(deltaL_mpfr - rhs_global_mpfr)
  list(available = TRUE, global_ok = (lhs_global <= 0), lhs_global = lhs_global)
}

### Function: recompute_cond_individual_mpfr_stablediff
### Purpose:
###   Recompute rowwise individual line-search condition for one row i in mpfr arithmetic:
###     [L_i(Y+eta d_i) - L_i(Y)] - c_ls * n * eta^2 * <grad_i, d_i>.
###   Rowwise Theta expansion used by this function:
###     Theta_row(eta) = Theta0 + eta * B_row,
###   and the check for node i uses the i-th row of Theta_row.
###   In current code, B_row is formed as:
###     B_row = grad_a %*% t(1_n) + grad_Z %*% t(Z0).
### Inputs:
###   i        : integer scalar in {1,...,n}, row/node index.
###   eta_ls   : numeric scalar, candidate step size eta.
###   Theta0   : n x n numeric matrix
###   B_row    : n x n numeric matrix, rowwise first-order update term.
###   A        : n x n numeric matrix, observed adjacency matrix.
###   family   : character in {'poisson','bernoulli','gaussian'}.
###   c_ls     : numeric scalar, line-search constant C_ls.
###   n        : integer scalar, number of nodes.
###   d_inner_i: numeric vector (length n), rowwise inner-product terms.
###   precBits : integer scalar, mpfr precision in bits.
### Output:
###   Numeric scalar (double), recomputed LHS value for row i (NA if Rmpfr unavailable).
recompute_cond_individual_mpfr_stablediff <- function(i, eta_ls, Theta0, B_row, A, family, c_ls, n, d_inner_i, precBits = 192){
  if(!requireNamespace("Rmpfr", quietly = TRUE)){
    return(NA_real_)
  }
  eta_mpfr <- Rmpfr::mpfr(eta_ls, precBits = precBits)
  theta0_i_mpfr <- Rmpfr::mpfr(Theta0[i, ], precBits = precBits)
  brow_i_mpfr <- Rmpfr::mpfr(B_row[i, ], precBits = precBits)
  theta_new_i_mpfr <- theta0_i_mpfr + eta_mpfr * brow_i_mpfr
  
  # nll_term_mpfr_vec() recomputes row-i NLL terms with mpfr precision.
  nll_new_i <- nll_term_mpfr_vec(theta_new_i_mpfr, A[i, ], family, precBits = precBits)
  nll0_i <- nll_term_mpfr_vec(theta0_i_mpfr, A[i, ], family, precBits = precBits)
  nll_new_i[i] <- Rmpfr::mpfr(0, precBits = precBits)
  nll0_i[i] <- Rmpfr::mpfr(0, precBits = precBits)
  
  deltaLi_i <- sum(nll_new_i - nll0_i)
  rhs_i_mpfr <- Rmpfr::mpfr(c_ls * n * eta_ls^2 * d_inner_i[i], precBits = precBits)
  as.numeric(deltaLi_i - rhs_i_mpfr)
}










##### other functions used in the simulations, including 
###   `eta_cornodiv6_from_init`: calculate eta_0 from initial estimate;
###   `DLii`: compute the i-th diagonal expected negative Hessian block
###   `build_asymptotic_record`: extract (Z,alpha) true, (Z,alpha) estimate, and calculate t(z_{q,11}) from the GD estimate

### Function: eta_cornodiv6_from_init
### Purpose:
###   Compute eta_0 in the corollary-style formula without the denominator 6:
###     eta_0 = [1 / max_{i<j}{1 - ell''(Theta^0_ij)}] * min{1/||Z0||_op^2, 1/n},
###   where (Z0, alpha0) comes from RASVT initialization.
### Inputs:
###   A     : n x n numeric matrix, observed adjacency matrix.
###   k     : integer scalar, latent dimension.
###   family: character in {'poisson','bernoulli','gaussian'}.
### Output:
###   Numeric scalar eta_0.
eta_cornodiv6_from_init <- function(A, k, family){
  n <- nrow(A)
  J <- diag(rep(1, n)) - rep(1, n) %*% t(rep(1, n)) / n
  # RASVT() provides initialization used to evaluate corollary-style eta_0.
  init <- RASVT(A, k, family = family)
  Z0 <- J %*% init$Z0
  alpha0 <- as.vector(init$alpha0)
  Theta0 <- alpha0 %*% t(rep(1, n)) + rep(1, n) %*% t(alpha0) + Z0 %*% t(Z0)

  if(family == "bernoulli"){
    ell_pp <- -exp(Theta0) / (1 + exp(Theta0))^2
  } else if(family == "poisson"){
    ell_pp <- -exp(Theta0)
  } else if(family == "gaussian"){
    ell_pp <- -matrix(1, n, n)
  } else {
    stop("Unsupported family in eta_cornodiv6_from_init.")
  }

  curvature_max <- max((1 - ell_pp)[upper.tri(ell_pp)])
  z_op <- svd(Z0, nu = 0, nv = 0)$d[1]
  if(!is.finite(z_op) || z_op <= 0) z_op <- sqrt(n)
  (1 / curvature_max) * min(1 / (z_op^2), 1 / n)
}



### Function: DLii
### Purpose:
###   Compute node-specific expected negative Hessian block wrt y_i = (z_i, alpha_i)
### Inputs:
###   Z     : n x k numeric matrix, latent positions.
###   alpha : numeric vector (length n), degree/intercept parameters.
###   family: character in {'poisson','bernoulli','gaussian'}.
###   i     : index of the node for which to compute the Hessian
### Output:
###   Numeric matrix of size (k+1) x (k+1).
DLii <- function(Z,alpha,family, i){
  n <- nrow(Z)
  k <- ncol(Z) + 1
  if(family == 'poisson')
    nu_pp <- function(x) exp(x)
  if(family == 'bernoulli')
    nu_pp <- function(x) exp(x)/(1+exp(x))^2
  if(family == 'gaussian')
    nu_pp <- function(x) !is.na(x)
  mu <- nu_pp(alpha %*% t(rep(1,n)) + rep(1,n) %*% t(alpha) + Z %*% t(Z))
  
  W <- cbind(Z,rep(1,n))
  DL11 <- matrix(0,k,k)
  
  DL11 <- - mu[i,i] * W[i,] %*% t(W[i,])
  for(jj in 1:n)
    DL11 <- DL11 + mu[i,jj] * W[jj,] %*% t(W[jj,])
  
  return(DL11)
}

### Function: build_asymptotic_record
### Purpose:
###   Build per-run records of final GD estimator:
###   err_vec = c(t(Z), alpha, t(est$Z), est$alpha, t(Delta_Z), Delta_alpha,
###               score_Y_max_abs, score_Y_max_which, score_Y_abs_mean,
###               num_of_steps, c(O), t11_hat, t11_star),
###   where
###     t11_hat  = (zhat_q11 - z11) / sqrt([Sigma_1(hatY_q)^{-1}]_{1,1}),
###     t11_star = (zhat_q11 - z11) / sqrt([Sigma_1(Y*)^{-1}]_{1,1}).
### Inputs:
###   Z_true   : n x k numeric matrix.
###   alpha_true: numeric vector (length n).
###   est      : output list from RASVT_PGDwBLS (must include Z, alpha, score stats).
###   family   : character in {'poisson','bernoulli','gaussian'}.
### Output:
###   Numeric vector err_vec in the exact concatenation order above.
build_asymptotic_record <- function(Z_true, alpha_true, est, family){
  n <- nrow(Z_true)
  k <- ncol(Z_true)
  alpha_true <- as.vector(alpha_true)
  Z_est <- est$Z
  alpha_est <- as.vector(est$alpha)
  
  if(any(!is.finite(Z_est)) || any(!is.finite(alpha_est))){
    O <- diag(k)
    Np <- n * k + n
    err_vec <- c(
      c(t(Z_true)),
      c(alpha_true),
      rep(NA_real_, 2 * Np + 1 + 2 + 1 + 1),
      c(O),
      NA_real_,
      NA_real_
    )
    return(err_vec)
  }
  
  sv <- svd(t(Z_est) %*% Z_true)
  O <- sv$v %*% t(sv$u)
  Z_q <- Z_est %*% t(O)   ## aligning Zhat to Zstar
  
  Delta_Z <- Z_q - Z_true
  Delta_alpha <- alpha_est - alpha_true
  
  # Robust matrix inverse helper:
  # try solve(M) first; if singular/ill-conditioned, fallback to qr.solve(M).
  # Returns NULL if both fail.
  solve_mat <- function(M){
    out <- tryCatch(solve(M), error = function(e) NULL)
    if(is.null(out)){
      out <- tryCatch(qr.solve(M), error = function(e) NULL)
    }
    out
  }
  # DLii() builds node-wise curvature blocks; solve_mat() inverts them robustly.
  Dhat1_inv <- solve_mat(DLii(Z_q, alpha_est, family, 1))
  Dstar1_inv <- solve_mat(DLii(Z_true, alpha_true, family, 1))
  num11 <- Z_q[1, 1] - Z_true[1, 1]
  
  t11_hat <- NA_real_
  if(!is.null(Dhat1_inv)){
    den_hat <- Dhat1_inv[1, 1]
    if(is.finite(den_hat) && den_hat > 0){
      t11_hat <- num11 / sqrt(den_hat)
    }
  }
  
  t11_star <- NA_real_
  if(!is.null(Dstar1_inv)){
    den_star <- Dstar1_inv[1, 1]
    if(is.finite(den_star) && den_star > 0){
      t11_star <- num11 / sqrt(den_star)
    }
  }
  
  err_vec <- c(
    c(t(Z_true)),                        ## nk dim
    c(alpha_true),                       ## n dim
    c(t(Z_est)),                         ## nk dim
    c(alpha_est),                        ## n dim
    c(t(Delta_Z)),                       ## nk dim
    c(Delta_alpha),                      ## n dim
    est$score_Y_max_abs,                 ## 1 dim
    c(est$score_Y_max_which),            ## 2 dim
    est$score_Y_abs_mean,                ## 1 dim
    est$num_of_steps,                    ## 1 dim
    c(O),                                ## k^2 dim
    t11_hat,                             ## 1 dim
    t11_star                             ## 1 dim
  )
  err_vec
}
