################################################################################
# simulation_functions.R
#
# Functions for the simulation study of the Bayesian joint latent class causal
# model (the "proposed" model) and the two-step inverse probability weighted
# latent class estimator of Lanza et al. (the "comparator").
#
# This file only DEFINES functions and constants. Run the study with
#   simulation_code.R          proposed model + comparator (bootstrap B = 50)
#   lanza_bootstrap_B1000.R    comparator re-run with B = 1000, same datasets
#
# Notation (same as the manuscript)
#   A      exposure, 1 = active ("physical_activity")
#   L      baseline confounders: education (Low/Medium/High), sex, age
#   C      latent class: 1 = HH, 2 = LL, 3 = HL, 4 = LH (memory / executive)
#   Y      six continuous outcomes: REY_I, REY_II, MAT, AFT, FAS, STP
#   gamma  coefficients of the exposure model  P(A | L)
#   alpha  coefficients of the class model     P(C | A, L)   (class 1 = reference)
#   theta  class-specific outcome means; sigma common outcome SD
#   psi_c(a) = P(C^a = c);   ATE_c = psi_c(1) - psi_c(0)
#
# Contents
#   1. Design: constants, scenario grid, seeds
#   2. Data generation
#   3. Proposed model: JAGS fit
#   4. Label switching
#   5. Posterior predictive standardization (g-computation)
#   6. True causal effect
#   7. Comparator: two-step IPW latent class estimator
#   8. One simulated dataset, end to end
#   9. Summaries (MSE and coverage)
#  10. Utilities
################################################################################

suppressPackageStartupMessages({
  library(truncnorm)
  library(clue)
  library(flexmix)
  library(WeightIt)
  library(nnet)
  library(future.apply)
})
# rjags (and the JAGS program) is only needed by fit_one(); it is loaded there,
# so the comparator script runs without JAGS.


# 1. Design ---------------------------------------------------------------------

SIM <- list(
  n               = 2000,             # subjects per dataset
  K               = 4L,               # true (and fitted) number of classes
  sigma           = 0.5,              # common outcome SD
  target_treat    = 0.3,              # marginal P(A = 1), fixed in every scenario
  D_per_scenario  = 100L,             # datasets per scenario
  master_seed     = 2025,
  outcome_names   = c("REY_I", "REY_II", "MAT", "AFT", "FAS", "STP"),
  class_names     = c("High-High", "Low-Low", "High-Low", "Low-High")
)

# Class proportions: imbalance level -> target P(C = c)
target_pi_from_label <- function(lbl) {
  switch(lbl,
         "balanced"   = c(0.25, 0.25, 0.25, 0.25),
         "one_rare"   = c(0.30, 0.30, 0.30, 0.10),
         "two_rare"   = c(0.40, 0.40, 0.10, 0.10),
         "severe_one" = c(0.30, 0.35, 0.30, 0.05),
         stop("unknown imbalance label: ", lbl))
}

# Class separation: level -> scale delta applied to the base class means
delta_from_label <- function(lbl) {
  switch(lbl,
         "low"  = 0.6,
         "high" = 1.2,
         stop("unknown separation label: ", lbl))
}

# Confounding strength: level -> scale rho applied to the confounder coefficients
# of the exposure model (gamma) and of the class model (alpha)
conf_scales_from_label <- function(lbl) {
  switch(lbl,
         "weak"   = list(gamma_scale = 0.6, alpha_scale = 0.6),
         "strong" = list(gamma_scale = 1.6, alpha_scale = 1.6),
         stop("unknown confounding label: ", lbl))
}

scenario_tag <- function(scn) {
  paste0("n", scn$n, "_conf-", scn$confounding,
         "_imb-", scn$imbalance, "_sep-", scn$separation)
}

# Scenario grid (2 x 4 x 2 = 16 scenarios) and one seed per (scenario, dataset).
#
# Seeds: during development the study also had a third, smaller separation level
# (delta = 0.3) that is not reported. The seed vector was drawn over that
# 3-level, 24-scenario grid, and the 8 scenarios with delta = 0.3 came first.
# To reproduce exactly the datasets analysed in the paper, we draw the same
# 24 x 100 seeds and drop the first 8 x 100. The seed of dataset d in scenario r
# is therefore seeds[(r - 1) * D_per_scenario + d].
make_design <- function() {
  scenarios <- expand.grid(
    confounding = c("weak", "strong"),
    imbalance   = c("balanced", "one_rare", "two_rare", "severe_one"),
    separation  = c("low", "high"),
    n           = SIM$n,
    stringsAsFactors = FALSE
  )
  scenarios$scenario_id <- seq_len(nrow(scenarios))

  n_dropped_blocks <- 8L
  set.seed(SIM$master_seed)
  all_seeds <- sample.int(1e9, SIM$D_per_scenario * (nrow(scenarios) + n_dropped_blocks))
  seeds     <- all_seeds[-seq_len(SIM$D_per_scenario * n_dropped_blocks)]

  list(scenarios = scenarios, seeds = seeds)
}

# helpers
rcat <- function(M, prob_vec) sample.int(length(prob_vec), size = M, replace = TRUE, prob = prob_vec)

softmax_vec <- function(M) {            # row-wise softmax of a matrix
  ex <- exp(M - apply(M, 1, max))
  ex / rowSums(ex)
}

# Exposure intercept gamma0 such that mean P(A = 1) equals the target (bisection)
tune_gamma0_for_treatshare <- function(eta_no_intercept, target = 0.5, tol = 1e-8, maxit = 60) {
  f  <- function(g0) mean(plogis(eta_no_intercept + g0)) - target
  lo <- -20; hi <- 20
  flo <- f(lo); fhi <- f(hi)
  itw <- 0
  while (flo * fhi > 0 && itw < 10) {
    lo <- lo - 10; hi <- hi + 10
    flo <- f(lo);  fhi <- f(hi)
    itw <- itw + 1
  }
  for (it in 1:maxit) {
    mid  <- 0.5 * (lo + hi)
    fmid <- f(mid)
    if (abs(fmid) < tol) return(mid)
    if (flo * fmid <= 0) { hi <- mid; fhi <- fmid } else { lo <- mid; flo <- fmid }
  }
  0.5 * (lo + hi)
}

# Class intercepts alpha0 such that the marginal class proportions equal target_pi
tune_alpha0_for_classmix <- function(X_class, B, target_pi, maxit = 500, reltol = 1e-8) {
  n   <- nrow(X_class)
  Km1 <- ncol(B)
  K   <- Km1 + 1
  loss <- function(par) {
    Lin <- matrix(0, n, K)
    Lin[, 2:K] <- sweep(X_class %*% B, 2, par, FUN = "+")
    sum((colMeans(softmax_vec(Lin)) - target_pi)^2)
  }
  opt <- optim(rep(0, Km1), loss, method = "L-BFGS-B",
               control = list(maxit = maxit, reltol = reltol))
  c(0, opt$par)
}

# Class means theta_c(delta) = theta_bar + delta * (theta_c0 - theta_bar)
scale_separation <- function(theta_matrix, delta = 1) {
  mu0 <- colMeans(theta_matrix)
  D   <- sweep(theta_matrix, 2, mu0, FUN = "-")
  theta_scaled <- sweep(delta * D, 2, mu0, FUN = "+")
  rownames(theta_scaled) <- rownames(theta_matrix)
  colnames(theta_scaled) <- colnames(theta_matrix)
  theta_scaled
}


# 2. Data generation ---------------------------------------------------------------

# Simulate one dataset.
#   delta        class separation scale
#   gamma_scale  rho applied to the confounder coefficients of the exposure model
#   alpha_scale  rho applied to the confounder coefficients of the class model
# Returns the data and the quantities needed to evaluate the truth.
simulate_one <- function(seed, n = 800,
                         target_treat = 0.5,
                         target_pi = c(0.25, 0.25, 0.25, 0.25),
                         delta = 1.0,
                         sigma = 0.5,
                         gamma_scale = 1.0,
                         alpha_scale = 1.0) {
  set.seed(seed)

  ## baseline confounders L
  edu_probs       <- c(0.2, 0.5, 0.3)
  X_edu           <- sample(1:3, size = n, replace = TRUE, prob = edu_probs)
  education_level <- factor(X_edu, levels = 1:3, labels = c("Low", "Medium", "High"))
  sex             <- factor(rbinom(n, 1, 0.5), levels = 0:1, labels = c("Female", "Male"))
  age             <- round(truncnorm::rtruncnorm(n, a = 45, b = 85, mean = 65, sd = 10))

  edu_med  <- as.integer(education_level == "Medium")
  edu_high <- as.integer(education_level == "High")
  sex_male <- as.integer(sex == "Male")
  age_c    <- as.numeric(scale(age))

  ## exposure A | L: logit P(A = 1) = gamma0 + gamma_scale * (base coefficients) * L
  gamma_base <- c(edu_med = 0.2, edu_high = 0.5, sex_male = 0.2, age_c = -0.4)
  gamma_used <- gamma_scale * gamma_base
  eta_treat  <- gamma_used["edu_med"]  * edu_med +
                gamma_used["edu_high"] * edu_high +
                gamma_used["sex_male"] * sex_male +
                gamma_used["age_c"]    * age_c
  gamma0   <- tune_gamma0_for_treatshare(eta_treat, target = target_treat)
  p_active <- plogis(eta_treat + gamma0)
  A        <- rbinom(n, 1, prob = p_active)
  physical_activity <- factor(A, levels = 0:1, labels = c("Inactive", "Active"))

  ## latent class C | A, L: multinomial logit, class 1 (High-High) is the reference.
  ## Each vector is c(intercept [not used, re-tuned below], A, edu_med, edu_high, sex_male, age_c)
  K <- 4
  alpha_class2 <- c( 1.4, -1.2, -0.4, -0.6,  0.2,  0.01)   # Low-Low
  alpha_class3 <- c(-0.6, -0.1,  0.5,  0.3, -0.2, -0.03)   # High-Low
  alpha_class4 <- c(-0.3, -1.5, -0.3,  0.1,  0.1, -0.02)   # Low-High

  X_class <- cbind(A, edu_med, edu_high, sex_male, age_c)
  B_base <- rbind(
    A        = c(alpha_class2[2], alpha_class3[2], alpha_class4[2]),
    edu_med  = c(alpha_class2[3], alpha_class3[3], alpha_class4[3]),
    edu_high = c(alpha_class2[4], alpha_class3[4], alpha_class4[4]),
    sex_male = c(alpha_class2[5], alpha_class3[5], alpha_class4[5]),
    age_c    = c(alpha_class2[6], alpha_class3[6], alpha_class4[6])
  )

  ## the exposure coefficients are kept fixed; only the confounder coefficients
  ## are scaled, so rho controls confounding and not the overall exposure effect
  B_used <- B_base
  B_used[c("edu_med", "edu_high", "sex_male", "age_c"), ] <-
    alpha_scale * B_base[c("edu_med", "edu_high", "sex_male", "age_c"), ]

  ## class intercepts calibrated to the target class proportions
  alpha0 <- tune_alpha0_for_classmix(X_class, B_used, target_pi)

  eta2 <- alpha0[2] + X_class %*% B_used[, 1]
  eta3 <- alpha0[3] + X_class %*% B_used[, 2]
  eta4 <- alpha0[4] + X_class %*% B_used[, 3]
  eta  <- cbind(rep(0, n), eta2, eta3, eta4)
  p_mat <- exp(eta - apply(eta, 1, max))
  p_mat <- p_mat / rowSums(p_mat)

  set.seed(seed + 1L)
  C <- apply(p_mat, 1, function(prob) sample.int(K, size = 1, prob = prob))
  latent_class <- factor(C, levels = 1:4, labels = SIM$class_names)

  ## outcomes Y | C: independent Normal(theta_c, sigma^2)
  ## base class means (rows HH, LL, HL, LH; columns REY_I ... STP)
  p <- 6
  theta_base <- matrix(c(
    1.5, 1.6, 1.7, 1.8, 0.9, 1.7,   # HH
    0.1, 0.2, 0.1, 0.2, 0.1, 0.2,   # LL
    1.8, 0.9, 0.2, 0.3, 0.2, 0.3,   # HL
    0.2, 0.3, 1.6, 1.7, 1.8, 1.9    # LH
  ), nrow = K, byrow = TRUE)
  colnames(theta_base) <- SIM$outcome_names
  rownames(theta_base) <- c("HH", "LL", "HL", "LH")
  theta_matrix <- scale_separation(theta_base, delta = delta)

  Y <- matrix(NA_real_, nrow = n, ncol = p)
  for (i in 1:n) Y[i, ] <- rnorm(p, mean = theta_matrix[C[i], ], sd = sigma)
  colnames(Y) <- SIM$outcome_names

  simulated_data <- data.frame(
    ID = 1:n,
    latent_class = latent_class,
    physical_activity = physical_activity,
    education_level = education_level,
    sex = sex,
    age = age,
    as.data.frame(Y)
  )

  list(
    simulated_data = simulated_data,
    theta_matrix   = theta_matrix,
    sigma          = sigma,
    gamma0_tuned   = gamma0,
    alpha0_tuned   = alpha0,
    gamma_used     = gamma_used,
    B_used         = B_used,
    gamma_scale    = gamma_scale,
    alpha_scale    = alpha_scale,
    achieved_treat = mean(A),
    achieved_pi    = as.numeric(prop.table(table(factor(C, levels = 1:K)))),
    target_treat   = target_treat,
    target_pi      = target_pi,
    delta          = delta
  )
}

# The data-generating call for one (scenario, seed) pair. Every script goes
# through this function, so the datasets are identical across scripts.
simulate_scenario_dataset <- function(scn, seed) {
  scales <- conf_scales_from_label(scn$confounding)
  simulate_one(
    seed         = seed,
    n            = scn$n,
    target_treat = SIM$target_treat,
    target_pi    = target_pi_from_label(scn$imbalance),
    delta        = delta_from_label(scn$separation),
    sigma        = SIM$sigma,
    gamma_scale  = scales$gamma_scale,
    alpha_scale  = scales$alpha_scale
  )
}


# 3. Proposed model: JAGS -----------------------------------------------------------

# Joint model for the covariates, exposure, latent class and outcomes (Section 2.4
# of the manuscript). Class 1 is the reference class of the multinomial logit.
JAGS_MODEL <- "
model {

  for (i in 1:N) {

    # (1) covariate model  P(L | beta)
    age_z[i] ~ dnorm(mu_age, tau_age)
    sex[i]   ~ dbern(pi_sex)
    edu[i]   ~ dcat(p_edu[1:Jedu])

    edu_med_i[i]  <- equals(edu[i], 2)
    edu_high_i[i] <- equals(edu[i], 3)

    # (2) exposure model  P(A | L, gamma)   (education: 'Low' is the reference)
    logit(pA[i]) <- gamma0
                + gamma_age*age_z[i]
                + gamma_sex*sex[i]
                + gamma_med*edu_med_i[i]
                + gamma_high*edu_high_i[i]
    A[i] ~ dbern(pA[i])

    # (3) class model  P(C | A, L, alpha)   (softmax; class 1 is the reference)
    for (k in 2:K) {
      eta[i,k] <- alpha0[k]
                + alphaA[k]*A[i]
                + alpha_age[k]*age_z[i] + alpha_sex[k]*sex[i]
                + alpha_med[k]*edu_med_i[i] + alpha_high[k]*edu_high_i[i]
      expeta[i,k] <- exp(eta[i,k])
    }
    expeta[i,1] <- 1

    denom[i] <- sum(expeta[i,1:K])
    for (k in 1:K) {
      p_class[i,k] <- expeta[i,k] / denom[i]
    }
    C[i] ~ dcat(p_class[i,1:K])

    # (4) outcome model  P(Y | C, theta)
    for (j in 1:p) {
      Y[i,j] ~ dnorm(theta[C[i], j], tau_y[j])
    }
  }

  # priors
  # beta: covariate model
  mu_age ~ dnorm(0, 1)
  tau_age ~ dgamma(2, 1)
  sd_age <- 1 / sqrt(tau_age)

  pi_sex ~ dbeta(1, 1)

  p_edu[1:Jedu] ~ ddirch(alpha_dirich[1:Jedu])

  # gamma: exposure model
  gamma0     ~ dnorm(0, 1)
  gamma_age  ~ dnorm(0, 1)
  gamma_sex  ~ dnorm(0, 1)
  gamma_med  ~ dnorm(0, 1)
  gamma_high ~ dnorm(0, 1)

  # alpha: class model (class 1 is the reference: all zeros)
  for (k in 2:K) {
    alpha0[k]     ~ dnorm(0, 1)
    alphaA[k]     ~ dnorm(0, 1)
    alpha_age[k]  ~ dnorm(0, 1)
    alpha_sex[k]  ~ dnorm(0, 1)
    alpha_med[k]  ~ dnorm(0, 1)
    alpha_high[k] ~ dnorm(0, 1)
  }

  # theta: class-specific outcome means; tau_y: outcome precisions
  for (k in 1:K) {
    for (j in 1:p) {
      theta[k,j] ~ dnorm(0, 1.0E-4)
    }
  }
  for (j in 1:p) {
    tau_y[j] ~ dgamma(0.1, 0.1)
    sigma_y[j] <- 1 / sqrt(tau_y[j])
  }
}
"

# Fit the joint model by MCMC (3 chains: 4,000 adaptation, 6,000 burn-in and 2,000
# retained iterations each, i.e. 6,000 posterior draws in total).
# Returns the draws (S x parameters) and the class means as an S x K x p array.
fit_one <- function(simdat, n_chains = 3, n_adapt = 4000, n_burn = 6000, n_iter = 2000) {
  suppressPackageStartupMessages(library(rjags))
  simulated_data <- simdat$simulated_data

  A     <- as.integer(simulated_data$physical_activity == "Active")
  sex   <- as.integer(simulated_data$sex == "Male")
  age_z <- as.numeric(scale(as.numeric(simulated_data$age)))

  edu_fac <- factor(simulated_data$education_level, levels = c("Low", "Medium", "High"))
  Jedu    <- nlevels(edu_fac)
  edu_idx <- as.integer(edu_fac)

  Y <- as.matrix(simulated_data[, SIM$outcome_names, drop = FALSE])
  p <- ncol(Y)

  data_jags <- list(
    N = nrow(simulated_data), K = 4, p = p, Jedu = Jedu,
    A = A, age_z = age_z, sex = sex, edu = edu_idx, Y = Y,
    alpha_dirich = rep(1, Jedu)
  )

  jags_model <- jags.model(textConnection(JAGS_MODEL), data = data_jags,
                           n.chains = n_chains, n.adapt = n_adapt)
  update(jags_model, n_burn)

  vars <- c(
    "mu_age", "tau_age", "pi_sex", "p_edu",                               # covariate model
    "gamma0", "gamma_age", "gamma_sex", "gamma_med", "gamma_high",        # exposure model
    "alpha0", "alphaA", "alpha_age", "alpha_sex", "alpha_med", "alpha_high",  # class model
    "theta", "tau_y"                                                      # outcome model
  )
  samples     <- coda.samples(jags_model, variable.names = vars, n.iter = n_iter)
  samples_mat <- as.matrix(samples)

  K <- 4; S <- nrow(samples_mat)
  theta_array <- array(NA, dim = c(S, K, p))
  for (k in 1:K) for (j in 1:p) {
    theta_array[, k, j] <- samples_mat[, paste0("theta[", k, ",", j, "]")]
  }

  list(samples_mat = samples_mat, theta_array = theta_array)
}


# 4. Label switching ------------------------------------------------------------------

# Relabel every posterior draw to the data-generating class labels. In the
# simulation the true class means (Theta_ref) are known, so each draw is matched
# to them with the Hungarian algorithm (clue::solve_LSAP), minimising the squared
# distance between the draw's class means and Theta_ref. The same permutation is
# applied to the class-membership coefficients.
# (This reference-based relabelling is only possible in a simulation; the CLSA
# analysis relabels with the ECR algorithm of the label.switching package.)
relabel_and_align_one <- function(fit, Theta_ref) {
  Theta <- fit$theta_array                              # S x K x p
  S <- dim(Theta)[1]; K <- dim(Theta)[2]

  solve_perm <- function(M_kp, Theta_ref) {
    cost <- outer(1:K, 1:K, Vectorize(function(k, h) sum((M_kp[k, ] - Theta_ref[h, ])^2)))
    as.integer(clue::solve_LSAP(cost))
  }

  perm_mat <- matrix(NA_integer_, nrow = S, ncol = K)
  for (s in 1:S) {
    perm <- solve_perm(Theta[s, , , drop = TRUE], Theta_ref)
    inv  <- integer(K); inv[perm] <- 1:K                # inv[h] = k
    perm_mat[s, ] <- inv
  }

  theta_relab <- array(NA_real_, dim = dim(Theta))
  for (s in 1:S) theta_relab[s, , ] <- Theta[s, perm_mat[s, ], , drop = FALSE]

  # class-model coefficients as S x K matrices (class 1 = reference = 0)
  extract_alpha <- function(base) {
    out <- matrix(0, nrow = S, ncol = K)
    out[, 2:K] <- as.matrix(fit$samples_mat[, paste0(base, "[", 2:K, "]"), drop = FALSE])
    colnames(out) <- paste0("k", 1:K)
    out
  }
  perm_rows <- function(M, perms) {
    out <- M
    for (s in 1:nrow(M)) out[s, ] <- M[s, perms[s, ]]
    out
  }
  alpha_names <- c("alpha0", "alphaA", "alpha_age", "alpha_sex", "alpha_med", "alpha_high")
  alpha_relab <- lapply(setNames(alpha_names, alpha_names),
                        function(nm) perm_rows(extract_alpha(nm), perm_mat))

  list(theta = theta_relab, alpha = alpha_relab, perm_mat = perm_mat)
}


# 5. Posterior predictive standardization -----------------------------------------------

# For every posterior draw s: simulate M_synth covariate vectors from the covariate
# model, evaluate the class probabilities under A = 1 and A = 0, and average.
# Returns psi1 = psi_c(1), psi0 = psi_c(0) and ATE (each S x K).
gcomp_one_dataset <- function(relabeled_global, samples_mat, M_synth = 1000, Jedu = 3) {
  S <- nrow(samples_mat); K <- ncol(relabeled_global$alpha$alpha0)

  mu_age_draw  <- samples_mat[, "mu_age"]
  tau_age_draw <- samples_mat[, "tau_age"]
  pi_sex_draw  <- samples_mat[, "pi_sex"]
  p_edu_draw   <- as.matrix(samples_mat[, paste0("p_edu[", 1:Jedu, "]"), drop = FALSE])

  a0     <- relabeled_global$alpha$alpha0
  aA     <- relabeled_global$alpha$alphaA
  a_age  <- relabeled_global$alpha$alpha_age
  a_sex  <- relabeled_global$alpha$alpha_sex
  a_med  <- relabeled_global$alpha$alpha_med
  a_high <- relabeled_global$alpha$alpha_high

  gcomp_one_draw <- function(s, A_val, M) {
    age_z_s <- rnorm(M, mu_age_draw[s], 1 / sqrt(tau_age_draw[s]))
    sex_s   <- rbinom(M, 1, pi_sex_draw[s])
    edu_s   <- rcat(M, p_edu_draw[s, ])
    med_i   <- as.integer(edu_s == 2L)
    high_i  <- as.integer(edu_s == 3L)

    eta <- outer(rep(1, M), a0[s, ]) +
      A_val * outer(rep(1, M), aA[s, ]) +
      outer(age_z_s, a_age[s, ]) +
      outer(sex_s,   a_sex[s, ]) +
      outer(med_i,   a_med[s, ]) +
      outer(high_i,  a_high[s, ])

    eta <- eta - apply(eta, 1, max)
    ex  <- exp(eta)
    colMeans(ex / rowSums(ex))
  }

  psi1 <- t(vapply(1:S, function(s) gcomp_one_draw(s, 1, M_synth), numeric(K)))
  psi0 <- t(vapply(1:S, function(s) gcomp_one_draw(s, 0, M_synth), numeric(K)))
  list(psi1 = psi1, psi0 = psi0, ATE = psi1 - psi0)
}


# 6. True causal effect -----------------------------------------------------------------

# True ATE_c: evaluate the data-generating class model with A set to 1 and to 0 for
# every subject in the generated dataset and average the difference.
true_ate_by_class <- function(sim) {
  df <- sim$simulated_data
  n  <- nrow(df)
  K  <- 4

  X_noA <- cbind(edu_med  = as.integer(df$education_level == "Medium"),
                 edu_high = as.integer(df$education_level == "High"),
                 sex_male = as.integer(df$sex == "Male"),
                 age_c    = as.numeric(scale(df$age)))
  alpha0 <- sim$alpha0_tuned
  B_used <- sim$B_used

  get_p <- function(A_val) {
    X_class <- cbind(rep(A_val, n), X_noA)
    Lin <- matrix(0, nrow = n, ncol = K)
    Lin[, 2:K] <- sweep(X_class %*% B_used, 2, alpha0[2:K], FUN = "+")
    Lin <- Lin - apply(Lin, 1, max)
    ex  <- exp(Lin)
    ex / rowSums(ex)
  }

  colMeans(get_p(1L) - get_p(0L))
}


# 7. Comparator: two-step IPW latent class estimator ---------------------------------------
# Step 1: stabilized ATE weights from a logistic exposure model.
# Step 2: weighted Gaussian finite mixture (flexmix) with a multinomial
#         concomitant model for the class given the exposure; ATE_c is the
#         difference of the fitted class probabilities at A = 1 and A = 0.

# Stabilized inverse probability of treatment weights
weighting <- function(dat, A_var, X_vars) {
  f  <- as.formula(paste(A_var, "~", paste(X_vars, collapse = " + ")))
  df <- dat[, c(A_var, X_vars), drop = FALSE]
  df[[A_var]] <- as.integer(df[[A_var]] == "Active")
  wobj <- WeightIt::weightit(f, data = df, estimand = "ATE", method = "glm", stabilize = TRUE)
  as.numeric(wobj$weights)
}

# flexmix needs integer frequency weights: rescale to mean 1, multiply by `scale`,
# round, and bound below at 1 and above at `cap` (relative weights 0.05 to 5)
as_rep_weights <- function(w, cap = 100, scale = 20) {
  w_norm <- w / mean(w)
  w_rep  <- pmax(1, as.integer(round(w_norm * scale)))
  w_rep[w_rep > cap] <- cap
  as.integer(w_rep)
}

# Weighted K-component Gaussian mixture, class ~ exposure (k-means initialisation)
fit_flexmix_once <- function(dat, outcome_names, K, w_rep) {
  Ydf <- as.data.frame(dat[, outcome_names, drop = FALSE])
  colnames(Ydf) <- outcome_names
  df_fm <- Ydf
  df_fm$A <- as.integer(dat$physical_activity == "Active")

  models <- lapply(colnames(Ydf), function(v)
    flexmix::FLXMRglm(as.formula(paste(v, "~ 1")), family = "gaussian"))
  concom <- flexmix::FLXPmultinom(~ A)

  km <- stats::kmeans(scale(Ydf), centers = K, nstart = 50)$cluster
  n  <- nrow(Ydf)
  Z  <- matrix(1e-6, n, K); Z[cbind(seq_len(n), km)] <- 1; Z <- Z / rowSums(Z)

  flexmix::flexmix(
    formula     = as.formula(paste(colnames(Ydf)[1], "~ 1")),
    data        = df_fm,
    k           = K,
    model       = models,
    concomitant = concom,
    weights     = w_rep,
    cluster     = Z,
    control     = list(minprior = 0)   # flexmix drops classes below 5% by default; keep all
  )
}

# Permutation matching the fitted classes to the true class means (Hungarian
# algorithm on the weighted class means), as in relabel_and_align_one()
perm_from_reference <- function(dat, outcome_names, K, w_ate, fit, Theta_ref) {
  Ydf       <- as.data.frame(dat[, outcome_names, drop = FALSE])
  class_hat <- as.integer(flexmix::clusters(fit))

  theta_hat <- matrix(NA_real_, nrow = K, ncol = ncol(Ydf))
  for (k in 1:K) {
    idx <- which(class_hat == k)
    if (length(idx) == 0L) next
    Yk <- as.matrix(Ydf[idx, , drop = FALSE])
    theta_hat[k, ] <- colSums(w_ate[idx] * Yk) / sum(w_ate[idx])
  }

  cost <- matrix(0, K, K)
  for (k in 1:K) for (h in 1:K) cost[k, h] <- sum((theta_hat[k, ] - Theta_ref[h, ])^2, na.rm = TRUE)

  assign <- clue::solve_LSAP(cost)
  perm <- integer(K)
  perm[assign] <- seq_len(K)
  perm
}

# ATE_c from the fitted concomitant (class | exposure) model, in true-class order
ate_from_concomitant <- function(fit, perm, class_names = SIM$class_names) {
  B_reduced <- fit@concomitant@coef
  P1 <- softmax_vec(cbind(1, 1) %*% B_reduced)
  P0 <- softmax_vec(cbind(1, 0) %*% B_reduced)
  ATE <- as.numeric(P1 - P0)[perm]
  names(ATE) <- class_names
  ATE
}

# Point estimate of ATE_c (both steps)
compute_ate_point <- function(dat, outcome_names, K, Theta_ref) {
  w_ate <- weighting(dat, A_var = "physical_activity", X_vars = c("education_level", "sex", "age"))
  fit   <- fit_flexmix_once(dat, outcome_names, K, as_rep_weights(w_ate))
  perm  <- perm_from_reference(dat, outcome_names, K, w_ate, fit, Theta_ref)
  ate_from_concomitant(fit, perm)
}

# Nonparametric bootstrap (both steps repeated in every resample);
# 95% interval = 2.5th and 97.5th percentiles of the bootstrap distribution.
# Resamples in which flexmix fails are dropped (their number is nrow(draws)).
bootstrap_ate <- function(dat, outcome_names = SIM$outcome_names, K = 4,
                          Theta_ref, B = 50, seed = 20251007) {
  set.seed(seed)
  ate_hat  <- compute_ate_point(dat, outcome_names, K, Theta_ref)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = K)
  colnames(boot_mat) <- names(ate_hat)
  n <- nrow(dat)
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    boot_mat[b, ] <- tryCatch(compute_ate_point(dat[idx, , drop = FALSE], outcome_names, K, Theta_ref),
                              error = function(e) rep(NA_real_, K))
  }
  boot_mat <- boot_mat[stats::complete.cases(boot_mat), , drop = FALSE]

  list(ate_hat = ate_hat,
       ci95    = cbind(`2.5%`  = apply(boot_mat, 2, stats::quantile, probs = 0.025),
                       `97.5%` = apply(boot_mat, 2, stats::quantile, probs = 0.975)),
       draws   = boot_mat)
}


# 8. One simulated dataset, end to end -----------------------------------------------------------

# Simulate one dataset, fit the proposed model and the comparator, and return the
# estimates, the truth and the 95% interval coverage of both methods.
#   M_synth   covariate vectors per posterior draw in the standardization
#   B_boot    bootstrap replicates of the comparator
run_one_in_scenario <- function(scn, seed, M_synth = 1000, B_boot = 50) {
  scales <- conf_scales_from_label(scn$confounding)
  sim    <- simulate_scenario_dataset(scn, seed)

  # proposed model: posterior of ATE_c and 95% credible interval
  fit      <- fit_one(sim)
  aligned  <- relabel_and_align_one(fit, Theta_ref = sim$theta_matrix)
  g        <- gcomp_one_dataset(aligned, fit$samples_mat, M_synth = M_synth)
  ATE_true <- true_ate_by_class(sim)

  ATE_mean  <- colMeans(g$ATE)
  ATE_lwr   <- apply(g$ATE, 2, quantile, probs = 0.025)
  ATE_upr   <- apply(g$ATE, 2, quantile, probs = 0.975)
  cover     <- as.integer(ATE_true >= ATE_lwr & ATE_true <= ATE_upr)
  names(cover) <- paste0("class", 1:4)

  # comparator: point estimate and bootstrap interval
  base_res   <- bootstrap_ate(sim$simulated_data, outcome_names = SIM$outcome_names, K = 4,
                              Theta_ref = sim$theta_matrix, B = B_boot, seed = seed + 10000)
  base_cover <- as.integer(ATE_true >= base_res$ci95[, "2.5%"] & ATE_true <= base_res$ci95[, "97.5%"])
  names(base_cover) <- paste0("class", 1:4)

  list(
    scenario_id = scn$scenario_id, tag = scenario_tag(scn),
    confounding = scn$confounding, imbalance = scn$imbalance,
    separation = scn$separation, n = scn$n,
    gamma_scale = scales$gamma_scale, alpha_scale = scales$alpha_scale,
    delta = sim$delta, sigma = sim$sigma,
    target_pi = sim$target_pi, target_treat = sim$target_treat,
    achieved_treat = sim$achieved_treat, achieved_pi = sim$achieved_pi,
    ATE_true = ATE_true,
    # proposed model
    ATE_mean = ATE_mean, ATE_lwr = ATE_lwr, ATE_upr = ATE_upr, cover = cover,
    # comparator
    base_mean = base_res$ate_hat,
    base_lwr = base_res$ci95[, "2.5%"], base_upr = base_res$ci95[, "97.5%"],
    base_cover = base_cover
  )
}

# Comparator only, one dataset (used by lanza_bootstrap_B1000.R). Uses the same
# data-generating call and the same bootstrap seed (seed + 10000) as
# run_one_in_scenario(), so the first B = 50 resamples reproduce that run.
run_comparator_one <- function(scn, seed, B_boot) {
  sim      <- simulate_scenario_dataset(scn, seed)
  ATE_true <- true_ate_by_class(sim)
  base_res <- bootstrap_ate(sim$simulated_data, outcome_names = SIM$outcome_names, K = 4,
                            Theta_ref = sim$theta_matrix, B = B_boot, seed = seed + 10000)
  base_lwr <- base_res$ci95[, "2.5%"]
  base_upr <- base_res$ci95[, "97.5%"]

  list(
    scenario_id = scn$scenario_id, tag = scenario_tag(scn), n = scn$n,
    confounding = scn$confounding, imbalance = scn$imbalance, separation = scn$separation,
    seed = seed, B_requested = B_boot,
    B_effective = nrow(base_res$draws),          # resamples in which flexmix converged
    base_mean = base_res$ate_hat, base_lwr = base_lwr, base_upr = base_upr,
    base_true = ATE_true,
    base_cover = as.integer(ATE_true >= base_lwr & ATE_true <= base_upr)
  )
}


# 9. Summaries -----------------------------------------------------------------------------------

# Monte Carlo summary of one scenario for both methods: mean squared error and
# coverage of the 95% interval for ATE_c across the datasets. `resL` is the list of
# run_one_in_scenario() results of the scenario.
summarize_scenario <- function(scn, resL) {
  summarize_method <- function(method_name) {
    if (method_name == "Bayesian") {
      est_mat <- do.call(rbind, lapply(resL, function(rr) rr$ATE_mean))
    } else {
      est_mat <- do.call(rbind, lapply(resL, function(rr) rr$base_mean))
    }
    cover_mat <- do.call(rbind, lapply(resL, function(rr)
      if (method_name == "Bayesian") rr$cover else rr$base_cover))
    true_mat  <- do.call(rbind, lapply(resL, function(rr) rr$ATE_true))
    data.frame(
      Scenario_ID = scn$scenario_id, n = scn$n,
      Conf = scn$confounding, Imb = scn$imbalance, Sep = scn$separation,
      Method   = method_name,
      Class    = paste0("class", 1:4),
      MSE      = as.numeric(colMeans((est_mat - true_mat)^2, na.rm = TRUE)),
      Coverage = as.numeric(colMeans(cover_mat, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  rbind(summarize_method("Bayesian"), summarize_method("Frequentist"))
}

# Scenario settings and the achieved exposure prevalence / class proportions
# (averaged over the datasets of the scenario)
scenario_index_row <- function(scn, resL) {
  achieved_pi <- colMeans(do.call(rbind, lapply(resL, function(rr) rr$achieved_pi)), na.rm = TRUE)
  data.frame(
    Scenario       = scenario_tag(scn),
    n              = scn$n,
    confounding    = scn$confounding,
    imbalance      = scn$imbalance,
    separation     = scn$separation,
    gamma_scale    = resL[[1]]$gamma_scale,
    alpha_scale    = resL[[1]]$alpha_scale,
    delta          = resL[[1]]$delta,
    sigma          = resL[[1]]$sigma,
    target_treat   = resL[[1]]$target_treat,
    achieved_treat = mean(vapply(resL, function(rr) rr$achieved_treat, numeric(1)), na.rm = TRUE),
    target_pi      = paste(resL[[1]]$target_pi, collapse = ";"),
    achieved_pi    = paste(round(achieved_pi, 4), collapse = ";"),
    stringsAsFactors = FALSE
  )
}


# 10. Utilities -----------------------------------------------------------------------------

# Parse a selector such as "1,3,5-8" into integer ids ("" selects all of all_ids)
parse_ids <- function(txt, all_ids) {
  if (!nzchar(txt)) return(all_ids)
  out <- unlist(lapply(strsplit(txt, ",")[[1]], function(tok) {
    tok <- trimws(tok)
    if (grepl("-", tok, fixed = TRUE)) {
      ab <- as.integer(strsplit(tok, "-")[[1]]); seq(ab[1], ab[2])
    } else {
      as.integer(tok)
    }
  }))
  out <- sort(unique(out))
  if (!all(out %in% all_ids)) stop("id(s) out of range: ", paste(setdiff(out, all_ids), collapse = ","))
  out
}

# Number of parallel workers: WORKERS, else the Slurm allocation, else 1
get_workers <- function() {
  as.integer(Sys.getenv("WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
}
