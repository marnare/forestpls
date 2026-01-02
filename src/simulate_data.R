simulate_policy_data <- function(
    
  N = 2000,
  p_cont = 5, p_bin = 5, p_cnt = 5,
  noise = c("t3", "cauchy", "mixture", "het_t"),
  effect_noise = c("none", "t3", "cauchy", "mixture", "het_t"),
  treat_mech = c("rct", "logit"),
  clip_ps = c(0.02, 0.98),
  seed = 1212
  
) {
  noise <- match.arg(noise)
  effect_noise <- match.arg(effect_noise)
  treat_mech <- match.arg(treat_mech)
  
  set.seed(seed)
  
  # --- covariates: continuous + binary + count ---
  X_cont <- sapply(1:p_cont, function(j) rnorm(N, mean = (j - 1)/10 - 1, sd = runif(1, 0.5, 2)))
  X_bin  <- sapply(1:p_bin,  function(j) rbinom(N, 1, prob = plogis(rnorm(1, 0, 1))))
  X_cnt  <- sapply(1:p_cnt,  function(j) rpois(N, lambda = runif(1, 0.5, 5)))
  
  X <- cbind(X_cont, X_bin, X_cnt)
  X <- as.matrix(X)
  colnames(X) <- c(paste0("Xc", 1:p_cont), paste0("Xb", 1:p_bin), paste0("Xn", 1:p_cnt))
  
  # helper (safe if column doesn't exist)
  get_col <- function(j) if (ncol(X) >= j) X[, j] else rep(0, N)
  
  x1 <- get_col(1); x2 <- get_col(2); x3 <- get_col(3); x4 <- get_col(4)
  x5 <- get_col(5); x6 <- get_col(6); x7 <- get_col(7)
  
  # --- treatment assignment ---
  if (treat_mech == "rct") {
    ps <- rep(0.5, N)
  } else {
    lp <- -0.3 +
      0.8*as.numeric(scale(x1)) - 0.6*as.numeric(scale(x2)) + 0.5*as.numeric(scale(x3)) +
      0.4*(x1*x4) - 0.3*abs(x5) +
      0.6*(x6 == 1) - 0.4*(x7 == 1)
    ps <- plogis(lp)
    ps <- pmin(pmax(ps, clip_ps[1]), clip_ps[2])
  }
  
  P <- rbinom(N, 1, ps)
  
  # --- baseline outcome and heterogeneous effect ---
  mu0   <- 100*x1 + 100*x2 + 10*x3 - 5*x4
  tau_x <- 2.0*x1 - 1.5*x2 + 1.0*x3*x4 + 0.8*sin(x5) + 0.5*x6
  
  # --- noise generators ---
  gen_noise <- function(type) {
    switch(
      none    = rep(0, N),
      type,
      t3      = rt(N, df = 3),
      mixture = {
        is_out <- rbinom(N, 1, 0.03)
        rnorm(N, 0, 1) + is_out * rnorm(N, 0, 20)
      },
      het_t   = (0.5 + abs(x4)) * rt(N, df = 4),
      stop("Unknown noise type: ", type)
    )
  }
  
  gen_effect_multiplier <- function(type) {
    switch(
      type,
      none    = rep(1, N),
      t3      = 1 + 0.3*rt(N, df = 3),                     # mean ~ 1, heavy tails
      mixture = {
        is_out <- rbinom(N, 1, 0.03)
        1 + rnorm(N, 0, 0.2) + is_out * rnorm(N, 0, 3)
      },
      het_t   = 1 + 0.2*((0.5 + abs(x4)) * rt(N, df = 4)),
      stop("Unknown effect_noise type: ", type)
    )
  }
  
  eps <- gen_noise(noise)
  U   <- gen_effect_multiplier(effect_noise)
  
  # --- observed outcome ---
  Y <- mu0 + P * (tau_x * U) + eps
  
  df <- data.frame(
    Y = as.numeric(Y),
    P = as.integer(P),
    ps = as.numeric(ps),
    tau_true = as.numeric(tau_x),
    mu0 = as.numeric(mu0),
    X
  )
  
  list(df = df, X = X, Y = Y, P = P, ps = ps, mu0 = mu0, tau_true = tau_x)
}

