#!/usr/bin/env Rscript
library(tidyverse)

dirg <- "D:/JS"
setwd(dirg)

# ----------------------------
# Read follow-up times
# ----------------------------
long.time <- read.csv("long.data_new600.csv")
first.tt <- long.time[, 2]
last.tt  <- long.time[, 3]

N  <- length(last.tt)
id <- 1:N

t  <- round(first.tt)
tt <- round(last.tt)

k.pa <- (tt - t) * 4
kk   <- max(k.pa)

# ----------------------------
# Number of simulated datasets
# ----------------------------
Int <- 501
set.seed(123)

# ----------------------------
# TRUE VALUES (match JAGS names)
# ----------------------------
# PA
c10 <- -4
c20 <- -2
c   <- c(0.4, 0.2, 0.3, -0.05)  # c[1], c[2], c[3], c[4]
cp1_mu_true <- 9
cp1_sd_true <- 1
pi_true <- c(0.5, 0.5)          # Pr(z=1), Pr(z=2) for PA component

# PE
a1_true  <- 1.1
a2_true  <- 0.9
b10_true <- -1.5
b20_true <- -2
b_true   <- c(0.2, 0.3)           # b[1], b[2]
ga10_true <- 0.3 
ga20_true <- 0.1
ga11_true <- -0.2
ga12_true <- 0.1
pi_r_true <- c(0.5, 0.5)          # Pr(z.r=1), Pr(z.r=2) for PE component

# Random effects truths (precisions)
u_tau1_true <- 4             # sd(u1)=0.5
u_tau2_true <- 4            # sd(u2)=0.5
w_tau1_true <- 25            # sd(w1)=0.2
w_tau2_true <- 25               # sd(w2)=0.2

# ----------------------------
# NHPP generator (0..T)
# Returns interval rows + a terminal censoring row at T with status=0
# ----------------------------
NHPP <- function(a, b, T) {
  mu <- b * T^a
  n  <- rpois(1, mu)
  
  if (n != 0) {
    u <- sort(runif(n, 0, 1))
    y <- T * u^(1/a)
    y[length(y) + 1] <- T
    
    y_0 <- rep(NA_real_, length(y))
    for (ii in 2:length(y_0)) y_0[ii] <- y[ii - 1]
    y_0[is.na(y_0)] <- 0
    
    cbind(y_0, y, c(rep(1, length(y) - 1), 0), n)
  } else {
    cbind(0, T, 0, n)
  }
}

# ----------------------------
# Simulation loop
# ----------------------------
for (r in 2:Int) {
  
  # Treatment vector (robust if N odd)
  X1 <- c(rep(1, floor(N/2)), rep(0, N - floor(N/2)))
  
  # Latent variables exactly as in JAGS
  u1  <- rnorm(N, 0, sqrt(1 / u_tau1_true))
  u2  <- rnorm(N, 0, sqrt(1 / u_tau2_true))
  w1  <- rnorm(N, 0, sqrt(1 / w_tau1_true))
  w2  <- rnorm(N, 0, sqrt(1 / w_tau2_true))
  cp1 <- rnorm(N, cp1_mu_true, cp1_sd_true)
  
  while(any(cp1 > 21)){
    ind <- which(cp1 > 21)
    cp1[ind] <- rnorm(length(ind), cp1_mu_true, cp1_sd_true)
  }
  
  # Second change point:
  # cp2 | cp1 ~ Uniform(cp1, 21)
  cp2 <- runif(
    N,
    min = cp1,
    max = 21
  )
  
  # Class indicators (truth)
  z   <- sample(1:2, size = N, prob = pi_true,  replace = TRUE)   # PA class
  z_r <- sample(1:2, size = N, prob = pi_r_true, replace = TRUE)  # PE class
  
  # Frailties
  cp1c <- cp1 - cp1_mu_true
  # Conditional mean of cp2 given cp1 under Uniform(cp1, 21)
  cp2_mu <- (cp1 + 21) / 2
  cp2c <- cp2 - cp2_mu
  
  # Shared longitudinal contribution
  eta_shared <- ifelse(
    z == 1,
    ga10_true * u1 + ga11_true * cp1c + ga12_true * cp2c,
    ga20_true * u2
  )
  
  # PE-component-specific frailties
  v1 <- exp(w1 + eta_shared)
  v2 <- exp(w2 + eta_shared)
  
  # ----------------------------
  # Build PA design time matrix X and responses Y
  # ----------------------------
  X <- matrix(NA_real_, nrow = N, ncol = kk)
  for (i in 1:N) {
    X[i, 1:k.pa[i]] <- seq(t[i], tt[i] - 0.25, by = 0.25)
  }
  
  Y <- matrix(NA_integer_, nrow = N, ncol = kk)
  
  for (i in 1:N) {
    for (j in 1:k.pa[i]) {
      Ind1 <- ifelse(X[i, j] < cp1[i], -1, 1)
      Ind2 <- ifelse(X[i, j] < cp2[i], -1, 1)
      
      logit1 <- c10 +
        c[1] * (X[i, j] - cp1[i]) +
        c[2] * (X[i, j] - cp1[i]) * Ind1 +
        c[3] * (X[i, j] - cp2[i]) * Ind2 +
        c[4] * X1[i] + u1[i]
      
      logit2 <- c20 +
        (c[1] - c[2] - c[3]) * X[i, j] +
        c[4] * X1[i] + u2[i]
      
      p1 <- plogis(logit1)
      p2 <- plogis(logit2)
      p  <- if (z[i] == 1) p1 else p2
      
      Y[i, j] <- rbinom(1, 1, p)
    }
  }
  
  # ----------------------------
  # Simulate PE recurrent events
  # ----------------------------
  tau <- tt - 0.25
  sim_pe_list <- vector("list", N)
  
  for (i in 1:N) {
    
    if (z_r[i] == 1) {
      b_scale <- v1[i] * exp(b10_true + b_true[1] * X1[i])
      times <- NHPP(a = a1_true, b = b_scale, T = tau[i])
    } else {
      b_scale <- v2[i] * exp(b20_true + b_true[2] * X1[i])
      times <- NHPP(a = a2_true, b = b_scale, T = tau[i])
    }
    
    dat_i <- data.frame(
      id     = i,
      xi     = X1[i],
      Tei    = tau[i],
      n.rec  = as.integer(times[, 4]),
      start  = as.numeric(times[, 1]),
      stop   = as.numeric(times[, 2]),
      status = as.integer(times[, 3])
    )
    
    sim_pe_list[[i]] <- dat_i
  }
  
  simdat.pe <- bind_rows(sim_pe_list) %>%
    arrange(id, stop)
  
  # ----------------------------
  # Write outputs
  # ----------------------------
  write.csv(as.data.frame(X),
            file = paste0("mixJM.X_newdata7.", r - 2, ".csv"),
            row.names = FALSE)
  
  write.csv(as.data.frame(Y),
            file = paste0("mixJM.Y_newdata7.", r - 2, ".csv"),
            row.names = FALSE)
  
  write.csv(simdat.pe,
            file = paste0("mixJM.rec_newdata7.", r - 2, ".csv"),
            row.names = FALSE)
  
  # ----------------------------
  # Write truth file for recovery checks
  # ----------------------------
  truth <- data.frame(
    c10 = c10, c20 = c20,
    c1 = c[1], c2 = c[2], c3 = c[3], c4 = c[4],
    B1 = c[1] - c[2] - c[3],
    B2 = c[1] + c[2] - c[3],
    B3 = c[1] + c[2] + c[3],
    cp1_mu = cp1_mu_true, cp1_tau = 1 / (cp1_sd_true^2), cp2_mean = (cp1_mu_true + 21) / 2,
    a1 = a1_true, a2 = a2_true,
    b10 = b10_true, b20 = b20_true, b1 = b_true[1], b2 = b_true[2],
    ga10 = ga10_true, ga20 = ga20_true,
    ga11 = ga11_true, ga12 = ga12_true,
    pi1 = pi_true[1], pi2 = pi_true[2],
    pi1r = pi_r_true[1], pi2r = pi_r_true[2],
    u_tau1 = u_tau1_true, u_tau2 = u_tau2_true,
    w_tau1 = w_tau1_true, w_tau2 = w_tau2_true
  )
  
  write.csv(truth,
            file = paste0("mixJM.truth7.", r - 2, ".csv"),
            row.names = FALSE)
}