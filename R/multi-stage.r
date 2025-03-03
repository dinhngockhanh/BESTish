if (!require("deSolve")) install.packages("deSolve")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("reshape2")) install.packages("reshape2")
if (!require("patchwork")) install.packages("patchwork")
if (!require("zoo")) install.packages("zoo") # for na.approx

library(deSolve)
library(ggplot2)
library(dplyr)
library(reshape2)
library(patchwork)
library(zoo)
library(numDeriv)

# 1) SINGLE-RUN TAU-LEAP SIMULATION
simulate_continuous_moran_tau <- function(
    r, # initial # of type-0 cells
    lambda_vec, # e.g. c(lambda0, lambda1, lambda2)
    u_vec, # e.g. c(u0, u1, u2)
    alpha, # e.g. 1.0
    max_time, # e.g. 100
    tau, # step for tau-leaping
    seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    n_types <- length(lambda_vec)
    state <- numeric(n_types)
    state[1] <- r # all cells in type-0

    t <- 0
    record_time <- c(t)
    record_state <- matrix(state, nrow = 1)

    birth_no_mut <- (1 - u_vec) * lambda_vec
    birth_mut <- u_vec * lambda_vec

    while (t < max_time) {
        N_total <- sum(state)
        if (N_total <= 0) break

        # (A) Divisions w/o mutation
        divisions <- rpois(n_types, birth_no_mut * state * tau)
        # (B) Mutations
        mutations <- rpois(n_types, birth_mut * state * tau)
        # (C) Deaths
        death_rate <- alpha * sum(birth_no_mut * state)
        deaths <- rpois(n_types, death_rate * (state / N_total) * tau)

        t <- t + tau
        new_state <- state + divisions - mutations - deaths

        # Shift mutated cells from type j -> j+1
        for (j in seq_len(n_types)) {
            if (j < n_types) {
                new_state[j + 1] <- new_state[j + 1] + mutations[j]
            } else {
                new_state[j] <- new_state[j] + mutations[j]
            }
        }

        new_state[new_state < 0] <- 0
        state <- new_state

        record_time <- c(record_time, t)
        record_state <- rbind(record_state, state)
    }

    out_df <- data.frame(time = record_time)
    for (j in seq_len(n_types)) {
        out_df[[paste0("N", j - 1)]] <- record_state[, j]
        out_df[[paste0("Nbar", j - 1)]] <- record_state[, j] / r
    }
    out_df
}

# My code does n_types < length(lambda(lambda_vec) but then in usage it is always 4.
# Your code also does n_types <- length(lambda_vec)
# but I think it was originally used for 3-compartment scenario (N0, N1, N2 in the code snippet)
# otherwise the structure inside is basically same.

# 2) MULTIPLE REPLICATES
run_replicates <- function(
    n_reps,
    r,
    lambda_vec,
    u_vec,
    alpha,
    max_time,
    tau,
    seed = 123) {
    set.seed(seed)
    sim_list <- vector("list", n_reps)
    for (i in seq_len(n_reps)) {
        sim_i <- simulate_continuous_moran_tau(r, lambda_vec, u_vec, alpha, max_time, tau)
        sim_i$replicate <- i
        sim_list[[i]] <- sim_i
    }
    dplyr::bind_rows(sim_list)
}

# I think this part is also same structure.

# 3) SOLVE THE ODE FOR THE MEAN (ALPHA=1)
two_hit_ode_mean <- function(t, state, params) {
    N0 <- state[1]
    N1 <- state[2]
    N2 <- state[3]

    w0 <- params$w[1]
    v0 <- params$v[1]
    w1 <- params$w[2]
    v1 <- params$v[2]
    w2 <- params$w[3]
    v2 <- params$v[3]

    sumN <- N0 + N1 + N2
    sum_wv <- (w0 + v0) * N0 + (w1 + v1) * N1 + (w2 + v2) * N2

    dN0 <- w0 * N0 - (N0 / sumN) * sum_wv
    dN1 <- w1 * N1 + v0 * N0 - (N1 / sumN) * sum_wv
    dN2 <- w2 * N2 + v1 * N1 - (N2 / sumN) * sum_wv
    list(c(dN0, dN1, dN2))
}

solve_mean_ode <- function(lambda_vec, u_vec, max_time, dt = 0.1) {
    w <- (1 - 2 * u_vec) * lambda_vec
    v <- u_vec * lambda_vec
    params <- list(w = w, v = v)

    state_init <- c(1, 0, 0)
    times <- seq(0, max_time, by = dt)

    ode_out <- deSolve::ode(
        y     = state_init,
        times = times,
        func  = two_hit_ode_mean,
        parms = params
    )

    df <- as.data.frame(ode_out)
    names(df)[2:4] <- c("Nbar0", "Nbar1", "Nbar2")
    df
}

ode_objective_var <- function(t, state, params) {
    f <- function(x) {
        output <- params$w * x + c(0, params$v[1:(params$nCompartments - 1)] * x[1:(params$nCompartments - 1)]) - params$alpha * x / sum(x) * sum((params$w + params$v) * x)
        return(output)
    }

    V <- matrix(state, nrow = params$nCompartments, byrow = TRUE)

    tmp <- which(params$mean_ode_df$time <= t)
    N <- as.numeric(params$mean_ode_df[tmp[length(tmp)], paste0("Nbar", 0:(params$nCompartments - 1))])

    S_diag <- (params$w + 2 * params$v) * N + params$alpha * sum((params$w + params$v) * N) * N / sum(N) + c(0, params$v[1:(params$nCompartments - 1)] * N[1:(params$nCompartments - 1)])
    S_offdiag <- -params$v[1:(params$nCompartments - 1)] * N[1:(params$nCompartments - 1)]
    S <- matrix(0, nrow = params$nCompartments, ncol = params$nCompartments)
    for (j in 1:params$nCompartments) S[j, j] <- S_diag[j]
    for (j in 1:(params$nCompartments - 1)) S[j, j + 1] <- S_offdiag[j]

    A <- jacobian(f, N)
    dVar <- A %*% V + t(A %*% V) + S

    list(dVar)
}

solve_variance_ode <- function(mean_ode_df, lambda_vec, u_vec, alpha) {
    w <- (1 - 2 * u_vec) * lambda_vec
    v <- u_vec * lambda_vec
    nCompartments <- ncol(mean_ode_df) - 1

    params <- list()
    params$w <- w
    params$v <- v
    params$nCompartments <- nCompartments
    params$alpha <- alpha
    params$mean_ode_df <- mean_ode_df


    ode_out <- deSolve::ode(
        y     = rep(0, nCompartments^2),
        times = mean_ode_df$time,
        func  = ode_objective_var,
        parms = params
    )

    variance_matrix <- as.data.frame(ode_out)
    df_names <- c()
    for (i in 0:(nCompartments - 1)) {
        for (j in 0:(nCompartments - 1)) {
            df_names[length(df_names) + 1] <- paste0("var_", i, "_", j)
        }
    }
    names(variance_matrix)[2:ncol(variance_matrix)] <- df_names


    res <- data.frame(time = mean_ode_df$time)
    for (i in 0:(nCompartments - 1)) {
        res[[paste0("varNhat", i, "_the")]] <- variance_matrix[[paste0("var_", i, "_", i)]]
    }
    print(res)
    return(res)
}

solve_variance_ode_OLD <- function(mean_ode_df, lambda_vec, u_vec, alpha) {
    w <- (1 - 2 * u_vec) * lambda_vec
    v <- u_vec * lambda_vec
    nCompartments <- ncol(mean_ode_df) - 1

    f <- function(x) {
        output <- w * x + c(0, v[1:(nCompartments - 1)] * x[1:(nCompartments - 1)]) - alpha * x / sum(x) * sum((w + v) * x)
        return(output)
    }

    variance_matrix <- data.frame(matrix(0, nrow = length(mean_ode_df$time), ncol = 0))
    for (i in 0:(nCompartments - 1)) {
        for (j in 0:(nCompartments - 1)) {
            variance_matrix[[paste0("var_", i, "_", j)]] <- rep(0, length(mean_ode_df$time))
        }
    }
    for (i in 1:(nrow(variance_matrix) - 1)) {
        V <- matrix(as.numeric(variance_matrix[i, ]), nrow = nCompartments, byrow = TRUE)
        N <- as.numeric(mean_ode_df[i, paste0("Nbar", 0:(nCompartments - 1))])

        S_diag <- (w + 2 * v) * N + alpha * sum((w + v) * N) * N / sum(N) + c(0, v[1:(nCompartments - 1)] * N[1:(nCompartments - 1)])
        S_offdiag <- -v[1:(nCompartments - 1)] * N[1:(nCompartments - 1)]
        S <- matrix(0, nrow = nCompartments, ncol = nCompartments)
        for (j in 1:nCompartments) S[j, j] <- S_diag[j]
        for (j in 1:(nCompartments - 1)) S[j, j + 1] <- S_offdiag[j]

        A <- jacobian(f, N)
        update <- A %*% V + t(A %*% V) + S

        next_V <- V + update * (mean_ode_df$time[i + 1] - mean_ode_df$time[i])
        variance_matrix[i + 1, ] <- as.vector(t(next_V))
    }
    res <- data.frame(time = mean_ode_df$time)
    for (i in 0:(nCompartments - 1)) {
        res[[paste0("varNhat", i, "_the")]] <- variance_matrix[[paste0("var_", i, "_", i)]]
    }
    return(res)
}
