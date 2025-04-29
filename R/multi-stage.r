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
simulate_continuous_moran_tau <- function(r,
                                          lambda_vec,
                                          u_vec,
                                          alpha,
                                          max_time,
                                          tau = 0.001,
                                          seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    #---Initiate the system
    n_types <- length(lambda_vec)
    if (length(r) == 1) {
        state <- numeric(n_types)
        state[1] <- r
    } else {
        state <- r
    }
    t <- 0
    #---Tau-leaping algorithm
    record_time <- seq(0, max_time, by = tau)
    record_state <- matrix(NA, nrow = length(record_time), ncol = n_types)
    record_state[1, ] <- state
    for (index in 2:length(record_time)) {
        N_total <- sum(state)
        if (N_total <= 0) {
            record_time <- record_time[1:index]
            record_state <- record_state[1:index, ]
            break
        }
        #   Simulate number of divisions into same compartments
        divisions <- rpois(n_types, (1 - u_vec) * lambda_vec * state * tau)
        #   Simulate number of mutations into new compartments
        mutations <- rpois(n_types, u_vec * lambda_vec * state * tau)
        #   Simulate death counts
        deaths <- rpois(n_types, alpha * sum((1 - u_vec) * lambda_vec * state) * state / sum(state) * tau)
        #   Update time
        t <- record_time[index]
        #   Update compartment cell counts
        state <- state - deaths + divisions - mutations +
            c(0, mutations[1:(length(mutations) - 1)])
        state[state < 0] <- 0
        #   Record state and time
        record_state[index, ] <- state
    }
    #   Output system states and times
    out_df <- data.frame(time = record_time)
    for (j in seq_len(n_types)) {
        out_df[[paste0("N", j - 1)]] <- record_state[, j]
        out_df[[paste0("Nbar", j - 1)]] <- record_state[, j] / r
    }
    return(out_df)
}

# My code does n_types < length(lambda(lambda_vec) but then in usage it is always 4.
# Your code also does n_types <- length(lambda_vec)
# but I think it was originally used for 3-compartment scenario (N0, N1, N2 in the code snippet)
# otherwise the structure inside is basically same.

# 2) MULTIPLE REPLICATES
run_replicates <- function(n_reps,
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
        print(i)
        sim_i <- simulate_continuous_moran_tau(r, lambda_vec, u_vec, alpha, max_time, tau)
        sim_i$replicate <- i
        sim_list[[i]] <- sim_i
    }
    dplyr::bind_rows(sim_list)
}

mean_ODE_equations <- function(t,
                               state,
                               params) {
    n_comp <- length(state)
    w <- params$w
    v <- params$v
    sumN <- sum(state)
    # sum_wv = sum( (w_j+v_j)*N_j )
    sum_wv <- 0
    for (j in seq_len(n_comp)) {
        sum_wv <- sum_wv + (w[j] + v[j]) * state[j]
    }
    dN <- numeric(n_comp)
    dN[1] <- w[1] * state[1] - (state[1] / sumN) * sum_wv
    for (j in 2:n_comp) {
        dN[j] <- w[j] * state[j] + v[j - 1] * state[j - 1] - (state[j] / sumN) * sum_wv
    }
    list(dN)
}

mean_ODE_solver <- function(lambda_vec,
                            u_vec,
                            max_time,
                            dt = 0.1) {
    w <- (1 - 2 * u_vec) * lambda_vec
    v <- u_vec * lambda_vec
    params <- list(w = w, v = v)

    n_comp <- length(lambda_vec)
    state_init <- numeric(n_comp)
    state_init[1] <- 1

    times <- seq(0, max_time, by = dt)

    ode_out <- deSolve::ode(
        y     = state_init,
        times = times,
        func  = mean_ODE_equations,
        parms = params
    )

    df <- as.data.frame(ode_out)
    for (j in seq_len(n_comp)) {
        names(df)[j + 1] <- paste0("Nbar", j - 1)
    }
    df
}

variance_ODE_equations <- function(t,
                                   state,
                                   params) {
    f <- function(x) {
        output <- params$w * x + c(0, params$v[1:(params$nCompartments - 1)] * x[1:(params$nCompartments - 1)]) -
            params$alpha * x / sum(x) * sum((params$w + params$v) * x)
        return(output)
    }

    V <- matrix(state, nrow = params$nCompartments, byrow = TRUE)

    tmp <- which(params$mean_ode_df$time <= t)
    N <- as.numeric(params$mean_ode_df[tmp[length(tmp)], paste0("Nbar", 0:(params$nCompartments - 1))])

    S_diag <- (params$w + 2 * params$v) * N +
        params$alpha * sum((params$w + params$v) * N) * N / sum(N) +
        c(0, params$v[1:(params$nCompartments - 1)] * N[1:(params$nCompartments - 1)])
    S_offdiag <- -params$v[1:(params$nCompartments - 1)] * N[1:(params$nCompartments - 1)]
    S <- matrix(0, nrow = params$nCompartments, ncol = params$nCompartments)
    for (j in 1:params$nCompartments) S[j, j] <- S_diag[j]
    for (j in 1:(params$nCompartments - 1)) S[j, j + 1] <- S_offdiag[j]

    A <- jacobian(f, N)
    dVar <- A %*% V + t(A %*% V) + S
    list(dVar)
}

variance_ODE_solver <- function(mean_ode_df, lambda_vec, u_vec, alpha) {
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
        func  = variance_ODE_equations,
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
    return(res)
}

variance_ODE_solver_OLD <- function(mean_ode_df, lambda_vec, u_vec, alpha) {
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
