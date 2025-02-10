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
    state <- numeric(n_types)
    state[1] <- r
    t <- 0
    record_time <- c(t)
    record_state <- matrix(state, nrow = 1)
    #---Tau-leaping algorithm
    while (t < max_time) {
        N_total <- sum(state)
        if (N_total <= 0) break
        #   Simulate number of divisions into same compartments
        divisions <- rpois(n_types, (1 - u_vec) * lambda_vec * state * tau)
        #   Simulate number of mutations into new compartments
        mutations <- rpois(n_types, u_vec * lambda_vec * state * tau)
        #   Simulate death counts
        deaths <- rpois(n_types, alpha * sum((1 - u_vec) * lambda_vec * state) * state / sum(state) * tau)
        #   Update time
        t <- t + tau
        #   Update compartment cell counts
        state <- state - deaths + divisions - mutations +
            c(0, mutations[1:(length(mutations) - 1)])
        state[state < 0] <- 0
        #   Record state and time
        record_time <- c(record_time, t)
        record_state <- rbind(record_state, state)
    }
    #   Output system states and times
    out_df <- data.frame(time = record_time)
    for (j in seq_len(n_types)) {
        out_df[[paste0("N", j - 1)]] <- record_state[, j]
        out_df[[paste0("Nbar", j - 1)]] <- record_state[, j] / r
    }
    return(out_df)
}
