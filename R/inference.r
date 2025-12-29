Jg <- function(N0, N1) {
    S <- N0 + N1
    matrix(c(
        -0.5 * N1 / (S * S),
        0.5 * N0 / (S * S)
    ), nrow = 1, ncol = 2)
}

A_analytic <- function(N0, N1, w0, w1, v0, alpha) {
    N0 <- max(0, N0)
    N1 <- max(0, N1)
    S <- N0 + N1
    a <- w0 + v0
    b <- w1
    ttl <- a * N0 + b * N1

    # derivatives of terms inside the crowding component
    dT0_dN0 <- (ttl * N1 + a * N0 * S) / (S * S)
    dT0_dN1 <- N0 * (b * S - ttl) / (S * S)
    dT1_dN0 <- N1 * (a * S - ttl) / (S * S)
    dT1_dN1 <- (ttl * N0 + b * N1 * S) / (S * S)

    A11 <- w0 - alpha * dT0_dN0
    A12 <- -alpha * dT0_dN1
    A21 <- v0 - alpha * dT1_dN0
    A22 <- w1 - alpha * dT1_dN1

    matrix(c(
        A11, A12,
        A21, A22
    ), nrow = 2, byrow = TRUE)
}

compute_expected_means_and_variances <- function(time_grids, h, w0, w1, v0, alpha) {
    #---Compute expected means
    # delta <- (w1 - w0) * time_grids
    # cval <- (v0 / (w1 - w0)) * expm1(delta)
    # N0_bar <- exp((1 - alpha) * w0 * time_grids) * (1 + cval)^(-alpha)
    # N1_bar <- N0_bar * cval
    # N0_bar[!is.finite(N0_bar)] <- 0
    # N1_bar[!is.finite(N1_bar)] <- 0
    # N1_bar <- pmax(0, N1_bar)
    # vaf_mu <- 0.5 * N1_bar / (N0_bar + N1_bar)

    if (w1 > w0) {
        delta <- (w1 - w0) * time_grids
        cval <- (v0 / (w1 - w0)) * expm1(delta)
        N0_bar <- exp((1 - alpha) * w0 * time_grids) * (1 + cval)^(-alpha)
        N1_bar <- N0_bar * cval
        N0_bar[!is.finite(N0_bar)] <- 0
        N1_bar[!is.finite(N1_bar)] <- 0
        N1_bar <- pmax(0, N1_bar)
    } else if (w1 == w0) {
        exponent <- (1 - alpha) * w0 * time_grids
        N0_bar <- (1 + v0 * time_grids)^(-alpha) * exp(exponent)
        N1_bar <- v0 * time_grids * N0_bar
    } else {
        stop("w1 must be greater than or equal to w0")
    }
    vaf_mu <- 0.5 * N1_bar / (N0_bar + N1_bar)

    #---Compute expected variances
    V11 <- V12 <- V22 <- numeric(length(time_grids))

    # Pre-compute constants
    a <- w0 + v0
    b <- w1

    # Vectorize preparatory calculations for all time points
    n <- length(time_grids) - 1
    n0_vec <- pmax(0, N0_bar[1:n])
    n1_vec <- pmax(0, N1_bar[1:n])
    sumN_vec <- n0_vec + n1_vec

    # Create mask for valid points (avoid division by zero)
    valid <- sumN_vec >= .Machine$double.eps

    # Pre-compute vectorized intermediate values
    ttl_vec <- a * n0_vec + b * n1_vec
    inv_sumN_vec <- ifelse(valid, 1 / sumN_vec, 0)
    inv_sumN2_vec <- inv_sumN_vec * inv_sumN_vec

    # Vectorized A matrix elements
    dT0_dN0_vec <- (ttl_vec * n1_vec + a * n0_vec * sumN_vec) * inv_sumN2_vec
    dT0_dN1_vec <- n0_vec * (b * sumN_vec - ttl_vec) * inv_sumN2_vec
    dT1_dN0_vec <- n1_vec * (a * sumN_vec - ttl_vec) * inv_sumN2_vec
    dT1_dN1_vec <- (ttl_vec * n0_vec + b * n1_vec * sumN_vec) * inv_sumN2_vec

    A11_vec <- w0 - alpha * dT0_dN0_vec
    A12_vec <- -alpha * dT0_dN1_vec
    A21_vec <- v0 - alpha * dT1_dN0_vec
    A22_vec <- w1 - alpha * dT1_dN1_vec

    # Vectorized Sigma matrix elements
    Sigma11_vec <- (w0 + 2 * v0) * n0_vec + alpha * n0_vec * ttl_vec * inv_sumN_vec
    Sigma12_vec <- -v0 * n0_vec
    Sigma22_vec <- w1 * n1_vec + v0 * n0_vec + alpha * n1_vec * ttl_vec * inv_sumN_vec

    # Loop is still needed for sequential dependency, but all preparatory work is vectorized
    for (i in 1:n) {
        if (!valid[i]) next

        v11 <- V11[i]
        v12 <- V12[i]
        v22 <- V22[i]

        # Compute A*V using pre-computed A elements
        AV11 <- A11_vec[i] * v11 + A12_vec[i] * v12
        AV12 <- A11_vec[i] * v12 + A12_vec[i] * v22
        AV21 <- A21_vec[i] * v11 + A22_vec[i] * v12
        AV22 <- A21_vec[i] * v12 + A22_vec[i] * v22

        # upd = A*V + (A*V)^T + Sigma (using pre-computed Sigma)
        upd11 <- 2 * AV11 + Sigma11_vec[i]
        upd12 <- AV12 + AV21 + Sigma12_vec[i]
        upd22 <- 2 * AV22 + Sigma22_vec[i]

        # Update V values
        V11[i + 1] <- v11 + upd11 * h
        V12[i + 1] <- v12 + upd12 * h
        V22[i + 1] <- v22 + upd22 * h
    }

    #---Return expected means and variances
    return(data.frame(
        time = time_grids,
        N0_bar = N0_bar,
        N1_bar = N1_bar,
        vaf_mu = vaf_mu,
        V11 = V11,
        V12 = V12,
        V22 = V22
    ))
}

loglikelihood_cohort <- function(param,
                                 patient_age,
                                 patient_vaf,
                                 h = 0.01) {
    w0 <- param[1]
    w1 <- param[2]
    v0 <- param[3]
    alpha <- param[4]
    R <- param[5]
    if (abs(w1 - w0) < .Machine$double.eps) {
        return(-Inf)
    }
    #---Create time grid from 0 to max age + 1 with input step size
    time_grids <- seq(0, max(patient_age, na.rm = TRUE) + 1, by = h)
    obs_indices <- pmax(1, findInterval(patient_age, time_grids))
    #---Compute closed-form mean and variance trajectory
    mean_df <- compute_expected_means_and_variances(time_grids, h, w0, w1, v0, alpha)
    N0_bar <- mean_df$N0_bar
    N1_bar <- mean_df$N1_bar
    vaf_mu <- mean_df$vaf_mu
    V11 <- mean_df$V11
    V12 <- mean_df$V12
    V22 <- mean_df$V22
    #---Find predicted means at the observed age points
    mu <- vaf_mu[obs_indices]
    #---Compute log-likelihood of observed VAFs under independent normal model
    n <- length(patient_age)
    ll <- 0.0
    for (i in 1:n) {
        ti <- obs_indices[i]
        x <- c(N0_bar[ti], N1_bar[ti])
        Vt <- matrix(c(
            V11[ti], V12[ti],
            V12[ti], V22[ti]
        ), nrow = 2, byrow = TRUE)
        J <- Jg(x[1], x[2])
        var_i <- as.numeric(J %*% Vt %*% t(J)) / R
        resid <- patient_vaf[i] - mu[i]
        ll <- ll + (-0.5) * (log(2 * pi * var_i) + (resid * resid) / var_i)
    }
    return(as.numeric(ll))
}

loglikelihood_timeseries <- function(param,
                                     patient_age,
                                     patient_vaf,
                                     h = 0.0005) {
    library(mvtnorm)
    library(Matrix)
    w0 <- param[1]
    w1 <- param[2]
    v0 <- param[3]
    alpha <- param[4]
    R <- param[5]
    if (abs(w1 - w0) < .Machine$double.eps) {
        return(-Inf)
    }
    #---Create time grid from 0 to max age + 1 with input step size
    time_grids <- seq(0, max(patient_age, na.rm = TRUE) + 1, by = h)
    obs_indices <- pmax(1, findInterval(patient_age, time_grids))
    #---Compute closed-form mean and variance trajectory
    mean_df <- compute_expected_means_and_variances(time_grids, h, w0, w1, v0, alpha)
    N0_bar <- mean_df$N0_bar
    N1_bar <- mean_df$N1_bar
    vaf_mu <- mean_df$vaf_mu
    V11 <- mean_df$V11
    V12 <- mean_df$V12
    V22 <- mean_df$V22
    #---Find predicted means at the observed age points
    mu <- vaf_mu[obs_indices]
    #---Compute Phi at observed age points
    Phi <- vector("list", max(obs_indices))
    Phi[[1]] <- diag(2)
    last_computed <- 1
    for (target_idx in obs_indices[-1]) {
        for (k in (last_computed + 1):target_idx) {
            A_mid <- A_analytic(
                (N0_bar[k - 1] + N0_bar[k]) / 2,
                (N1_bar[k - 1] + N1_bar[k]) / 2,
                w0, w1, v0, alpha
            )
            Phi[[k]] <- Phi[[k - 1]] + A_mid %*% Phi[[k - 1]] * h
        }
        last_computed <- target_idx
    }
    #---Define autocovariance kernel K(s,t)
    Kst <- function(si, ti) {
        Vs <- matrix(c(
            V11[si], V12[si],
            V12[si], V22[si]
        ), nrow = 2, byrow = TRUE)
        Js <- Jg(N0_bar[si], N1_bar[si])
        Jt <- Jg(N0_bar[ti], N1_bar[ti])
        Ps <- Phi[[si]]
        Pt <- Phi[[ti]]
        # More stable computation: solve(A, B) instead of solve(A) %*% B
        # Compute: Js %*% Vs %*% inv(t(Ps)) %*% t(Pt) %*% t(Jt)
        # Rewrite as: Js %*% Vs %*% solve(t(Ps), t(Pt) %*% t(Jt))
        temp <- solve(t(Ps), t(Pt) %*% t(Jt))
        as.numeric(Js %*% Vs %*% temp)
    }
    #---Find predicted covariance matrix at the observed age points
    n <- length(patient_age)
    S <- matrix(0, n, n)
    for (i in 1:n) {
        for (j in i:n) {
            kij <- Kst(obs_indices[i], obs_indices[j]) / R
            S[i, j] <- kij
            S[j, i] <- kij
        }
    }
    S <- (S + t(S)) / 2
    #---Return log-likelihood of observed VAFs under multivariate normal model
    S <- as.matrix(nearPD(S, corr = FALSE)$mat)
    return(dmvnorm(patient_vaf, mean = mu, sigma = S, log = TRUE))
}
