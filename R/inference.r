Jg <- function(N0, N1) {
    eps <- 1e-15

    S <- N0 + N1 + eps
    matrix(c(
        -0.5 * N1 / (S * S),
        0.5 * N0 / (S * S)
    ), nrow = 1, ncol = 2)
}

A_analytic <- function(N0, N1, w0, w1, v0, alpha) {
    eps <- 1e-15

    N0 <- max(eps, N0)
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
    eps <- 1e-15
    #---Compute expected means
    delta <- (w1 - w0) * time_grids
    cval <- (v0 / (w1 - w0)) * expm1(delta)
    N0_bar <- exp((1 - alpha) * w0 * time_grids) * (1 + cval)^(-alpha)
    N1_bar <- N0_bar * cval
    N0_bar[!is.finite(N0_bar)] <- eps
    N1_bar[!is.finite(N1_bar)] <- 0
    N1_bar <- pmax(0, N1_bar)
    vaf_mu <- 0.5 * N1_bar / pmax(eps, N0_bar + N1_bar)
    #---Compute expected variances
    V11 <- V12 <- V22 <- numeric(length(time_grids))
    for (i in 1:(length(time_grids) - 1)) {
        n0 <- max(eps, N0_bar[i])
        n1 <- max(0, N1_bar[i])
        sumN <- n0 + n1
        ttl <- (w0 + v0) * n0 + w1 * n1
        Sigma <- matrix(c(
            (w0 + 2 * v0) * n0 + alpha * n0 * ttl / sumN,  -v0 * n0,
            -v0 * n0,                                      w1 * n1 + v0 * n0 + alpha * n1 * ttl / sumN
        ), 2, 2, byrow = TRUE)
        A <- A_analytic(n0, n1, w0, w1, v0, alpha)
        V <- matrix(c(
            V11[i], V12[i],
            V12[i], V22[i]
        ), nrow = 2, byrow = TRUE)
        upd <- A %*% V + t(A %*% V) + Sigma # = A V + V A^T + Sigma
        V11[i + 1] <- V11[i] + upd[1, 1] * h
        V12[i + 1] <- V12[i] + upd[1, 2] * h
        V22[i + 1] <- V22[i] + upd[2, 2] * h
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

loglikelihood_cohort <- function(param, patient_age, patient_vaf) {
    eps <- 1e-15

    w0 <- param[1]
    w1 <- param[2]
    v0 <- param[3]
    alpha <- param[4]
    R <- param[5]

    # guard (as in Ren-Yi; here log-lik => -Inf)
    if (abs(w1 - w0) < 1e-12) {
        return(-Inf)
    }

    # ---- time grid (0 .. max(age)+1) ----
    h <- 0.01
    tmax <- max(patient_age, na.rm = TRUE) + 1
    time_grids <- seq(0, tmax, by = h)

    # ---- mean trajectory (closed form) ----
    mean_df <- compute_expected_means_and_variances(time_grids, h, w0, w1, v0, alpha)
    N0_bar <- mean_df$N0_bar
    N1_bar <- mean_df$N1_bar
    vaf_mu <- mean_df$vaf_mu
    V11 <- mean_df$V11
    V12 <- mean_df$V12
    V22 <- mean_df$V22

    # Predicted means at the observed ages
    mu <- vapply(patient_age, function(a) {
        idx <- max(1, sum(time_grids <= a))
        vaf_mu[idx]
    }, numeric(1))

    # (No time-series covariance across different ages in cohort mode)
    # Each observation i contributes a scalar variance via the delta method
    # var_i = Jg(x_t) V_t Jg(x_t)^T / R
    n <- length(patient_age)
    ll <- 0.0
    for (i in 1:n) {
        ti <- max(1, sum(time_grids <= patient_age[i]))
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

    as.numeric(ll)
}

loglikelihood_timeseries <- function(param, patient_age, patient_vaf) {
    library(mvtnorm)
    library(Matrix)

    eps <- .Machine$double.eps

    w0 <- param[1]
    w1 <- param[2]
    v0 <- param[3]
    alpha <- param[4]
    R <- param[5]

    # guard (as in Ren-Yi; here log-lik => -Inf)
    if (abs(w1 - w0) < 1e-12) {
        return(-Inf)
    }
    # ---- time grid (0 .. max(age)+1) ----
    h <- 0.0001
    # h <- 0.01
    tmax <- max(patient_age, na.rm = TRUE) + 1
    time_grids <- seq(0, tmax, by = h)

    # ---- mean trajectory (closed form) ----
    mean_df <- compute_expected_means_and_variances(time_grids, h, w0, w1, v0, alpha)
    N0_bar <- mean_df$N0_bar
    N1_bar <- mean_df$N1_bar
    vaf_mu <- mean_df$vaf_mu
    V11 <- mean_df$V11
    V12 <- mean_df$V12
    V22 <- mean_df$V22

    # Predicted means at the observed ages
    mu <- vapply(patient_age, function(a) {
        idx <- max(1, sum(time_grids <= a))
        vaf_mu[idx]
    }, numeric(1))
    # ---- Fundamental matrix Phi (A_analytic) ----
    Phi <- vector("list", length(time_grids))
    Phi[[1]] <- diag(2)
    for (k in 2:length(time_grids)) {
        A <- A_analytic(N0_bar[k - 1], N1_bar[k - 1], w0, w1, v0, alpha)
        Phi[[k]] <- Phi[[k - 1]] + A %*% Phi[[k - 1]] * h
    }

    # ---- Autocovariance kernel K(s, t) via analytic Jg ----
    Kst <- function(s, t) {
        si <- max(1, sum(time_grids <= s))
        ti <- max(1, sum(time_grids <= t))
        Vs <- matrix(c(
            V11[si], V12[si],
            V12[si], V22[si]
        ), nrow = 2, byrow = TRUE)
        Js <- Jg(N0_bar[si], N1_bar[si])
        Jt <- Jg(N0_bar[ti], N1_bar[ti])
        Ps <- Phi[[si]]
        Pt <- Phi[[ti]]
        # No Cholesky; use solve() directly in this small 2x2 context
        as.numeric(Js %*% Vs %*% solve(t(Ps)) %*% t(Pt) %*% t(Jt))
    }

    # ---- Build observation mean/cov for given ages ----
    n <- length(patient_age)

    S <- matrix(0, n, n)
    for (i in 1:n) {
        for (j in i:n) {
            kij <- Kst(patient_age[i], patient_age[j]) / R
            S[i, j] <- kij
            S[j, i] <- kij
        }
    }

    S <- as.matrix(nearPD(S, corr = FALSE)$mat)

    return(dmvnorm(patient_vaf, mean = mu, sigma = S, log = TRUE))
}
