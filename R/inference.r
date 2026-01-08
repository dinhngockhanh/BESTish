make_BESTish_grid <- function(w1 = NULL,
                              w1_min = NULL,
                              w1_max = NULL,
                              w1_nbins = NULL,
                              log10v0 = NULL,
                              log10v0_min = NULL,
                              log10v0_max = NULL,
                              log10v0_nbins = NULL,
                              alpha = NULL,
                              alpha_min = NULL,
                              alpha_max = NULL,
                              alpha_nbins = NULL,
                              w0 = NULL,
                              w0_min = NULL,
                              w0_max = NULL,
                              w0_nbins = NULL,
                              R = NULL,
                              R_min = NULL,
                              R_max = NULL,
                              R_nbins = NULL) {
    make_centers <- function(minv, maxv, nb) minv + ((0:(nb - 1)) + 0.5) * ((maxv - minv) / nb)
    if (is.null(w0)) {
        if (is.null(w0_min) || is.null(w0_max) || is.null(w0_nbins)) {
            stop("w0_min, w0_max, and w0_nbins must all be provided when w0 is NULL")
        }
        w0 <- make_centers(w0_min, w0_max, w0_nbins)
    }
    if (is.null(w1)) {
        if (is.null(w1_min) || is.null(w1_max) || is.null(w1_nbins)) {
            stop("w1_min, w1_max, and w1_nbins must all be provided when w1 is NULL")
        }
        w1 <- make_centers(w1_min, w1_max, w1_nbins)
    }
    if (is.null(log10v0)) {
        if (is.null(log10v0_min) || is.null(log10v0_max) || is.null(log10v0_nbins)) {
            stop("log10v0_min, log10v0_max, and log10v0_nbins must all be provided when log10v0 is NULL")
        }
        log10v0 <- make_centers(log10v0_min, log10v0_max, log10v0_nbins)
    }
    if (is.null(alpha)) {
        if (is.null(alpha_min) || is.null(alpha_max) || is.null(alpha_nbins)) {
            stop("alpha_min, alpha_max, and alpha_nbins must all be provided when alpha is NULL")
        }
        alpha <- make_centers(alpha_min, alpha_max, alpha_nbins)
    }
    if (is.null(R)) {
        if (is.null(R_min) || is.null(R_max) || is.null(R_nbins)) {
            stop("R_min, R_max, and R_nbins must all be provided when R is NULL")
        }
        R <- make_centers(R_min, R_max, R_nbins)
    }
    BESTish_grid <- expand.grid(
        w0 = w0,
        w1 = w1,
        log10v0 = log10v0,
        alpha = alpha,
        R = R,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    BESTish_grid$v0 <- 10^BESTish_grid$log10v0
    return(BESTish_grid)
}

BESTish_inference <- function(data,
                              grid = NULL,
                              mode = NULL,
                              time_step = 0.0005,
                              parallel = FALSE,
                              progress = TRUE) {
    suppressPackageStartupMessages(library(crayon))
    if (!is.data.frame(data)) {
        stop("data must be a data frame")
    }
    if (!("Age" %in% names(data)) || !("VAF" %in% names(data))) {
        stop("data must contain 'Age' and 'VAF' columns")
    }
    if (is.null(grid)) {
        stop("BESTish grid must be provided. Check function make_BESTish_grid() for more information.")
    }
    if (!is.null(mode) && !(mode %in% c("longitudinal", "cohort"))) {
        stop("mode must be either 'longitudinal' or 'cohort'")
    }
    data <- data[order(data$Age), ]
    #---Compute log-likelihood at each grid point
    if (parallel) {
        suppressPackageStartupMessages(library(parallel))
        suppressPackageStartupMessages(library(pbapply))
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(
            cl,
            varlist = c(
                "loglikelihood_longitudinal", "loglikelihood_cohort",
                "compute_expected_means_and_variances", "A_analytic", "Jg"
            ),
            envir = environment()
        )
        if (mode == "cohort") {
            cat(bold(red("BESTish for cohort data in parallel computing mode...\n")))
        } else if (mode == "longitudinal") {
            cat(bold(red("BESTish for longitudinal data in parallel computing mode...\n")))
        }
        if (mode == "cohort") {
            if (progress) {
                loglikelihoods <- pblapply(
                    cl = cl, X = 1:nrow(grid),
                    FUN = function(i) {
                        loglikelihood_cohort(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        )
                    }
                )
            } else {
                loglikelihoods <- parLapply(
                    cl = cl, 1:nrow(grid),
                    function(i) {
                        loglikelihood_cohort(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        )
                    }
                )
            }
        } else if (mode == "longitudinal") {
            if (progress) {
                loglikelihoods <- pblapply(
                    cl = cl, X = 1:nrow(grid),
                    FUN = function(i) {
                        loglikelihood_longitudinal(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        )
                    }
                )
            } else {
                loglikelihoods <- parLapply(
                    cl = cl, 1:nrow(grid),
                    function(i) {
                        loglikelihood_longitudinal(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        )
                    }
                )
            }
        }
        stopCluster(cl)
        grid$loglikelihood <- unlist(loglikelihoods)
    } else {
        grid$loglikelihood <- NA_real_
        if (mode == "cohort") {
            cat(bold(red("BESTish for cohort data in sequential computing mode...\n")))
        } else if (mode == "longitudinal") {
            cat(bold(red("BESTish for longitudinal data in sequential computing mode...\n")))
        }
        if (progress) {
            pb <- txtProgressBar(
                min = 1, max = nrow(grid),
                style = 3, width = 50, char = "+"
            )
        }
        if (mode == "cohort") {
            for (i in seq_len(nrow(grid))) {
                if (progress) {
                    setTxtProgressBar(pb, i)
                    grid$loglikelihood[i] <- tryCatch(
                        loglikelihood_cohort(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        ),
                        error = function(e) {
                            message("loglikelihood_cohort error at row ", i, " of BESTish grid system: ", e$message)
                            -Inf
                        }
                    )
                }
            }
        } else if (mode == "longitudinal") {
            for (i in seq_len(nrow(grid))) {
                if (progress) {
                    setTxtProgressBar(pb, i)
                    grid$loglikelihood[i] <- tryCatch(
                        loglikelihood_longitudinal(
                            parameter = grid[i, ],
                            data = data,
                            time_step = time_step
                        ),
                        error = function(e) {
                            message("loglikelihood_longitudinal error at row ", i, " of BESTish grid system: ", e$message)
                            -Inf
                        }
                    )
                }
            }
        }
        if (progress == TRUE) {
            cat("\n")
        }
    }
    #---Compute scaled posterior distribution across whole grid
    grid$loglikelihood <- as.numeric(grid$loglikelihood)
    grid$posterior <- exp(grid$loglikelihood - max(grid$loglikelihood, na.rm = TRUE))
    return(grid)
}

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

loglikelihood_cohort <- function(parameter,
                                 data,
                                 time_step) {
    w0 <- parameter$w0
    w1 <- parameter$w1
    v0 <- parameter$v0
    alpha <- parameter$alpha
    R <- parameter$R
    if (abs(w1 - w0) < .Machine$double.eps) {
        return(-Inf)
    }
    #---Create time grid from 0 to max age + 1 with input step size
    time_grids <- seq(0, max(data$Age, na.rm = TRUE) + 1, by = time_step)
    obs_indices <- pmax(1, findInterval(data$Age, time_grids))
    #---Compute closed-form mean and variance trajectory
    mean_df <- compute_expected_means_and_variances(time_grids, time_step, w0, w1, v0, alpha)
    N0_bar <- mean_df$N0_bar
    N1_bar <- mean_df$N1_bar
    vaf_mu <- mean_df$vaf_mu
    V11 <- mean_df$V11
    V12 <- mean_df$V12
    V22 <- mean_df$V22
    #---Find predicted means at the observed age points
    mu <- vaf_mu[obs_indices]
    #---Compute log-likelihood of observed VAFs under independent normal model
    n <- length(data$Age)
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
        resid <- data$VAF[i] - mu[i]
        ll <- ll + (-0.5) * (log(2 * pi * var_i) + (resid * resid) / var_i)
    }
    return(as.numeric(ll))
}

loglikelihood_longitudinal <- function(parameter,
                                       data,
                                       time_step) {
    suppressPackageStartupMessages(library(mvtnorm))
    suppressPackageStartupMessages(library(Matrix))
    w0 <- parameter$w0
    w1 <- parameter$w1
    v0 <- parameter$v0
    alpha <- parameter$alpha
    R <- parameter$R
    if (abs(w1 - w0) < .Machine$double.eps) {
        return(-Inf)
    }
    #---Create time grid from 0 to max age + 1 year with input step size
    time_grids <- seq(0, max(data$Age, na.rm = TRUE) + 1, by = time_step)
    obs_indices <- pmax(1, findInterval(data$Age, time_grids))
    #---Compute closed-form mean and variance trajectory
    mean_df <- compute_expected_means_and_variances(time_grids, time_step, w0, w1, v0, alpha)
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
            Phi[[k]] <- Phi[[k - 1]] + A_mid %*% Phi[[k - 1]] * time_step
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
    n <- length(data$Age)
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
    return(dmvnorm(data$VAF, mean = mu, sigma = S, log = TRUE))
}
