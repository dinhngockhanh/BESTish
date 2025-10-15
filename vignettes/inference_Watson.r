suppressPackageStartupMessages({
    library(pracma)
    library(MASS)
    library(ggplot2)
    library(dplyr)
    library(viridis)
    library(grid)
    library(gridExtra)
    library(parallel)
    library(pbapply)
    # current_wd <- getwd()
    current_wd <- "/Users/kndinh/DINH LAB/RESULTS/2025-10-09.Keito Taketomi.Inference for DNMT3A VAF trajectories/"
    path_R <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R"
    setwd(path_R)
    files_sources <- list.files(pattern = "\\.[rR]$")
    sapply(files_sources, source)
    setwd(current_wd)
})

options(stringsAsFactors = FALSE)
set.seed(2025)
#-----------------------------------------------Parameter specifications
plot_nsimulations <- 100 ##############################################
plot_tau <- 0.01 #######################################################
make_centers <- function(minv, maxv, nb) minv + ((0:(nb - 1)) + 0.5) * ((maxv - minv) / nb)
w1_min <- 1
w1_max <- 1.3
w1_nbins <- 5 ########################################################
log10v0_min <- -7
log10v0_max <- -4
log10v0_nbins <- 5 ###################################################
alpha_min <- 0.5
alpha_max <- 1
alpha_nbins <- 5 #####################################################
w0 <- 1
R <- 1e5
grid <- expand.grid(
    w1 = make_centers(w1_min, w1_max, w1_nbins),
    log10v0 = make_centers(log10v0_min, log10v0_max, log10v0_nbins),
    alpha = make_centers(alpha_min, alpha_max, alpha_nbins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
grid$v0 <- 10^grid$log10v0
#------------------------------------------------------Input Watson data
path <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/data/watson.csv"
watson <- read.csv(path, stringsAsFactors = FALSE)
watson$Age <- as.numeric(watson$Age)
watson$VAF <- as.numeric(watson$VAF)
watson <- watson[order(watson$Age), ]
watson <- subset(watson, is.finite(Age) & is.finite(VAF))
patient_age <- watson$Age
patient_vaf <- watson$VAF
outdir <- "Results_Watson"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#----------------------------------------------------Necessary functions
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

lik_log_analytic_Watson <- function(param, patient_age, patient_vaf) {
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
# #-------------------------------------Compute log posterior distribution
# timing <- system.time({
#     #---Parallel version
#     cl <- makePSOCKcluster(detectCores() - 1)
#     clusterExport(
#         cl,
#         varlist = c(
#             "patient_age", "patient_vaf",
#             "grid", "w0", "R",
#             "lik_log_analytic_Watson", "Jg", "A_analytic", "compute_expected_means_and_variances"
#         ),
#         envir = environment()
#     )
#     loglikelihoods <- parLapply(
#         cl = cl, 1:nrow(grid),
#         function(i) {
#             p <- c(w0, grid$w1[i], grid$v0[i], grid$alpha[i], R)
#             loglikelihood <- lik_log_analytic_Watson(p, patient_age, patient_vaf)
#         }
#     )
#     stopCluster(cl)
#     grid$loglikelihood <- unlist(loglikelihoods)
#     # #---Sequential version
#     # grid$loglikelihood <- NA_real_
#     # for (i in seq_len(nrow(grid))) {
#     #     p <- c(w0, grid$w1[i], grid$v0[i], grid$alpha[i], R)
#     #     print(p)
#     #     loglikelihood <- tryCatch(lik_log_analytic_Watson(p, patient_age, patient_vaf),
#     #         error = function(e) {
#     #             message("lik_log_analytic_Watson error at row ", i, ": ", e$message)
#     #             -Inf
#     #         }
#     #     )
#     #     grid$loglikelihood[i] <- as.numeric(loglikelihood)
#     #     if (i %% 10000 == 0) cat(sprintf("  progress: %d / %d\n", i, nrow(grid)))
#     # }
# })
# print(timing)
# #-----------------------------------------Compute posterior distribution
# grid$posterior <- exp(grid$loglikelihood - max(grid$loglikelihood))
# #--------------------------------------------Save posterior distribution
# saveRDS(grid, file = file.path(outdir, "posterior_Watson.rds"))
# write.csv(grid,
#     file = file.path(outdir, "posterior_Watson.csv"),
#     row.names = FALSE
# )
########################################################################
########################################################################
########################################################################
grid <- readRDS(file.path(outdir, "posterior_Watson.rds"))
#-------------------------------------Find MAP in posterior distribution
#---Get MAP parameters from posterior distribution
map_idx <- which.max(grid$posterior)
map_row <- grid[map_idx, ]
#---Simulate trajectories with MAP parameters
lambda0 <- w0 + 2 * map_row$v0
lambda1 <- map_row$w1
u0 <- map_row$v0 / lambda0
u1 <- 0
map_simulations <- run_replicates(
    n_reps = plot_nsimulations,
    r = R,
    lambda_vec = c(lambda0, lambda1),
    u_vec = c(u0, u1),
    alpha = map_row$alpha,
    max_time = max(patient_age, na.rm = TRUE) + 1,
    tau = plot_tau,
    seed = 123
)
map_simulations$vaf <- 0.5 * map_simulations$N1 / (map_simulations$N0 + map_simulations$N1)
map_simulations_mean <- map_simulations %>%
    group_by(time) %>%
    summarise(across(where(is.numeric) & !matches("replicate"), mean, na.rm = TRUE))
#---Compute mean and variance trajectories with MAP parameters
map_mean_var_df <- compute_expected_means_and_variances(
    time_grids = seq(0, max(patient_age, na.rm = TRUE) + 1, by = plot_tau),
    h = plot_tau,
    w0 = w0,
    w1 = map_row$w1,
    v0 = map_row$v0,
    alpha = map_row$alpha
)
map_mean_var_df$vaf_var <- NA_real_
for (i in 1:nrow(map_mean_var_df)) {
    ti <- map_mean_var_df$time[i]
    x <- c(map_mean_var_df$N0_bar[i], map_mean_var_df$N1_bar[i])
    Vt <- matrix(c(
        map_mean_var_df$V11[i], map_mean_var_df$V12[i],
        map_mean_var_df$V12[i], map_mean_var_df$V22[i]
    ), nrow = 2, byrow = TRUE)
    J <- Jg(x[1], x[2])
    map_mean_var_df$vaf_var[i] <- as.numeric(J %*% Vt %*% t(J)) / R
}
map_mean_var_df$vaf_lowerCI <- map_mean_var_df$vaf_mu - 1.96 * sqrt(map_mean_var_df$vaf_var)
map_mean_var_df$vaf_upperCI <- map_mean_var_df$vaf_mu + 1.96 * sqrt(map_mean_var_df$vaf_var)
#-----------------------------------------------------------Plot results
FORCE_WHITE <- theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white")
)
#---Plot marginal distributions
plot_marginal <- function(df, col) {
    df <- df %>%
        dplyr::group_by_at(col) %>%
        dplyr::summarise(prob = sum(posterior), .groups = "drop") %>%
        dplyr::arrange_at(col) %>%
        dplyr::mutate(prob = prob / sum(prob))
    ggplot(df, aes(x = !!sym(col), y = prob)) +
        geom_area(alpha = 0.3) +
        geom_line(linewidth = 1.05) +
        labs(x = col, y = NULL, title = NULL) +
        theme_minimal(base_size = 30) +
        FORCE_WHITE
}
p1 <- plot_marginal(grid, "w1")
p2 <- plot_marginal(grid, "log10v0")
p3 <- plot_marginal(grid, "alpha")
pdf(file.path(outdir, "Marginal_distributions_Watson.pdf"), width = 30, height = 10)
p <- gridExtra::arrangeGrob(p1, p2, p3,
    ncol = 3,
    top = "Marginal posterior distributions (Watson)"
)
grid.draw(p)
dev.off()
#---Plot joint distributions for each pair of parameters
plot_joint <- function(df, xcol, ycol) {
    grid_step <- function(x) {
        ux <- sort(unique(x))
        if (length(ux) < 2) 0.01 else median(diff(ux))
    }
    df <- df %>%
        dplyr::group_by_at(c(xcol, ycol)) %>%
        dplyr::summarise(prob = sum(posterior), .groups = "drop")
    df$v1 <- df[[xcol]]
    df$v2 <- df[[ycol]]
    w1 <- grid_step(df$v1)
    w2 <- grid_step(df$v2)
    ggplot(df, aes(x = v1, y = v2)) +
        geom_tile(aes(fill = prob), width = w1, height = w2) +
        scale_fill_viridis(option = "C", name = "Posterior") +
        labs(x = xcol, y = ycol, title = NULL) +
        theme_minimal(base_size = 30) +
        FORCE_WHITE
}
p12 <- plot_joint(grid, "w1", "log10v0")
p13 <- plot_joint(grid, "w1", "alpha")
p23 <- plot_joint(grid, "log10v0", "alpha")
pdf(file.path(outdir, "Joint_distributions_Watson.pdf"), width = 30, height = 10)
p <- gridExtra::arrangeGrob(p12, p13, p23,
    ncol = 3,
    top = "Joint posterior distributions (Watson)"
)
grid.draw(p)
dev.off()
#---Plot VAF dynamics
p <- ggplot() +
    geom_ribbon(
        data = map_mean_var_df,
        aes(x = time, ymin = vaf_lowerCI, ymax = vaf_upperCI),
        fill = "#BC3C29", alpha = 0.2
    ) +
    geom_line(data = map_simulations, aes(x = time, y = vaf, group = replicate), color = "#BC3C29", size = 1, alpha = 0.2) +
    geom_line(data = map_mean_var_df, aes(x = time, y = vaf_mu), color = "#BC3C29", size = 5) +
    # geom_line(data = map_simulations_mean, aes(x = time, y = vaf), color = "#0072B5", linetype = "dotted", size = 5) +
    geom_point(data = watson, aes(x = Age, y = VAF), color = "white", size = 16) +
    geom_point(data = watson, aes(x = Age, y = VAF), color = "#BC3C29", size = 13) +
    labs(
        x = "Age",
        y = "Variant Allele Frequency (VAF)",
        title = NULL
    ) +
    theme_minimal(base_size = 50) +
    FORCE_WHITE +
    theme(legend.position = "none")
pdf(file.path(outdir, "VAF_trajectories_Watson.pdf"), width = 30, height = 15)
print(p)
dev.off()
png(file.path(outdir, "VAF_trajectories_Watson.png"), width = 30, height = 15, units = "in", res = 600)
print(p)
dev.off()
