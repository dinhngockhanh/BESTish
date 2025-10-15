suppressPackageStartupMessages({
    need <- c("pracma", "MASS", "ggplot2", "dplyr", "viridis", "gridExtra")
    have <- rownames(installed.packages())
    miss <- setdiff(need, have)
    if (length(miss) > 0) {
        stop(
            "Missing R packages: ", paste(miss, collapse = ", "),
            "\nInstall to R_LIBS_USER: ", Sys.getenv("R_LIBS_USER")
        )
    }
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
    current_wd <- "/Users/kndinh/Downloads/2025-10-05.Keito Taketomi"
    path_R <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R"
    setwd(path_R)
    files_sources <- list.files(pattern = "\\.[rR]$")
    sapply(files_sources, source)
    setwd(current_wd)
})

options(stringsAsFactors = FALSE)
set.seed(2025)
#-----------------------------------------------Parameter specifications
patient_ID <- "4240" ###################################################
plot_nsimulations <- 20 ###############################################
plot_tau <- 0.001 #######################################################
make_centers <- function(minv, maxv, nb) minv + ((0:(nb - 1)) + 0.5) * ((maxv - minv) / nb)
w1_min <- 1
w1_max <- 1.3
w1_nbins <- 20 ########################################################
log10v0_min <- -7
log10v0_max <- -4
log10v0_nbins <- 20 ###################################################
alpha <- 1 #############################################################
w0 <- 1
R <- 1e5
grid <- expand.grid(
    w1 = make_centers(w1_min, w1_max, w1_nbins),
    log10v0 = make_centers(log10v0_min, log10v0_max, log10v0_nbins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
grid$v0 <- 10^grid$log10v0
#-------------------------------------------------------Input Fabre data
path <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/data/ALLvariants_exclSynonymous_Xadj.txt"
fabre <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
fabre$Age <- as.numeric(fabre$Age)
fabre$VAF <- as.numeric(fabre$VAF)
fabre <- subset(
    fabre,
    Gene == "DNMT3A" &
        SardID == patient_ID &
        grepl("\\bR882H\\b", AAChange.refGene, ignore.case = TRUE)
)
fabre <- fabre[order(fabre$Age), ]
patient_age <- fabre$Age
patient_vaf <- fabre$VAF
outdir <- paste0("Results_Fabre_", patient_ID)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#----------------------------------------------------Necessary functions
Jg <- function(N0, N1) {
    eps <- .Machine$double.eps

    S <- N0 + N1 + eps
    matrix(c(
        -0.5 * N1 / (S * S),
        0.5 * N0 / (S * S)
    ), nrow = 1, ncol = 2)
}

A_analytic <- function(N0, N1, w0, w1, v0, alpha) {
    eps <- .Machine$double.eps

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
    eps <- .Machine$double.eps
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

lik_log_analytic_Fabre <- function(param, patient_age, patient_vaf) {
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
# print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# print(lik_log_analytic_Fabre(
#     c(w0, 1.1155, 10^-4.615, alpha, R),
#     patient_age, patient_vaf
# ))
# print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# #-------------------------------------Compute log posterior distribution
# timing <- system.time({
#     #---Parallel version
#     cl <- makePSOCKcluster(detectCores() - 1)
#     clusterExport(
#         cl,
#         varlist = c(
#             "patient_age", "patient_vaf",
#             "grid", "w0", "R", "alpha",
#             "lik_log_analytic_Fabre", "Jg", "A_analytic", "compute_expected_means_and_variances"
#         ),
#         envir = environment()
#     )
#     loglikelihoods <- parLapply(
#         cl = cl, 1:nrow(grid),
#         function(i) {
#             p <- c(w0, grid$w1[i], grid$v0[i], alpha, R)
#             loglikelihood <- lik_log_analytic_Fabre(p, patient_age, patient_vaf)
#         }
#     )
#     stopCluster(cl)
#     grid$loglikelihood <- unlist(loglikelihoods)
#     # #---Sequential version
#     # grid$loglikelihood <- NA_real_
#     # for (i in seq_len(nrow(grid))) {
#     #     p <- c(w0, grid$w1[i], grid$v0[i], alpha, R)
#     #     print(p)
#     #     loglikelihood <- tryCatch(lik_log_analytic_Fabre(p, patient_age, patient_vaf),
#     #         error = function(e) {
#     #             message("lik_log_analytic_Fabre error at row ", i, ": ", e$message)
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
# saveRDS(grid, file = file.path(outdir, paste0("posterior_Fabre_", patient_ID, ".rds")))
# write.csv(grid,
#     file = file.path(outdir, paste0("posterior_Fabre_", patient_ID, ".csv")),
#     row.names = FALSE
# )
########################################################################
########################################################################
########################################################################
outdir <- "/Users/kndinh/DINH LAB/RESULTS/2025-10-09.Keito Taketomi.Inference for DNMT3A VAF trajectories/Results_Fabre_4240/"


grid <- readRDS(file.path(outdir, paste0("posterior_Fabre_", patient_ID, ".rds")))
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
    alpha = alpha,
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
    alpha = alpha
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
pdf(file.path(outdir, paste0("Marginal_distributions_Fabre_", patient_ID, ".pdf")), width = 20, height = 10)
p <- gridExtra::arrangeGrob(p1, p2,
    ncol = 2,
    top = paste0("Marginal posterior distributions (Fabre - ", patient_ID, ")")
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
p <- plot_joint(grid, "w1", "log10v0")
pdf(file.path(outdir, paste0("Joint_distributions_Fabre_", patient_ID, ".pdf")), width = 10, height = 10)
p <- gridExtra::arrangeGrob(p,
    ncol = 1,
    top = paste0("Joint posterior distributions (Fabre - ", patient_ID, ")")
)
grid.draw(p)
dev.off()
#---Plot VAF dynamics
p <- ggplot() +
    geom_ribbon(
        data = map_mean_var_df,
        aes(x = time, ymin = vaf_lowerCI, ymax = vaf_upperCI),
        fill = "#BC3C29", alpha = 0.3
    ) +
    geom_line(data = map_simulations, aes(x = time, y = vaf, group = replicate), color = "#0072B5", size = 1, alpha = 0.2) +
    geom_line(data = map_mean_var_df, aes(x = time, y = vaf_mu), color = "#BC3C29", size = 5) +
    geom_line(data = map_simulations_mean, aes(x = time, y = vaf), color = "#0072B5", linetype = "dotted", size = 5) +
    geom_point(data = fabre, aes(x = Age, y = VAF), color = "white", size = 14) +
    geom_point(data = fabre, aes(x = Age, y = VAF), color = "#00A087", size = 10) +
    labs(
        x = "Age",
        y = "Variant Allele Frequency (VAF)",
        title = NULL
    ) +
    theme_minimal(base_size = 50) +
    FORCE_WHITE +
    theme(legend.position = "none")
pdf(file.path(outdir, paste0("VAF_trajectories_Fabre_", patient_ID, ".pdf")), width = 30, height = 15)
print(p)
dev.off()
