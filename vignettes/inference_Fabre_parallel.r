args <- commandArgs(trailingOnly = TRUE)
inference_index <- as.numeric(args[1])



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
    path_R <- "/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R"
    setwd(path_R)
    files_sources <- list.files(pattern = "\\.[rR]$")
    sapply(files_sources, source)
})
########################################################################
setwd("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")
########################################################################
options(stringsAsFactors = FALSE)
set.seed(2025)
#-----------------------------------------------Parameter specifications
plot_nsimulations <- 10 ################################################
plot_tau <- 0.001 ######################################################
make_centers <- function(minv, maxv, nb) minv + ((0:(nb - 1)) + 0.5) * ((maxv - minv) / nb)
w1_min <- 1
w1_max <- 1.3
w1_nbins <- 10 #########################################################
# w1_nbins <- 500 ########################################################
log10v0_min <- -6
log10v0_max <- -3
log10v0_nbins <- 10 ####################################################
# log10v0_nbins <- 500 ###################################################
alpha <- 1 #############################################################
w0 <- 1
R <- 11000
basic_grid <- expand.grid(
    w1 = make_centers(w1_min, w1_max, w1_nbins),
    log10v0 = make_centers(log10v0_min, log10v0_max, log10v0_nbins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
basic_grid$v0 <- 10^basic_grid$log10v0
#---------Input list of MutationID x SardID pairs in Fabre for inference
inference_list <- read.csv("fabre_inference_list_filtered.csv")
#-------------------------------------------------------Input Fabre data
fabre <- read.csv("fabre_master.csv")
fabre$Age <- as.numeric(fabre$Age)
fabre$VAF <- as.numeric(fabre$VAF)




patient_ID <- inference_list$SardID[inference_index]
mutation_ID <- inference_list$MutationID[inference_index]
fabre_short <- subset(
    fabre,
    MutationID == mutation_ID & SardID == patient_ID
)
fabre_short <- fabre_short[order(fabre_short$Age), ]
patient_age <- fabre_short$Age
patient_vaf <- fabre_short$VAF
outdir <- paste0("Results_", mutation_ID, "_", patient_ID)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#---Compute log posterior distribution
grid <- basic_grid
timing <- system.time({
    #---Parallel version
    cl <- makePSOCKcluster(detectCores() - 1)
    clusterExport(
        cl,
        varlist = c(
            "patient_age", "patient_vaf",
            "grid", "w0", "R", "alpha",
            "loglikelihood_timeseries", "Jg", "A_analytic", "compute_expected_means_and_variances"
        ),
        envir = environment()
    )
    loglikelihoods <- parLapply(
        cl = cl, 1:nrow(grid),
        function(i) {
            p <- c(w0, grid$w1[i], grid$v0[i], alpha, R)
            loglikelihood <- loglikelihood_timeseries(p, patient_age, patient_vaf)
        }
    )
    stopCluster(cl)
    grid$loglikelihood <- unlist(loglikelihoods)
    # #---Sequential version
    # grid$loglikelihood <- NA_real_
    # for (i in seq_len(nrow(grid))) {
    #     p <- c(w0, grid$w1[i], grid$v0[i], alpha, R)
    #     print(p)
    #     loglikelihood <- tryCatch(loglikelihood_timeseries(p, patient_age, patient_vaf),
    #         error = function(e) {
    #             message("loglikelihood_timeseries error at row ", i, ": ", e$message)
    #             -Inf
    #         }
    #     )
    #     grid$loglikelihood[i] <- as.numeric(loglikelihood)
    #     if (i %% 10000 == 0) cat(sprintf("  progress: %d / %d\n", i, nrow(grid)))
    # }
})
print(timing)
#---Compute posterior distribution
grid$posterior <- exp(grid$loglikelihood - max(grid$loglikelihood))
#---Save posterior distribution
saveRDS(grid, file = file.path(outdir, paste0("posterior_", mutation_ID, "_", patient_ID, ".rds")))
write.csv(grid,
    file = file.path(outdir, paste0("posterior_", mutation_ID, "_", patient_ID, ".csv")),
    row.names = FALSE
)
#---Find MAP in posterior distribution
#   Get MAP parameters from posterior distribution
map_idx <- which.max(grid$posterior)
map_row <- grid[map_idx, ]
#   Simulate trajectories with MAP parameters
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
pdf(file.path(outdir, paste0("Marginal_distributions_", mutation_ID, "_", patient_ID, ".pdf")), width = 20, height = 10)
p <- gridExtra::arrangeGrob(p1, p2,
    ncol = 2,
    top = paste0("Marginal posterior distributions (", mutation_ID, "; ", patient_ID, ")")
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
pdf(file.path(outdir, paste0("Joint_distributions_", mutation_ID, "_", patient_ID, ".pdf")), width = 10, height = 10)
p <- gridExtra::arrangeGrob(p,
    ncol = 1,
    top = paste0("Joint posterior distributions (", mutation_ID, "; ", patient_ID, ")")
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
    geom_point(data = fabre_short, aes(x = Age, y = VAF), color = "white", size = 14) +
    geom_point(data = fabre_short, aes(x = Age, y = VAF), color = "#00A087", size = 10) +
    labs(
        x = "Age",
        y = "Variant Allele Frequency (VAF)",
        title = NULL
    ) +
    theme_minimal(base_size = 50) +
    FORCE_WHITE +
    theme(legend.position = "none")
pdf(file.path(outdir, paste0("VAF_trajectories_", mutation_ID, "_", patient_ID, ".pdf")), width = 30, height = 15)
print(p)
dev.off()
