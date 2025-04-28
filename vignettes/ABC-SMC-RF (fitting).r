setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/abcsmcrf/R")
# setwd("/Users/keitotaketomi/Downloads/abc-smc-rf 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R")
# setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")
# setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep 2/vignettes")



###############################################################################
# 2) LOAD THE REAL PATIENT TIMING TABLE
###############################################################################
real_df <- read.delim(
    "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/evolution_and_heterogeneity/2018-07-24-wgdMrcaTiming.txt",
    #   "/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.csv",
    stringsAsFactors = FALSE, header = TRUE
)

CANCER_TISSUE <- "Liver-HCC"
WGD_STATUS <- FALSE

# pick the right MRCA column + rename 'age' → 'diag_age'
if (!WGD_STATUS) {
    real_sub <- real_df %>%
        filter(tissue == CANCER_TISSUE, WGD == FALSE) %>%
        rename(
            MRCA_age = MRCA.time.linear,
            diag_age = age
        ) %>%
        filter(!is.na(MRCA_age), !is.na(diag_age))
    real_sub$t <- real_sub$MRCA_age
} else {
    real_sub <- real_df %>%
        filter(tissue == CANCER_TISSUE, WGD == TRUE) %>%
        rename(
            MRCA_age = MRCA.time.linear,
            diag_age = age
        )
}

histogram_x <- 0:(ceiling(max(real_sub$t)))
histogram_y <- hist(real_sub$t, breaks = histogram_x, plot = FALSE)$counts
histogram_y <- c(histogram_y / sum(histogram_y), 0)

statistics_target <- data.frame(matrix(histogram_y, nrow = 1))
colnames(statistics_target) <- c(paste0("Age_group_", histogram_x[2:length(histogram_x)]), "Age_group_NA")

model <- function(parameters, parallel = TRUE) {
    one_parameter <- function(parameter) {
        lambda <- as.numeric(parameter$lambda)
        #---Model parameters
        r_initial <- c(10000, 1)
        lambda_vec <- c(1, lambda)
        u_vec <- c(0, 0)
        alpha <- 1
        max_time <- histogram_x[length(histogram_x)]
        tau <- 0.01
        threshold_diagnosis <- 0.2
        n_simulations <- 100
        #---Simulate MRCA ages & diagnosis ages
        diagnosis_ages <- rep(NA, n_simulations)
        for (i in 1:n_simulations) {
            output <- simulate_continuous_moran_tau(
                r = r_initial,
                lambda_vec = lambda_vec,
                u_vec = u_vec,
                alpha = alpha,
                max_time = max_time,
                tau = tau
            )
            output$N1_frequency <- output$N1 / (output$N0 + output$N1)
            diagnosis_age <- output$time[which(output$N1_frequency > threshold_diagnosis)[1]]
            diagnosis_ages[i] <- diagnosis_age
        }
        simulation_histogram_y <- hist(diagnosis_ages, breaks = histogram_x, plot = FALSE)$counts
        simulation_histogram_y <- c(simulation_histogram_y, n_simulations - sum(simulation_histogram_y))
        simulation_histogram_y <- simulation_histogram_y / sum(simulation_histogram_y)
        statistics_output <- data.frame(matrix(simulation_histogram_y, nrow = 1))
        colnames(statistics_output) <- c(paste0("Age_group_", histogram_x[2:length(histogram_x)]), "Age_group_NA")
        return(statistics_output)
    }
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(
            cl,
            varlist = c("one_parameter", "simulate_continuous_moran_tau", "histogram_x"),
            envir = environment()
        )
        stats <- parLapply(
            cl = cl, 1:nrow(parameters),
            function(i) {
                sub_parameters <- as.data.frame(parameters[i, ])
                colnames(sub_parameters) <- colnames(parameters)
                one_parameter(sub_parameters)
                # one_parameter(parameters[i, ])
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            sub_parameters <- as.data.frame(parameters[i, ])
            colnames(sub_parameters) <- colnames(parameters)
            stats <- rbind(
                stats,
                one_parameter(sub_parameters)
                # one_parameter(parameters[i, ])
            )
        }
    }
    #   Add column names
    data <- cbind(parameters, as.data.frame(stats))
    return(data)
}

parameters <- data.frame(lambda = 1:10)
print(model(parameters))








###############################################################################
# 3) FUNCTION TO INVERT λ FROM Mutant‐fraction = 20%
###############################################################################
fit_lambda_by_threshold <- function(
    t,
    N0_init = 9999,
    N1_init = 1,
    lambda0 = 1, # wild‐type growth rate
    detect = 0.20, # 20% detection threshold
    grid_min = -5,
    grid_max = +5,
    n_grid = 101) {
    if (t <= 0) {
        return(NA_real_)
    }

    # the equation whose root we seek:
    f_detect <- function(lambda1) {
        N0t <- N0_init * exp(lambda0 * t)
        N1t <- N1_init * exp(lambda1 * t)
        (N1t / (N0t + N1t)) - detect
    }

    # 1) build a grid of lambdas
    lambdas <- seq(grid_min, grid_max, length.out = n_grid)
    vals <- vapply(lambdas, f_detect, numeric(1))

    # 2) drop any NA points
    ok <- is.finite(vals)
    lambdas <- lambdas[ok]
    vals <- vals[ok]

    # 3) look for the first place where sign flips
    sgn <- sign(vals)
    change_pts <- which(diff(sgn) != 0)

    if (length(change_pts) == 0) {
        # no crossing found
        return(NA_real_)
    }

    # bracket for uniroot:
    i <- change_pts[1]
    lower <- lambdas[i]
    upper <- lambdas[i + 1]

    # 4) solve:
    sol <- tryCatch(
        uniroot(f_detect, lower = lower, upper = upper)$root,
        error = function(e) NA_real_
    )

    sol
}


# apply to each patient
real_sub$lambda_est <- vapply(real_sub$t, fit_lambda_by_threshold, numeric(1))

# summarise per‐cancer
lambda_cancer <- median(real_sub$lambda_est, na.rm = TRUE)
message(
    "Fitted λ for ", CANCER_TISSUE,
    "  WGD=", WGD_STATUS, "  → ", round(lambda_cancer, 3)
)

# override your 'true' λ
parameters_truth$lambda <- lambda_cancer
###############################################################################
# 2) FUNCTIONS FOR ABC-SMC-RF
###############################################################################
simulate_one_run <- function(parameter_set) {
    # Convert from log10(u) to actual u
    u <- 10^as.numeric(parameter_set$log_u)
    lambda <- as.numeric(parameter_set$lambda)
    if (is.na(u) || is.na(lambda) || u <= 0 || lambda <= 0) {
        output <- data.frame(matrix(NA, nrow = 1, ncol = 2 * n_patients))
        output_colnames <- c()
        for (i in 1:n_patients) {
            output_colnames <- c(output_colnames, paste0("MRCA_age_", i))
            output_colnames <- c(output_colnames, paste0("Diagnosis_age_", i))
        }
        colnames(output) <- output_colnames
        return(output)
    }
    #---Model parameters
    r_initial <- 10000
    lambda_vec <- c(0.2, lambda)
    u_vec <- c(u, 0)
    alpha <- 1
    max_time <- 100
    tau <- 0.01

    threshold_MRCA <- 1 / r_initial
    threshold_diagnosis <- 0.2

    #---Simulate MRCA ages & diagnosis ages
    MRCA_ages <- c()
    diagnosis_ages <- c()
    for (i in 1:n_patients) {
        output <- simulate_continuous_moran_tau(
            r = r_initial,
            lambda_vec = lambda_vec,
            u_vec = u_vec,
            alpha = alpha,
            max_time = max_time,
            tau = tau
        )
        output$N1_frequency <- output$N1 / (output$N0 + output$N1)
        MRCA_age <- output$time[which(output$N1_frequency > threshold_MRCA)[1]]
        diagnosis_age <- output$time[which(output$N1_frequency > threshold_diagnosis)[1]]
        MRCA_ages <- c(MRCA_ages, MRCA_age)
        diagnosis_ages <- c(diagnosis_ages, diagnosis_age)
    }
    sorted_indices <- order(diagnosis_ages)
    diagnosis_ages <- diagnosis_ages[sorted_indices]
    MRCA_ages <- MRCA_ages[sorted_indices]

    output <- data.frame(matrix(NA, nrow = 1, ncol = 0))
    for (i in 1:n_patients) {
        output[[paste0("MRCA_age_", i)]] <- MRCA_ages[i]
        output[[paste0("Diagnosis_age_", i)]] <- diagnosis_ages[i]
    }
    output[which(is.na(output))] <- 2 * max_time
    return(output)
}

model <- function(parameters, parallel = TRUE) {
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(
            cl,
            varlist = c("simulate_one_run", "simulate_continuous_moran_tau", "n_patients"),
            envir = environment()
        )
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) simulate_one_run(parameters[i, ])
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            print(i)
            stats_one_simulation <- simulate_one_run(parameters[i, ])
            stats <- rbind(stats, stats_one_simulation)
        }
    }
    output <- cbind(parameters, stats)
    return(output)
}

###############################################################################
# 3) FUNCTIONS FOR PLOTTING
###############################################################################
# Note: annotation_text now automatically uses the hyperparameter variables.
annotation_text <- paste(
    "Truth: lambda =", parameters_truth$lambda,
    ", u =", format(10^(parameters_truth$log_u), scientific = TRUE),
    ", NUM_PARTICLES =", NUM_PARTICLES,
    ", NUM_ITERATIONS =", NUM_ITERATIONS
)

plot_ages <- function(smcdrf_results) {
    library(dplyr)
    library(ggplot2)
    statistics_target <- smcdrf_results$statistics_target
    final_iter_stats <- smcdrf_results[[paste0("Iteration_", smcdrf_results$nIterations + 1)]]$statistics
    # Reshape the single-row target data into a data frame of (patient, MRCA_age, Diagnosis_age)
    target_long <- data.frame(
        patient = seq_len(n_patients),
        MRCA_age = as.numeric(statistics_target[1, grep("^MRCA_age_", names(statistics_target))]),
        Diagnosis_age = as.numeric(statistics_target[1, grep("^Diagnosis_age_", names(statistics_target))])
    )
    target_long$set <- "Target"
    # For each row in final_iter_stats, reshape it to the same long format
    all_final_list <- list()
    for (i in seq_len(nrow(final_iter_stats))) {
        mrca_cols <- grep("^MRCA_age_", names(final_iter_stats), value = TRUE)
        diag_cols <- grep("^Diagnosis_age_", names(final_iter_stats), value = TRUE)
        df_i <- data.frame(
            patient = seq_len(n_patients),
            MRCA_age = as.numeric(final_iter_stats[i, mrca_cols]),
            Diagnosis_age = as.numeric(final_iter_stats[i, diag_cols])
        )
        df_i$set <- "Final iteration"
        all_final_list[[i]] <- df_i
    }
    final_long <- do.call(rbind, all_final_list)
    p <- ggplot() +
        geom_density_2d_filled(data = final_long, aes(x = Diagnosis_age, y = MRCA_age)) +
        geom_point(data = target_long, aes(x = Diagnosis_age, y = MRCA_age), color = "#A50021", size = 3) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E0FFFF", size = 2) +
        theme_minimal(base_size = 14) +
        scale_fill_grey() +
        xlim(c(0, 80)) +
        ylim(c(0, 80)) +
        labs(
            title = "Simulated Diagnosis vs. MRCA Ages",
            subtitle = annotation_text,
            x = "Diagnosis age",
            y = "MRCA age"
        ) +
        theme(legend.position = "none")
    file_name <- paste0("Posterior_ages_", smcdrf_results$method, ".png")
    png(file_name, res = 150, width = 10, height = 10, units = "in", pointsize = 12)
    print(p)
    dev.off()
}

###############################################################################
# 4) TEST WITH SYNTHETIC DATA
###############################################################################

statistics_target <- simulate_one_run(parameters_truth)

###############################################################################
# 5) HYPERPARAMETERS FOR ABC-SMC-RF
###############################################################################
parameters_labels <- data.frame(
    parameter = c("log_u", "lambda"),
    label = c("expression(log[10](u))", "expression(lambda)")
)
log_u_range <- c(-5, -2)
lambda_range <- c(0, 1)
rprior <- function(Nparameters) {
    data.frame(
        log_u = runif(Nparameters, log_u_range[1], log_u_range[2]),
        lambda = runif(Nparameters, lambda_range[1], lambda_range[2])
    )
}
dprior <- function(parameters) {
    probs <- dunif(parameters$log_u, log_u_range[1], log_u_range[2]) *
        dunif(parameters$lambda, lambda_range[1], lambda_range[2])
    return(probs)
}

###############################################################################
# 6) ABC-SMC-RF FOR MULTIPLE PARAMETERS
###############################################################################
smcdrf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(NUM_PARTICLES, NUM_ITERATIONS),
    num.trees = NUM_TREES,
    parallel = TRUE
)

#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)

#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    abc_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)

#---Plot simulated diagnosis & MRCA ages against statistics target
plot_ages(smcdrf_results = smcdrf_results)

###############################################################################
# 7) ABC-SMC-RF FOR SINGLE PARAMETERS
###############################################################################
smcdrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(NUM_PARTICLES, NUM_ITERATIONS),
    ntrees = NUM_TREES,
    parallel = TRUE
)

#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)

#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = smcdrf_results,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)

#---Plot simulated diagnosis & MRCA ages against statistics target
plot_ages(smcdrf_results = smcdrf_results)

###############################################################################
# 8) OPTIONAL: Save additional PDF files with the log(u) plots and annotations
###############################################################################
pdf(file.path(output_dir, "Figure1_Multi_Marginal_LogU.pdf"), width = 10, height = 8)
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
dev.off()

pdf(file.path(output_dir, "Figure2_Single_Marginal_LogU.pdf"), width = 10, height = 8)
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
dev.off()

pdf(file.path(output_dir, "Figure3_Multi_Compare_LogU.pdf"), width = 10, height = 8)
plot_compare_marginal(
    abc_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)
dev.off()

pdf(file.path(output_dir, "Figure4_Single_Compare_LogU.pdf"), width = 10, height = 8)
plot_compare_marginal(
    plots = plots,
    abc_results = smcdrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)
dev.off()
