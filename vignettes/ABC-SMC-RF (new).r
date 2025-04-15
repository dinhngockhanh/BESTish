library(invgamma)
library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(1)

setwd("/Users/keitotaketomi/Downloads/abc-smc-rf 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/keitotaketomi/Downloads/DriverSelectionSweep 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/keitotaketomi/Downloads/DriverSelectionSweep 2/vignettes")

CANCER_TISSUE <- "Liver-HCC"
WGD_STATUS <- FALSE
output_dir <- paste0(gsub("-", "_", tolower(CANCER_TISSUE)), "_wgd_", tolower(WGD_STATUS))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# FUNCTIONS FOR ABC-SMC-RF
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
    # ---Model parameters
    r_initial <- 10000
    lambda_vec <- c(0.2, lambda)
    u_vec <- c(u, 0)
    alpha <- 1
    max_time <- 100
    tau <- 0.01

    threshold_MRCA <- 1 / r_initial
    threshold_diagnosis <- 0.2

    # ---Simulate MRCA ages & diagnosis ages
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

# LOAD ACTUAL DATA AND BUILD OBSERVED STATISTICS 
pcawg_df <- read.csv("/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.csv", 
                     header = TRUE, stringsAsFactors = FALSE)
pcawg_df <- subset(pcawg_df, tissue == "Liver-HCC" & WGD == FALSE)

pcawg_df<-pcawg_df[which(!is.na(pcawg_df$age) & !is.na(pcawg_df$MRCA.time.linear)),]

set.seed(1)  
pcawg_df <- pcawg_df[sample(nrow(pcawg_df), 30), ]
n_patients <- nrow(pcawg_df)

order_index <- order(pcawg_df$age)
diagnosis_ages <- pcawg_df$age[order_index]
MRCA_ages <- diagnosis_ages - pcawg_df$MRCA.time.linear[order_index]

# Build the observed statistics_target as a single-row data frame
statistics_target <- data.frame(matrix(NA, nrow = 1, ncol = 0))
for (i in 1:n_patients) {
    statistics_target[[paste0("MRCA_age_", i)]] <- MRCA_ages[i]
    statistics_target[[paste0("Diagnosis_age_", i)]] <- diagnosis_ages[i]
}

parameters_labels <- data.frame(
    parameter = c("log_u", "lambda"),
    label = c("expression(log[10](u))", "expression(lambda)")
)
log_u_range <- c(-8, -2)
lambda_range <- c(0, 5)
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
NUM_TREES <- 500
NUM_PARTICLES <- 5000
NUM_ITERATIONS <- 3

# FUNCTIONS FOR PLOTTING
plot_ages <- function(smcdrf_results) {
    library(dplyr)
    library(ggplot2)
    statistics_target <- smcdrf_results$statistics_target
    final_iter_stats <- smcdrf_results[[paste0("Iteration_", smcdrf_results$nIterations + 1)]]$statistics
    target_long <- data.frame(
        patient = seq_len(n_patients),
        MRCA_age = as.numeric(statistics_target[1, grep("^MRCA_age_", names(statistics_target))]),
        Diagnosis_age = as.numeric(statistics_target[1, grep("^Diagnosis_age_", names(statistics_target))])
    )
    target_long$set <- "Target"
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
            x = "Diagnosis age",
            y = "MRCA age"
        ) +
        theme(
            legend.position = "none"
        )
    file_name <- paste0("Posterior_ages_", smcdrf_results$method, ".png")
    png(file_name, res = 150, width = 10, height = 10, units = "in", pointsize = 12)
    print(p)
    dev.off()
}

# ABC-SMC-RF FOR MULTIPLE PARAMETERS
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
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)

#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    abc_results = smcdrf_results,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)

#---Plot simulated diagnosis & MRCA ages against statistics target
plot_ages(smcdrf_results = smcdrf_results)

# ABC-SMC-RF FOR SINGLE PARAMETERS
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
    plot_hist = TRUE
)

#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = smcdrf_results,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)

#---Plot simulated diagnosis & MRCA ages against statistics target
plot_ages(smcdrf_results = smcdrf_results)

# Save additional PDF files with the log(u) plots 
pdf(file.path(output_dir, "Figure1_Multi_Marginal_LogU.pdf"), width = 10, height = 8)
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
dev.off()

pdf(file.path(output_dir, "Figure2_Single_Marginal_LogU.pdf"), width = 10, height = 8)
plot_smcrf_marginal(
    smcrf_results = smcdrf_results,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
dev.off()

pdf(file.path(output_dir, "Figure3_Multi_Compare_LogU.pdf"), width = 10, height = 8)
plot_compare_marginal(
    abc_results = smcdrf_results,
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
    parameters_labels = parameters_labels,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
)
dev.off()
