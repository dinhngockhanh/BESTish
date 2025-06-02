setwd("/Users/keitotaketomi/Downloads/abc-smc-rf 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/vignettes")

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(parallel)
library(pbapply)
# ============================================Input MRCA-WGD timing data
real_df <- read.table(
    "/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.txt",
    header = TRUE, stringsAsFactors = FALSE
)
# =======================================================Objective model
model <- function(parameters, statistic = "distance", parallel = TRUE) {
    one_parameter <- function(parameter) {
        #---Extract parameters to be inferred
        log10_r0 <- as.numeric(parameter$log10_r0) # log10(r₀)
        lambda1 <- as.numeric(parameter$lambda1) # λ₁
        lambda0 <- as.numeric(parameter$lambda0) # λ₀
        log10_thresh <- as.numeric(parameter$log10_threshold) # log10(threshold)
        alpha <- as.numeric(parameter$alpha) # α
        r0 <- 10^log10_r0
        threshold <- 10^log10_thresh
        #---Confirm fixed parameters
        u_vec <- c(0, 0)
        tau <- 0.01
        n_sim <- 100
        ################################################################
        ################################################################
        ################################################################
        # (b) Run a “pilot” simulation to discover the simulation’s time grid length
        pilot_out <- simulate_continuous_moran_tau(
            r = c(r0, lambda0),
            lambda_vec = c(lambda0, lambda1),
            u_vec = u_vec,
            alpha = alpha,
            max_time = histogram_x[length(histogram_x)],
            tau = tau
        )
        time_grid <- pilot_out$time
        n_time <- length(time_grid)
        # Pre‐allocate a matrix to store N₁(t) for each of n_sim replicates
        N1_mat <- matrix(NA_real_, nrow = n_sim, ncol = n_time)
        N1_mat[1, ] <- pilot_out$N1
        # Run the remaining (n_sim − 1) replicates
        for (i in 2:n_sim) {
            sim_out <- simulate_continuous_moran_tau(
                r = c(r0, lambda0),
                lambda_vec = c(lambda0, lambda1),
                u_vec = u_vec,
                alpha = alpha,
                max_time = histogram_x[length(histogram_x)],
                tau = tau
            )
            # Guard against variable‐length outputs:
            if (length(sim_out$N1) == n_time) {
                N1_mat[i, ] <- sim_out$N1
            } else if (length(sim_out$N1) > n_time) {
                N1_mat[i, ] <- sim_out$N1[seq_len(n_time)]
            } else {
                # If sim_out$N1 is shorter, copy what we have, pad the rest with last known N₁ (or 0)
                tmp_len <- length(sim_out$N1)
                if (tmp_len > 0) {
                    N1_mat[i, 1:tmp_len] <- sim_out$N1
                    pad_val <- if (!is.na(sim_out$N1[tmp_len])) sim_out$N1[tmp_len] else 0
                } else {
                    pad_val <- 0
                }
                if (tmp_len < n_time) {
                    N1_mat[i, (tmp_len + 1):n_time] <- pad_val
                }
            }
        }
        # (c) Compute “old diagnosis time” from the MEAN curve N₁̄(t)
        mean_N1 <- colMeans(N1_mat, na.rm = TRUE)
        idx_old <- which(mean_N1 > threshold)
        old_time <- if (length(idx_old) > 0) time_grid[idx_old[1]] else NA_real_
        # (d) Compute each replicate’s “new diagnosis time” (first t where N₁ > threshold)
        new_times <- rep(NA_real_, n_sim)
        for (i in seq_len(n_sim)) {
            idx_new <- which(N1_mat[i, ] > threshold)
            new_times[i] <- if (length(idx_new) > 0) time_grid[idx_new[1]] else NA_real_
        }
        new_times <- new_times[!is.na(new_times)]
        # (e) Build the CDF of “new_times” over the 10 uniform bins (histogram_x)
        if (length(new_times) > 0) {
            cnts <- hist(new_times, breaks = histogram_x, plot = FALSE)$counts
            freq <- cnts / sum(cnts)
            sim_cdf <- cumsum(freq)
        } else {
            # If no replicate ever crosses threshold, treat all bins as zero
            sim_cdf <- rep(0, length(histogram_x) - 1)
        }
        # (f) Return distance or the raw simulated‐CDF
        if (statistic == "distance") {
            d <- sum(abs(sim_cdf - cdf_obs))
            if (is.na(d) || !is.finite(d)) d <- Inf
            return(data.frame(distance = d))
        } else {
            df <- as.data.frame(t(sim_cdf))
            colnames(df) <- paste0("Age_bin_", seq_along(sim_cdf))
            return(df)
        }
        ################################################################
        ################################################################
        ################################################################
    }
    #---Apply one_parameter() to every row of “parameters”
    if (parallel) {
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl,
            c(
                "one_parameter",
                "simulate_continuous_moran_tau",
                "histogram_x", "cdf_obs", "statistic"
            ),
            envir = environment()
        )
        stats_list <- parLapply(
            cl, seq_len(nrow(parameters)), function(i) {
                one_parameter(as.data.frame(parameters[i, , drop = FALSE]))
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats_list)
    } else {
        stats <- rbindlist(
            lapply(seq_len(nrow(parameters)), function(i) {
                one_parameter(as.data.frame(parameters[i, , drop = FALSE]))
            })
        )
    }
    #---Output df containing parameters and corresponding statistics
    cbind(parameters, stats)
}
# ====================================================Prior distribution
rprior <- function(Nparameters) {
    #---log10(r₀) ∼ Uniform(2,6)
    log10_r0 <- runif(Nparameters, min = 2, max = 6)
    #---λ₁ ∼ Uniform(0,20)
    lambda1 <- runif(Nparameters, min = 0, max = 20)
    #---λ₀ ∼ Uniform(0, λ₁)
    lambda0 <- runif(Nparameters, min = 0, max = lambda1)
    #---log10(threshold) ∼ Uniform(1, log10(r₀))
    log10_threshold <- runif(Nparameters, min = 1, max = log10_r0)
    #---α ∼ Uniform(0.5,1)
    alpha <- runif(Nparameters, min = 0.5, max = 1)
    #---Return sampled parameters
    data.frame(
        log10_r0        = log10_r0,
        lambda1         = lambda1,
        lambda0         = lambda0,
        log10_threshold = log10_threshold,
        alpha           = alpha
    )
}
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "log10_r0")) {
        probs <- probs * dunif(parameters[["log10_r0"]], min = 2, max = 6)
    }
    if (parameter_id %in% c("all", "lambda1")) {
        probs <- probs * dunif(parameters[["lambda1"]], min = 0, max = 20)
    }
    if (parameter_id %in% c("all", "lambda0")) {
        probs <- probs * dunif(parameters[["lambda0"]], min = 0, max = parameters[["lambda1"]])
    }
    if (parameter_id %in% c("all", "log10_threshold")) {
        probs <- probs * dunif(parameters[["log10_threshold"]], min = 1, max = parameters[["log10_r0"]])
    }
    if (parameter_id %in% c("all", "alpha")) {
        probs <- probs * dunif(parameters[["alpha"]], min = 0.5, max = 1)
    }
    return(probs)
}

NUM_PARTICLES <- 5000
NUM_ITERATIONS <- 5
NUM_TREES <- 500

# ─────────────────────────────────────────────────────────────────────────────
# 3) LOOP OVER EACH CANCER TISSUE & WGD_STATUS
# ─────────────────────────────────────────────────────────────────────────────
CANCER_TISSUES <- unique(real_df$tissue)

for (CANCER_TISSUE in CANCER_TISSUES) {
    for (WGD_STATUS in c(FALSE, TRUE)) {
        message("——————————————————————————————")
        message(" Processing tissue: ", CANCER_TISSUE, "  |  WGD = ", WGD_STATUS)
        message("——————————————————————————————")

        # 3a) Pick & preprocess the subset of real data
        real_sub <- real_df %>%
            filter(tissue == CANCER_TISSUE, WGD == WGD_STATUS) %>%
            rename(
                MRCA_age = MRCA.time.linear,
                diag_age = age
            ) %>%
            filter(!is.na(MRCA_age))

        # Skip if too few samples:
        if (nrow(real_sub) < 10) {
            message("  → Too few samples (", nrow(real_sub), "), skipping.")
            next
        }
        real_sub$t <- real_sub$MRCA_age

        # (1) Build 10 uniform bins from t=0 to t=max_t:
        max_t <- ceiling(max(real_sub$t, na.rm = TRUE))
        if (max_t < 1) max_t <- 1
        histogram_x <- seq(0, max_t, length.out = 11) # 11 breakpoints → 10 bins
        histogram_x <- unique(histogram_x)
        if (length(histogram_x) < 2) histogram_x <- c(0, 1)

        # (2) Compute observed CDF over those 10 bins
        counts_obs <- hist(real_sub$t, breaks = histogram_x, plot = FALSE)$counts
        freq_obs <- counts_obs / sum(counts_obs)
        cdf_obs <- cumsum(freq_obs)

        # Sanity‐check: does “distance” vary across random draws?
        test_params <- rprior(10)
        test_ds <- model(test_params, statistic = "distance", parallel = FALSE)$distance
        if (length(unique(test_ds)) < 2) {
            message("  → Summary is constant across λ, skipping ABC.")
            next
        }

        # 3b) Run ABC‐SMC‐RF (wrapped in tryCatch so we skip on errors)
        smcrf_results <- tryCatch(
            {
                smcrf(
                    method            = "smcrf-single-param",
                    statistics_target = data.frame(distance = 0),
                    model             = model,
                    rprior            = rprior,
                    dprior            = dprior,
                    nParticles        = rep(NUM_PARTICLES, NUM_ITERATIONS),
                    ntrees            = NUM_TREES,
                    parallel          = TRUE
                )
            },
            error = function(e) {
                message("  → smcrf failed: ", conditionMessage(e))
                return(NULL)
            }
        )
        if (is.null(smcrf_results)) next

        # 4) Plot posterior marginal over λ₁
        plot_compare_marginal(
            abc_results     = smcrf_results,
            plot_statistics = FALSE,
            plot_hist       = TRUE,
            plot_prior      = TRUE
        )
        file.rename(
            "comparison-marginal-parameter=lambda1.png",
            paste0(CANCER_TISSUE, "_WGD=", WGD_STATUS, "_lambda1.png")
        )

        # 5) Simulated vs Observed CDF
        final_ps <- smcrf_results[[paste0("Iteration_", NUM_ITERATIONS + 1)]]$parameters
        sims_df <- model(
            final_ps[1:100, , drop = FALSE],
            statistic = "histogram",
            parallel = FALSE
        )

        sim_long <- sims_df %>%
            mutate(simID = row_number()) %>%
            pivot_longer(
                cols      = starts_with("Age_bin_"),
                names_to  = "bin",
                values_to = "cdf"
            ) %>%
            mutate(
                bin = as.integer(sub("Age_bin_", "", bin)),
                mid = (histogram_x[bin] + histogram_x[bin + 1]) / 2
            )

        target_df <- tibble(
            bin = seq_along(cdf_obs),
            mid = (histogram_x[-length(histogram_x)] + histogram_x[-1]) / 2,
            cdf = cdf_obs
        )

        p <- ggplot() +
            geom_line(
                data = sim_long,
                aes(x = mid, y = cdf, group = simID),
                color = "steelblue",
                alpha = 0.2
            ) +
            geom_line(
                data  = target_df,
                aes(x = mid, y = cdf),
                color = "firebrick",
                size  = 1.2
            ) +
            theme_minimal(base_size = 14) +
            labs(
                title = paste0(CANCER_TISSUE, "  (WGD=", WGD_STATUS, ")"),
                x     = "MRCA age",
                y     = "Cumulative frequency"
            )

        ggsave(
            filename = paste0(CANCER_TISSUE, "_WGD=", WGD_STATUS, "_cdf.png"),
            plot     = p,
            dpi      = 150,
            width    = 12,
            height   = 6,
            units    = "in"
        )
    }
}
