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

# ─────────────────────────────────────────────────────────────────────────────
# 1) LOAD THE REAL PATIENT TIMING TABLE
# ─────────────────────────────────────────────────────────────────────────────
real_df <- read.table(
    "/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.txt",
    header = TRUE, stringsAsFactors = FALSE
)

# ─────────────────────────────────────────────────────────────────────────────
# 2) DEFINE MODEL, PRIORS & ABC‐SMC‐RF WRAPPER
# ─────────────────────────────────────────────────────────────────────────────
model <- function(parameters, statistic = "distance", parallel = TRUE) {
    one_parameter <- function(parameter) {
        #---Parameters for fitting
        λ <- as.numeric(parameter$lambda)
        #---Parameters that are currently fixed >>> add these to the fitting
        r_initial <- c(1e8, 1) # r_0
        lambda_vec <- c(1e-5, λ) # lambda_0
        # threshold <- 0.001
        threshold_new <- 1000 # threshold for diagnosis
        alpha <- 1 # expansion rate
        #---Parameters that are fixed
        u_vec <- c(0, 0)
        tau <- 0.01
        n_sim <- 100

        # simulate n_sim ages
        ages <- replicate(n_sim, {
            out <- simulate_continuous_moran_tau(
                r = r_initial,
                lambda_vec = lambda_vec,
                u_vec = u_vec,
                alpha = alpha,
                max_time = histogram_x[length(histogram_x)],
                tau = tau
            )
            # f <- out$N1 / (out$N0 + out$N1)
            # out$time[which(f > threshold)[1]]
            out$time[which(out$N1 > threshold_new)[1]]
        })
        ages <- ages[!is.na(ages)]

        # CDF over the 10 uniform bins
        cnts <- hist(ages, breaks = histogram_x, plot = FALSE)$counts
        freq <- cnts / sum(cnts)
        sim_cdf <- cumsum(freq)

        if (statistic == "distance") {
            d <- sum(abs(sim_cdf - cdf_obs))
            if (is.na(d)) d <- Inf
            return(data.frame(distance = d))
        } else {
            df <- as.data.frame(t(sim_cdf))
            colnames(df) <- paste0("Age_bin_", seq_along(sim_cdf))
            return(df)
        }
    }

    if (parallel) {
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl,
            c(
                "one_parameter", "simulate_continuous_moran_tau",
                "histogram_x", "cdf_obs", "statistic"
            ),
            envir = environment()
        )
        out <- parLapply(cl, seq_len(nrow(parameters)), function(i) {
            one_parameter(as.data.frame(parameters[i, , drop = FALSE]))
        })
        stopCluster(cl)
        stats <- rbindlist(out)
    } else {
        stats <- rbindlist(lapply(seq_len(nrow(parameters)), function(i) {
            one_parameter(as.data.frame(parameters[i, , drop = FALSE]))
        }))
    }

    cbind(parameters, stats)
}

lambda_min <- 0
lambda_max <- 20
rprior <- function(Nparameters) {
    data.frame(lambda = runif(Nparameters, lambda_min, lambda_max))
}
dprior <- function(parameters, parameter_id = "all") {
    dunif(parameters$lambda, lambda_min, lambda_max)
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
        message(
            " Processing tissue: ", CANCER_TISSUE,
            " | WGD = ", WGD_STATUS
        )
        message("——————————————————————————————")

        # 3a) pick & preprocess
        real_sub <- real_df %>%
            filter(tissue == CANCER_TISSUE, WGD == WGD_STATUS) %>%
            rename(MRCA_age = MRCA.time.linear, diag_age = age) %>%
            filter(!is.na(MRCA_age))

        # **skip tiny datasets**
        if (nrow(real_sub) < 10) {
            message(" → too few samples (", nrow(real_sub), "), skipping.")
            next
        }
        real_sub$t <- real_sub$MRCA_age

        # **create 10 uniform bins** (force at least two distinct breaks)
        max_t <- ceiling(max(real_sub$t, na.rm = TRUE))
        if (max_t < 1) max_t <- 1
        histogram_x <- seq(0, max_t, length.out = 11)
        histogram_x <- unique(histogram_x)
        if (length(histogram_x) < 2) histogram_x <- c(0, 1)

        # observed CDF
        counts <- hist(real_sub$t, breaks = histogram_x, plot = FALSE)$counts
        freq <- counts / sum(counts)
        cdf_obs <- cumsum(freq)

        # sanity check for distance variation
        test_ds <- model(rprior(10), "distance", FALSE)$distance
        if (length(unique(test_ds)) < 2) {
            message(" → summary constant across λ, skipping ABC.")
            next
        }

        # 3b) run ABC‐SMC‐RF (**wrapped in tryCatch**)
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
                message(" → smcrf failed: ", e$message)
                NULL
            }
        )
        if (is.null(smcrf_results)) next

        # 4) plot posterior marginal
        plot_compare_marginal(
            abc_results     = smcrf_results,
            plot_statistics = FALSE,
            plot_hist       = TRUE,
            plot_prior      = TRUE
        )
        file.rename(
            "comparison-marginal-parameter=lambda.png",
            paste0(
                CANCER_TISSUE,
                "_WGD=", WGD_STATUS, "_lambda.png"
            )
        )

        # 5) sim vs target CDF
        finals <- smcrf_results[[paste0(
            "Iteration_",
            NUM_ITERATIONS + 1
        )]]$parameters
        sims_df <- model(
            finals[1:100, , drop = FALSE],
            statistic = "histogram",
            parallel = FALSE
        )

        sim_long <- sims_df %>%
            mutate(simID = row_number()) %>%
            pivot_longer(starts_with("Age_bin_"),
                names_to = "bin", values_to = "cdf"
            ) %>%
            mutate(
                bin = as.integer(sub("Age_bin_", "", bin)),
                mid = (histogram_x[bin] + histogram_x[bin + 1]) / 2
            )

        target_df <- tibble(
            bin = seq_along(cdf_obs),
            mid = (histogram_x[-length(histogram_x)] +
                histogram_x[-1]) / 2,
            cdf = cdf_obs
        )

        p <- ggplot() +
            geom_line(
                data = sim_long,
                aes(mid, cdf, group = simID),
                color = "steelblue", alpha = 0.2
            ) +
            geom_line(
                data = target_df,
                aes(mid, cdf),
                color = "firebrick", size = 1.2
            ) +
            theme_minimal(base_size = 14) +
            labs(
                title = paste0(
                    CANCER_TISSUE,
                    " (WGD=", WGD_STATUS, ")"
                ),
                x = "MRCA age", y = "Cumulative frequency"
            )

        ggsave(
            paste0(
                CANCER_TISSUE,
                "_WGD=", WGD_STATUS, ".jpg"
            ),
            p,
            dpi = 150, width = 12, height = 6, units = "in"
        )
    }
}
