#setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/abcsmcrf/R")
setwd("/Users/keitotaketomi/Downloads/abc-smc-rf 2/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

#setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R")
setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)

#setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")
setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/vignettes")



###############################################################################
# 2) LOAD THE REAL PATIENT TIMING TABLE
###############################################################################
real_df <- read.table(
    #"/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/DATASETS/PCAWG/evolution_and_heterogeneity/2018-07-24-wgdMrcaTiming.txt",
    "/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.txt",
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
    print(nrow(parameters))
    one_parameter <- function(parameter) {
        lambda <- as.numeric(parameter$lambda)
        #---Model parameters
        r_initial <- c(100000000, 1)
        lambda_vec <- c(0.00001, lambda)
        threshold_diagnosis <- 0.001
        alpha <- 1
        #---Fixed parameters
        u_vec <- c(0, 0)
        tau <- 0.01
        n_simulations <- 100
        max_time <- histogram_x[length(histogram_x)]
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

###############################################################################
# 3) SET UP PRIOR, DENSITY, AND A UNIVERSAL smcrf‐COMPATIBLE WRAPPER
###############################################################################

# (a) extract bounds from the lambda grid
lambda_min <- 0
lambda_max <- 10

# (b) prior sampler: must be called rprior(Nparameters)
rprior <- function(Nparameters) {
    data.frame(
        lambda = runif(Nparameters, min = lambda_min, max = lambda_max)
    )
}

# (c) density of that prior: smcrf will pass a data.frame with column 'lambda'
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "lambda")) {
        probs <- probs * dunif(parameters[["lambda"]], min = lambda_min, max = lambda_max)
    }
    return(probs)
}

###############################################################################
# 4) RUN ABC‐SMC‐RF & PLOT POSTERIOR MARGINAL OVER λ
###############################################################################

# hyperparameters for the ABC‐SMC‐RF
NUM_PARTICLES <- 1000
NUM_ITERATIONS <- 5
NUM_TREES <- 500

# # run the single‐parameter ABC‐SMC‐RF
# smcrf_results <- smcrf(
#     method            = "smcrf-single-param",
#     statistics_target = statistics_target,
#     model             = model,
#     rprior            = rprior,
#     dprior            = dprior,
#     nParticles        = rep(NUM_PARTICLES, NUM_ITERATIONS),
#     ntrees            = NUM_TREES,
#     parallel          = TRUE
# )

# plot_compare_marginal(
#     abc_results = smcrf_results,
#     # parameters_labels = parameters_labels,
#     plot_statistics = FALSE,
#     plot_hist = TRUE,
#     plot_prior = TRUE
# )

# final_parameters <- smcrf_results[[paste0("Iteration_", NUM_ITERATIONS + 1)]]$parameters
# final_statistics <- smcrf_results[[paste0("Iteration_", NUM_ITERATIONS + 1)]]$statistics

final_statistics<-model(
    parameters = data.frame(lambda=rep(14,10))
)

p <- plot_age_distribution(
    statistics_simulated = final_statistics,
    statistics_target = statistics_target,
    plot_title = paste0(CANCER_TISSUE, " (WGD=", WGD_STATUS, "): sim mean±SD vs. target")
)
png(paste0(CANCER_TISSUE, "_WGD=", WGD_STATUS, "_results.png"), res = 150, width = 30, height = 15, units = "in", pointsize = 12)
print(p)
dev.off()

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) pivot into long form, drop the NA‐bin, and parse age
sim_long <- final_statistics %>%
  mutate(sim_id = row_number()) %>% 
  pivot_longer(
    cols     = starts_with("Age_group_"),
    names_to = "bin",
    values_to= "freq"
  ) %>%
  filter(bin != "Age_group_NA") %>%
  mutate(
    age = as.integer(sub("Age_group_", "", bin))
  )

# 2) compute mean ± SD at each age
sim_summary <- sim_long %>%
  group_by(age) %>%
  dplyr::summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    sd_freq   = sd(freq,   na.rm = TRUE),
    .groups   = "drop"
  )

# 3) reshape target in the same way
target_long <- statistics_target %>%
  pivot_longer(
    cols     = everything(),
    names_to = "bin",
    values_to= "freq"
  ) %>%
  filter(bin != "Age_group_NA") %>%
  mutate(
    age = as.integer(sub("Age_group_", "", bin))
  )

# 4) plot
p <- ggplot() +
  geom_ribbon(
    data = sim_summary,
    aes(x = age,
        ymin = mean_freq - sd_freq,
        ymax = mean_freq + sd_freq),
    fill  = "steelblue",
    alpha = 0.3
  ) +
  geom_line(
    data = sim_summary,
    aes(x = age, y = mean_freq),
    color = "steelblue",
    size  = 1
  ) +
  geom_line(
    data  = target_long,
    aes(x = age, y = freq),
    color = "firebrick",
    size  = 1.2
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x     = "Age of diagnosis (years)",
    y     = "Normalized frequency",
    title = paste0(
      CANCER_TISSUE,
      " (WGD=", WGD_STATUS, 
      "): Simulated mean ± SD vs. Target"
    )
  )

print(p)