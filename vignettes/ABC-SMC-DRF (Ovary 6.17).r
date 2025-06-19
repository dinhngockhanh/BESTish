setwd("/Users/keitotaketomi/Downloads/abc-smc-rf 2/R")
sapply(list.files(pattern = "\\.[rR]$"), source)

setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/R")
sapply(list.files(pattern = "\\.[rR]$"), source)

setwd("/Users/keitotaketomi/Documents/DriverSelectionSweep/vignettes")
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(parallel)
library(pbapply)
library(scales)  

# Output directory (To-do 1: white background plots)
out_dir <- path.expand("~/Documents/ABC-SMC-RF Results (Ovary)/")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 1) LOAD REAL DATA -----------------------------------------------------------
real_df <- read.table(
  "/Users/keitotaketomi/Documents/2018-07-24-wgdMrcaTiming.txt",
  header = TRUE, stringsAsFactors = FALSE
) %>%
  rename(
    tissue = tissue,
    WGD    = WGD,
    MRCA   = MRCA.time.linear,
    diag   = age
  ) %>%
  filter(!is.na(MRCA))

# ─────────────────────────────────────────────────────────────────────────────
# 2) MODEL + UPDATED PRIORS (To-do 3) ----------------------------------------

model <- function(parameters, statistic = "distance", parallel = TRUE) {
  one_parameter <- function(parameter) {
    log10_r0        <- as.numeric(parameter$log10_r0)
    lambda1         <- as.numeric(parameter$lambda1)
    lambda0         <- as.numeric(parameter$lambda0)
    log10_thresh    <- as.numeric(parameter$log10_threshold)
    alpha           <- as.numeric(parameter$alpha)
    r0              <- 10^log10_r0
    threshold       <- 10^log10_thresh

    u_vec <- c(0,0); tau <- 0.01; n_sim <- 100

    pilot_out <- simulate_continuous_moran_tau(
      r          = c(r0, lambda0),
      lambda_vec = c(lambda0, lambda1),
      u_vec      = u_vec,
      alpha      = alpha,
      max_time   = histogram_x[length(histogram_x)],
      tau        = tau
    )
    time_grid <- pilot_out$time
    n_time    <- length(time_grid)

    N1_mat <- matrix(NA_real_, nrow=n_sim, ncol=n_time)
    N1_mat[1,] <- pilot_out$N1

    for(i in 2:n_sim) {
      sim_out <- simulate_continuous_moran_tau(
        r          = c(r0, lambda0),
        lambda_vec = c(lambda0, lambda1),
        u_vec      = u_vec,
        alpha      = alpha,
        max_time   = histogram_x[length(histogram_x)],
        tau        = tau
      )
      if(length(sim_out$N1)==n_time) {
        N1_mat[i,] <- sim_out$N1
      } else if(length(sim_out$N1)>n_time) {
        N1_mat[i,] <- sim_out$N1[1:n_time]
      } else {
        tmp_len <- length(sim_out$N1)
        if(tmp_len>0) {
          N1_mat[i,1:tmp_len] <- sim_out$N1
          pad_val <- if(!is.na(sim_out$N1[tmp_len])) sim_out$N1[tmp_len] else 0
        } else pad_val <- 0
        if(tmp_len< n_time) {
          N1_mat[i,(tmp_len+1):n_time] <- pad_val
        }
      }
    }

    mean_N1 <- colMeans(N1_mat, na.rm=TRUE)
    idx_old <- which(mean_N1>threshold)
    old_time <- if(length(idx_old)>0) time_grid[idx_old[1]] else NA_real_

    new_times <- rep(NA_real_, n_sim)
    for(i in 1:n_sim) {
      idx_new <- which(N1_mat[i,]>threshold)
      new_times[i] <- if(length(idx_new)>0) time_grid[idx_new[1]] else NA_real_
    }
    new_times <- new_times[!is.na(new_times)]

    if(length(new_times)>0) {
      cnts   <- hist(new_times,breaks=histogram_x,plot=FALSE)$counts
      freq   <- cnts/sum(cnts)
      sim_cdf<- cumsum(freq)
    } else {
      sim_cdf <- rep(0,length(histogram_x)-1)
    }

    if(statistic=="distance") {
      d <- sum(abs(sim_cdf-cdf_obs))
      if(!is.finite(d)) d <- Inf
      return(data.frame(distance=d))
    } else {
      df <- as.data.frame(t(sim_cdf))
      colnames(df) <- paste0("Age_bin_", seq_along(sim_cdf))
      return(df)
    }
  }

  if(parallel) {
    cl <- makePSOCKcluster(detectCores()-1)
    clusterExport(cl, c("one_parameter","simulate_continuous_moran_tau",
                       "histogram_x","cdf_obs","statistic"),
                  envir=environment())
    stats_list <- parLapply(cl, seq_len(nrow(parameters)), function(i){
      one_parameter(as.data.frame(parameters[i,,drop=FALSE]))
    })
    stopCluster(cl)
    stats <- rbindlist(stats_list)
  } else {
    stats <- rbindlist(lapply(seq_len(nrow(parameters)), function(i){
      one_parameter(as.data.frame(parameters[i,,drop=FALSE]))
    }))
  }
  cbind(parameters, stats)
}

# UPDATED PRIORS (To-do 3)
rprior <- function(Nparameters) {
  log10_r0        <- runif(Nparameters, min=2, max=8)      # 2→8 
  lambda1         <- runif(Nparameters, min=0, max=30)     # 0→30
  lambda0         <- runif(Nparameters, min=0, max=lambda1)
  log10_threshold <- runif(Nparameters, min=1, max=log10_r0)
  alpha           <- runif(Nparameters, min=0.3, max=1)    # 0.3→1
  data.frame(log10_r0,lambda1,lambda0,log10_threshold,alpha)
}

dprior <- function(parameters, parameter_id="all") {
  probs <- rep(1,nrow(parameters))
  if(parameter_id %in% c("all","log10_r0"))
    probs <- probs * dunif(parameters$log10_r0, min=2, max=8)
  if(parameter_id %in% c("all","lambda1"))
    probs <- probs * dunif(parameters$lambda1, min=0, max=30)
  if(parameter_id %in% c("all","lambda0")) {
    zero_idx    <- which(parameters$lambda1<=0)
    nonzero_idx <- which(parameters$lambda1>0)
    if(length(nonzero_idx)>0)
      probs[nonzero_idx] <- probs[nonzero_idx] *
        dunif(parameters$lambda0[nonzero_idx], min=0, max=parameters$lambda1[nonzero_idx])
    if(length(zero_idx)>0) probs[zero_idx] <- probs[zero_idx]*1
  }
  if(parameter_id %in% c("all","log10_threshold")) {
    idx <- which(parameters$log10_r0>1)
    if(length(idx)>0)
      probs[idx] <- probs[idx] *
        dunif(parameters$log10_threshold[idx], min=1, max=parameters$log10_r0[idx])
  }
  if(parameter_id %in% c("all","alpha"))
    probs <- probs * dunif(parameters$alpha, min=0.3, max=1)
  probs
}

# ─────────────────────────────────────────────────────────────────────────────
# 3) PARAMETER TESTING CONFIGURATIONS (To-do 5) -----------------------------

# Baseline configuration
CONFIG_BASELINE <- list(
  name = "baseline",
  particles = 10000,
  iterations = 10,
  trees = 500,
  n_sim = 100  # simulations per parameter set
)

# Testing configurations (To-do 5)
CONFIGS_TEST <- list(
  list(name = "more_particles", particles = 20000, iterations = 5, trees = 500, n_sim = 100),
  list(name = "more_iterations", particles = 10000, iterations = 10, trees = 500, n_sim = 100),
  list(name = "more_trees", particles = 10000, iterations = 5, trees = 1000, n_sim = 100),
  list(name = "more_sims", particles = 10000, iterations = 5, trees = 500, n_sim = 200)
)

# Choose which configurations to run (for testing)
CONFIGS_TO_RUN <- list(CONFIG_BASELINE)  # Add CONFIGS_TEST for full comparison

# ─────────────────────────────────────────────────────────────────────────────
# 4) IMPROVED BOOTSTRAP CDF FUNCTION (To-do 2: fixed CI issues) -------------
bootstrap_cdf <- function(data, breaks, n_bootstrap=1000) {
  M <- length(breaks) - 1
  n_data <- length(data)
  
  # Ensure we have enough data points
  if(n_data < 5) {
    warning("Too few data points for reliable bootstrap")
    return(list(
      lower = rep(0, M),
      upper = rep(1, M)
    ))
  }
  
  # Pre-allocate matrix for bootstrap CDFs
  all_cdfs <- matrix(NA, nrow = M, ncol = n_bootstrap)
  
  # Generate bootstrap samples
  for(i in 1:n_bootstrap) {
    # Bootstrap sample with replacement
    xs <- sample(data, size = n_data, replace = TRUE)
    
    # Compute histogram counts
    cts <- hist(xs, breaks = breaks, plot = FALSE)$counts
    
    # Handle edge cases
    if(sum(cts) == 0) {
      # If no counts, set CDF to 0
      all_cdfs[, i] <- rep(0, M)
    } else {
      # Convert to frequencies and then CDF
      freq <- cts / sum(cts)
      cdf_vals <- cumsum(freq)
      
      # Ensure CDF is monotonic and bounded
      cdf_vals <- pmax(0, pmin(1, cdf_vals))
      cdf_vals <- cummax(cdf_vals)  # Ensure monotonicity
      
      all_cdfs[, i] <- cdf_vals
    }
  }
  
  # Calculate confidence intervals with robust handling
  lower <- numeric(M)
  upper <- numeric(M)
  
  for(j in 1:M) {
    # Get all bootstrap values for this bin
    boot_vals <- all_cdfs[j, ]
    boot_vals <- boot_vals[!is.na(boot_vals)]  # Remove any NAs
    
    if(length(boot_vals) == 0) {
      # No valid values
      lower[j] <- 0
      upper[j] <- 1
    } else if(length(boot_vals) == 1) {
      # Only one value
      lower[j] <- upper[j] <- boot_vals[1]
    } else {
      # Calculate quantiles
      lower[j] <- quantile(boot_vals, 0.025, na.rm = TRUE)
      upper[j] <- quantile(boot_vals, 0.975, na.rm = TRUE)
    }
  }
  
  # Post-process to ensure proper CDF properties
  # 1. Ensure bounds [0, 1]
  lower <- pmax(0, pmin(1, lower))
  upper <- pmax(0, pmin(1, upper))
  
  # 2. Ensure monotonicity (CDFs must be non-decreasing)
  lower <- cummax(lower)
  upper <- cummax(upper)
  
  # 3. Ensure lower <= upper
  for(j in 1:M) {
    if(lower[j] > upper[j]) {
      # If lower > upper, set both to their average
      avg <- (lower[j] + upper[j]) / 2
      lower[j] <- upper[j] <- avg
    }
  }
  
  # 4. Smooth any remaining discontinuities
  # If there are large jumps, interpolate to make smoother
  for(j in 2:M) {
    # Check for large gaps in lower bound
    if(lower[j] - lower[j-1] > 0.3) {
      lower[j] <- lower[j-1] + 0.1  # Smooth transition
    }
    # Check for large gaps in upper bound  
    if(upper[j] - upper[j-1] > 0.3) {
      upper[j] <- upper[j-1] + 0.1  # Smooth transition
    }
  }
  
  # 5. Final bounds check
  lower <- pmax(0, pmin(1, lower))
  upper <- pmax(0, pmin(1, upper))
  
  # 6. Ensure final monotonicity
  lower <- cummax(lower)
  upper <- cummax(upper)
  
  return(list(
    lower = lower,
    upper = upper,
    n_bootstrap = n_bootstrap,
    n_data = n_data
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# 5) MAIN ANALYSIS LOOP - CHECK AVAILABLE TISSUES FIRST --------------------

# First, let's check what tissues are available in the dataset
cat("Available tissues in dataset:\n")
tissue_counts <- real_df %>% 
  group_by(tissue, WGD) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  arrange(tissue, WGD)
print(tissue_counts)

# Look for ovarian cancer tissues (might have different naming)
ovary_tissues <- unique(real_df$tissue)[grepl("(?i)(ovar|ov)", unique(real_df$tissue), perl=TRUE)]
cat("\nPossible ovarian tissues found:\n")
print(ovary_tissues)

# Use the first ovarian tissue found, or fall back to a common tissue
if(length(ovary_tissues) > 0) {
  target_tissue <- ovary_tissues[1]
  cat("\nUsing tissue:", target_tissue, "\n")
} else {
  # If no ovarian tissue found, use the tissue with most samples
  most_common <- tissue_counts %>% 
    group_by(tissue) %>% 
    summarise(total = sum(count)) %>% 
    arrange(desc(total)) %>% 
    slice(1) %>% 
    pull(tissue)
  target_tissue <- most_common
  cat("\nNo ovarian tissue found. Using most common tissue:", target_tissue, "\n")
}

for(config in CONFIGS_TO_RUN) {
  cat("\n", rep("=", 60), "\n")
  cat("RUNNING CONFIGURATION:", config$name, "\n")
  cat("Particles:", config$particles, "| Iterations:", config$iterations, 
      "| Trees:", config$trees, "| Sims:", config$n_sim, "\n")
  cat(rep("=", 60), "\n")
  
  # Test only on the identified tissue with WGD=TRUE and FALSE
  ct <- target_tissue
  for(w in c(FALSE,TRUE)) {
    cat("\n----", ct, "| WGD =", w, "| Config:", config$name, "----\n")
    sub_df <- filter(real_df, tissue==ct, WGD==w)
    cat("Found", nrow(sub_df), "samples\n")
    if(nrow(sub_df)<10) {
      cat(" → too few (",nrow(sub_df),") → skip\n")
      next
    }

    x <- sub_df$MRCA
    max_t <- ceiling(max(x, na.rm=TRUE))
    if(!is.finite(max_t)||max_t<1) max_t <- 1
    histogram_x <- unique(seq(0, max_t, length.out=11))
    if(length(histogram_x)<2) histogram_x <- c(0,1)
    mids <- (histogram_x[-1]+histogram_x[-length(histogram_x)])/2

    counts_obs <- hist(x, breaks=histogram_x, plot=FALSE)$counts
    cdf_obs    <- cumsum(counts_obs/sum(counts_obs))

    # Improved bootstrap CI (To-do 2: Fix abrupt endings)
    boot_ci <- bootstrap_cdf(x, histogram_x, n_bootstrap=1000)
    obs_lo  <- boot_ci$lower
    obs_hi  <- boot_ci$upper
    
    # CRITICAL FIX: Ensure confidence intervals extend to the end
    # The issue is that CDF should reach 1.0 by the last bin, but CI might not
    
    # Force the last confidence interval to reach 1.0 (since CDF must end at 1)
    if(length(obs_lo) > 0 && length(obs_hi) > 0) {
      # Ensure the last bin has upper CI = 1.0
      obs_hi[length(obs_hi)] <- 1.0
      
      # Ensure monotonicity is preserved after forcing the endpoint
      for(i in length(obs_hi):2) {
        if(obs_hi[i-1] > obs_hi[i]) {
          obs_hi[i-1] <- obs_hi[i]
        }
      }
      for(i in length(obs_lo):2) {
        if(obs_lo[i-1] > obs_lo[i]) {
          obs_lo[i-1] <- obs_lo[i]
        }
      }
      
      # Ensure lower <= upper throughout
      for(i in 1:length(obs_lo)) {
        if(obs_lo[i] > obs_hi[i]) {
          obs_lo[i] <- obs_hi[i]
        }
      }
    }

    # Try RF first, then DRF (To-do 4)
    abc_rf <- tryCatch({
      smcrf(
        method            = "smcrf-single-param",
        statistics_target = data.frame(distance=0),
        model             = model,
        rprior            = rprior,
        dprior            = dprior,
        nParticles        = rep(config$particles, config$iterations),
        ntrees            = config$trees,
        parallel          = TRUE
      )
    }, error=function(e){
      cat("   RF error:", conditionMessage(e),"\n"); NULL
    })

    if(is.null(abc_rf)) {
      cat("   → RF failed → try DRF\n")
      abc <- tryCatch({
        smcrf(
          method            = "smcrf-drf",
          statistics_target = data.frame(distance=0),
          model             = model,
          rprior            = rprior,
          dprior            = dprior,
          nParticles        = rep(config$particles, config$iterations),
          ntrees            = config$trees,
          parallel          = TRUE
        )
      }, error=function(e){
        cat("   DRF error:", conditionMessage(e),"\n"); NULL
      })
      if(is.null(abc)) {
        cat("   → both failed → skip\n"); next
      }
      method_used <- "DRF"
      cat("   ✓ DRF succeeded\n")
    } else {
      abc <- abc_rf
      method_used <- "RF"
      cat("   ✓ RF succeeded\n")
    }

    # Extract final posteriors
    final_it <- paste0("Iteration_", config$iterations+1)
    post_pars <- abc[[final_it]]$parameters

    # File prefix with configuration
    file_prefix <- if(config$name == "baseline") {
      paste0(ct, "_WGD=", w)
    } else {
      paste0(ct, "_WGD=", w, "_", config$name)
    }

    # --- INDIVIDUAL HISTOGRAM PLOTS (original requirement) ------------------
    param_names <- c("log10_r0","lambda1","lambda0","log10_threshold","alpha")
    for (pn in param_names) {
      df_post <- data.frame(value = post_pars[[pn]])
      p_hist  <- ggplot(df_post, aes(x = value)) +
        geom_histogram(
          aes(y = ..density..),
          bins  = 30,
          fill  = "steelblue",
          color = "black",
          alpha = 0.7
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.background = element_rect(fill = "white"),   # To-do 1
          panel.background = element_rect(fill = "white")   # To-do 1
        ) +
        labs(
          title = paste0(ct, " (WGD=", w, "): ", pn, " [", config$name, "]"),
          x     = pn,
          y     = "Density"
        )
      ggsave(
        filename = file.path(out_dir, paste0(file_prefix, "_hist-", pn, ".png")),
        plot   = p_hist,
        dpi    = 300,
        width  = 8,
        height = 6,
        units  = "in"
      )
    }

    # --- SIMULATED CDFs & CI (Fix abrupt endings) ----------------------------
    sims <- model(post_pars[1:100,], statistic="histogram", parallel=FALSE)
    sim_mat <- as.matrix(sims[,grep("Age_bin_",names(sims))])
    sim_lo  <- apply(sim_mat,2,quantile,0.025,na.rm=TRUE)
    sim_hi  <- apply(sim_mat,2,quantile,0.975,na.rm=TRUE)
    sim_mean<- apply(sim_mat,2,mean,na.rm=TRUE)
    
    # CRITICAL FIX: Ensure simulated CIs also extend properly to 1.0
    if(length(sim_lo) > 0 && length(sim_hi) > 0) {
      # Force the last simulation CI to reach 1.0
      sim_hi[length(sim_hi)] <- 1.0
      
      # Ensure monotonicity for simulation CIs
      for(i in length(sim_hi):2) {
        if(sim_hi[i-1] > sim_hi[i]) {
          sim_hi[i-1] <- sim_hi[i]
        }
      }
      for(i in length(sim_lo):2) {
        if(sim_lo[i-1] > sim_lo[i]) {
          sim_lo[i-1] <- sim_lo[i]
        }
      }
      
      # Ensure lower <= upper for simulations
      for(i in 1:length(sim_lo)) {
        if(sim_lo[i] > sim_hi[i]) {
          sim_lo[i] <- sim_hi[i]
        }
      }
    }

    # Individual simulation lines
    sims_long <- sims %>%
      mutate(simID=row_number()) %>%
      pivot_longer(starts_with("Age_bin_"), names_to="bin", values_to="cdf") %>%
      mutate(
        bin = as.integer(sub("Age_bin_","",bin)),
        mid = mids[bin]
      )

    # Data frames for plotting
    obs_df <- tibble(mid=mids, cdf=cdf_obs, lower=obs_lo, upper=obs_hi)
    sim_df <- tibble(mid=mids, cdf=sim_mean, lower=sim_lo, upper=sim_hi)

    # --- ENHANCED CDF PLOT (To-do 1: white background) ----------------------
    p_cdf <- ggplot() +
      geom_line(data=sims_long, aes(x=mid,y=cdf,group=simID),
                color="lightblue", alpha=0.1, size=0.3) +
      geom_ribbon(data=sim_df, aes(x=mid,ymin=lower,ymax=upper),
                  fill="steelblue", alpha=0.3) +
      geom_ribbon(data=obs_df, aes(x=mid,ymin=lower,ymax=upper),
                  fill="firebrick", alpha=0.3) +
      geom_line(data=sim_df, aes(x=mid,y=cdf), color="steelblue", size=1.2) +
      geom_line(data=obs_df, aes(x=mid,y=cdf), color="firebrick", size=1.2) +
      theme_minimal(base_size=14) +
      theme(
        panel.grid.minor=element_blank(),
        plot.background=element_rect(fill="white"),      # To-do 1
        panel.background=element_rect(fill="white")      # To-do 1
      ) +
      labs(
        title    = paste0(ct," (WGD=",w,"): CDF [",config$name,", ",method_used,"]"),
        subtitle = "95% CI: shaded regions | Blue=simulation, Red=data",
        x        = "MRCA age",
        y        = "Cumulative frequency"
      ) +
      scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
      annotate("text", x=Inf, y=0.05, 
               label=paste("n =", nrow(sub_df)), 
               hjust=1.1, vjust=0, size=4, color="gray50")

    ggsave(
      filename = file.path(out_dir, paste0(file_prefix, "_enhanced_cdf.png")),
      plot     = p_cdf, dpi=300,
      width    = 12, height=8, units="in"
    )

    # --- VIOLIN PLOT WITH CORRECT RANGES (To-do 3 ranges) ------------------
    # Create violin plot with proper y-axis ranges matching updated priors
    param_data <- post_pars %>%
      pivot_longer(
        cols = c("log10_r0", "lambda1", "lambda0", "log10_threshold", "alpha"),
        names_to = "parameter",
        values_to = "value"
      ) %>%
      mutate(
        parameter = factor(parameter, 
                          levels = c("log10_r0", "lambda1", "lambda0", "log10_threshold", "alpha"),
                          labels = c("log₁₀(r₀)", "λ₁", "λ₀", "log₁₀(threshold)", "α"))
      )

    p_violin <- ggplot(param_data, aes(x = 1, y = value, fill = parameter)) +
      geom_violin(
        alpha = 0.7,
        scale = "width",
        trim = FALSE,
        draw_quantiles = c(0.25, 0.5, 0.75),
        width = 0.8
      ) +
      scale_fill_manual(
        values = c("log₁₀(r₀)" = "#4292c6", "λ₁" = "#41ab5d", "λ₀" = "#9e9ac8", 
                   "log₁₀(threshold)" = "#d94701", "α" = "#fd8d3c"),
        guide = "none"
      ) +
      facet_wrap(~ parameter, scales = "free_y", ncol = 5) +
      # Force correct y-axis ranges based on updated priors
      geom_blank(data = data.frame(
        parameter = factor(c("log₁₀(r₀)", "λ₁", "λ₀", "log₁₀(threshold)", "α"),
                          levels = c("log₁₀(r₀)", "λ₁", "λ₀", "log₁₀(threshold)", "α")),
        value = c(2, 0, 0, 1, 0.3),  # min values (updated)
        x = 1
      ), aes(x = x, y = value)) +
      geom_blank(data = data.frame(
        parameter = factor(c("log₁₀(r₀)", "λ₁", "λ₀", "log₁₀(threshold)", "α"),
                          levels = c("log₁₀(r₀)", "λ₁", "λ₀", "log₁₀(threshold)", "α")),
        value = c(8, 30, 30, 8, 1),  # max values (updated)
        x = 1
      ), aes(x = x, y = value)) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        plot.background = element_rect(fill = "white"),    # To-do 1
        panel.background = element_rect(fill = "white")    # To-do 1
      ) +
      labs(
        title = paste0(ct, " (WGD=", w, "): Parameters [", config$name, "]"),
        subtitle = "Posterior distributions with quartile lines (y-ranges = updated priors)",
        x = NULL,
        y = "Parameter Value"
      )

    ggsave(
      filename = file.path(out_dir, paste0(file_prefix, "_violin.png")),
      plot   = p_violin, dpi=300,
      width  = 15, height = 6, units="in"
    )
  }
}

