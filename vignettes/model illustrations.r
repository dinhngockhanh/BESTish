setwd("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R/")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd("/Users/kndinh/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes/")

plot_nsimulations <- 100
plot_tau <- 0.01
max_time <- 100
parameters <- data.frame(
    w0 = 1,
    w1 = c(1, 1.05, 1, 1.05),
    v0 = 5 * 10^-4,
    alpha = c(1, 1, 0.985, 0.985),
    R = 11000,
    ID = c("neutral_stable", "selective_stable", "neutral_expanding", "selective_expanding")
)
FORCE_WHITE <- theme(
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "white")
)
for (i in 1:nrow(parameters)) {
    w0 <- parameters$w0[i]
    w1 <- parameters$w1[i]
    v0 <- parameters$v0[i]
    alpha <- parameters$alpha[i]
    R <- parameters$R[i]
    ID <- parameters$ID[i]
    print(ID)
    #---Simulate trajectories with MAP parameters
    lambda0 <- w0 + 2 * v0
    lambda1 <- w1
    u0 <- v0 / lambda0
    u1 <- 0
    map_simulations <- run_replicates(
        n_reps = plot_nsimulations,
        r = R,
        lambda_vec = c(lambda0, lambda1),
        u_vec = c(u0, u1),
        alpha = alpha,
        max_time = max_time,
        tau = plot_tau,
        seed = 123
    )
    map_simulations$vaf <- 0.5 * map_simulations$N1 / (map_simulations$N0 + map_simulations$N1)
    map_simulations_mean <- map_simulations %>%
        group_by(time) %>%
        summarise(across(where(is.numeric) & !matches("replicate"), mean, na.rm = TRUE))
    #---Compute mean and variance trajectories with MAP parameters
    map_mean_var_df <- compute_expected_means_and_variances(
        time_grids = seq(0, max_time, by = plot_tau),
        h = plot_tau,
        w0 = w0,
        w1 = w1,
        v0 = v0,
        alpha = alpha
    )
    map_mean_var_df$N0_lowerCI <- map_mean_var_df$N0_bar - 1.96 * sqrt(map_mean_var_df$V11 / R)
    map_mean_var_df$N0_upperCI <- map_mean_var_df$N0_bar + 1.96 * sqrt(map_mean_var_df$V11 / R)
    map_mean_var_df$N1_lowerCI <- map_mean_var_df$N1_bar - 1.96 * sqrt(map_mean_var_df$V22 / R)
    map_mean_var_df$N1_upperCI <- map_mean_var_df$N1_bar + 1.96 * sqrt(map_mean_var_df$V22 / R)
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
    #---Plot population size trajectories
    p <- ggplot() +
        geom_ribbon(
            data = map_mean_var_df,
            aes(x = time, ymin = pmax(0, N0_lowerCI), ymax = pmin(3, N0_upperCI)),
            fill = "#3182BD", alpha = 0.3
        ) +
        geom_line(data = map_simulations, aes(x = time, y = Nbar0, group = replicate), color = "#3182BD", size = 1, alpha = 0.2) +
        geom_line(data = map_mean_var_df, aes(x = time, y = N0_bar), color = "#3182BD", size = 5) +
        geom_ribbon(
            data = map_mean_var_df,
            aes(x = time, ymin = pmax(0, N1_lowerCI), ymax = pmin(3, N1_upperCI)),
            fill = "#E6550D", alpha = 0.3
        ) +
        geom_line(data = map_simulations, aes(x = time, y = Nbar1, group = replicate), color = "#E6550D", size = 1, alpha = 0.2) +
        geom_line(data = map_mean_var_df, aes(x = time, y = N1_bar), color = "#E6550D", size = 5) +
        labs(
            x = "Age",
            y = "Scaled population size",
            title = NULL
        ) +
        xlim(0, max_time) +
        ylim(0, 3) +
        theme_minimal(base_size = 50) +
        FORCE_WHITE
    png(file.path(paste0("Illustration_", ID, "_population_size.png")), width = 30, height = 10, units = "in", res = 300)
    print(p)
    dev.off()
    #---Plot VAF trajectories
    p <- ggplot() +
        geom_ribbon(
            data = map_mean_var_df,
            aes(x = time, ymin = pmax(0, vaf_lowerCI), ymax = pmin(0.5, vaf_upperCI)),
            fill = "#666666", alpha = 0.3
        ) +
        geom_line(data = map_simulations, aes(x = time, y = vaf, group = replicate), color = "#666666", size = 1, alpha = 0.2) +
        geom_line(data = map_mean_var_df, aes(x = time, y = vaf_mu), color = "#666666", size = 5) +
        labs(
            x = "Age",
            y = "Variant Allele Frequency",
            title = NULL
        ) +
        xlim(0, max_time) +
        ylim(0, 0.5) +
        theme_minimal(base_size = 50) +
        FORCE_WHITE
    png(file.path(paste0("Illustration_", ID, "_vaf.png")), width = 30, height = 10, units = "in", res = 300)
    print(p)
    dev.off()
}
