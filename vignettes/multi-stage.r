path_R <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R"
path_vignettes <- "/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes"
setwd(path_R)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(path_vignettes)

# ## 5) EXAMPLE: PRODUCE PLOTS [A] & [B]

## Example parameters for 3 compartments
## e.g. lambda0=0.5, lambda1=0.4, lambda2=0.8
lambda_vec <- c(0.3, 0.1, 0.5)
u_vec <- c(0.01, 0.01, 0.0) # last compartment has 0 mutation out
alpha <- 1.0
r_initial <- 500000
max_time <- 50
n_reps <- 100
tau <- 0.1

## (A) RUN REPLICATES
all_sims <- run_replicates(
    n_reps, r_initial, lambda_vec, u_vec,
    alpha, max_time, tau,
    seed = 123
)

## (B) SOLVE ODE FOR THE MEAN
mean_ode_df <- solve_mean_ode(lambda_vec, u_vec, max_time, dt = 0.1)

## (C) SOLVE THE VARIANCE ODE
# the_var_df <- solve_variance_ode(mean_ode_df, lambda_vec, u_vec, alpha)
the_var_df <- solve_variance_ode_OLD(mean_ode_df, lambda_vec, u_vec, alpha)

## 6) MAKE PLOT (A): REPLICATES vs. MEAN vs. ODE
rep_melt <- melt(
    all_sims,
    id.vars = c("time", "replicate"),
    measure.vars = c("Nbar0", "Nbar1", "Nbar2"),
    variable.name = "compartment",
    value.name = "value"
)
rep_melt$source <- "replicate"

avg_sims <- all_sims %>%
    group_by(time) %>%
    summarize(
        meanNbar0 = mean(Nbar0),
        meanNbar1 = mean(Nbar1),
        meanNbar2 = mean(Nbar2),
        .groups = "drop"
    )
mean_melt <- melt(
    avg_sims,
    id.vars = "time",
    measure.vars = c("meanNbar0", "meanNbar1", "meanNbar2"),
    variable.name = "compartment",
    value.name = "value"
)
mean_melt$source <- "mean"
mean_melt$compartment <- sub("meanNbar", "Nbar", mean_melt$compartment)

ode_melt <- melt(
    mean_ode_df,
    id.vars = "time",
    measure.vars = c("Nbar0", "Nbar1", "Nbar2"),
    variable.name = "compartment",
    value.name = "value"
)
ode_melt$source <- "ODE"

lines_df <- bind_rows(rep_melt, mean_melt, ode_melt)
lines_df$compartment <- sub("Nbar", "", lines_df$compartment)
lines_df$compartment <- factor(
    lines_df$compartment,
    levels = c("0", "1", "2"),
    labels = c("N0", "N1", "N2")
)

p_top <- ggplot() +
    geom_line(
        data = subset(lines_df, source == "replicate"),
        aes(x = time, y = value, group = replicate, color = source),
        alpha = 0.2
    ) +
    geom_line(
        data = subset(lines_df, source == "mean"),
        aes(x = time, y = value, color = source),
        size = 2
    ) +
    geom_line(
        data = subset(lines_df, source == "ODE"),
        aes(x = time, y = value, color = source),
        size = 2, linetype = "dashed"
    ) +
    facet_wrap(~compartment, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("replicate" = "#FFAAAA", "mean" = "#AA0000", "ODE" = "#0080FF")) +
    labs(
        title = "[A] Scaled population size",
        x = "Time", y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")

## 7) MAKE PLOT (B): Variance of difference in scaled population size
fN0 <- approxfun(mean_ode_df$time, mean_ode_df$Nbar0, rule = 2)
fN1 <- approxfun(mean_ode_df$time, mean_ode_df$Nbar1, rule = 2)
fN2 <- approxfun(mean_ode_df$time, mean_ode_df$Nbar2, rule = 2)

all_sims_nhat <- all_sims %>%
    group_by(replicate) %>%
    mutate(
        Nbar0_ode = fN0(time),
        Nbar1_ode = fN1(time),
        Nbar2_ode = fN2(time),
        Nhat0 = sqrt(r_initial) * (Nbar0 - Nbar0_ode),
        Nhat1 = sqrt(r_initial) * (Nbar1 - Nbar1_ode),
        Nhat2 = sqrt(r_initial) * (Nbar2 - Nbar2_ode)
    )

var_nhat_emp <- all_sims_nhat %>%
    group_by(time) %>%
    summarize(
        varNhat0_emp = var(Nhat0),
        varNhat1_emp = var(Nhat1),
        varNhat2_emp = var(Nhat2),
        .groups = "drop"
    )

var_long <- data.frame(time = c(), compartment = c(), source = c(), value = c())
for (i in 1:length(lambda_vec)) {
    var_long <- rbind(
        var_long,
        data.frame(
            time = var_nhat_emp$time,
            compartment = paste0("N", i - 1),
            source = "Empirical",
            value = var_nhat_emp[[paste0("varNhat", i - 1, "_emp")]]
        ),
        data.frame(
            time = the_var_df$time,
            compartment = paste0("N", i - 1),
            source = "Theory",
            value = the_var_df[[paste0("varNhat", i - 1, "_the")]]
        )
    )
}


p_bottom <- ggplot() +
    geom_line(
        data = subset(var_long, source == "Empirical"),
        aes(x = time, y = value, color = source), size = 2
    ) +
    geom_line(
        data = subset(var_long, source == "Theory"),
        aes(x = time, y = value, color = source), size = 2, linetype = "dashed"
    ) +
    facet_wrap(~compartment, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("Empirical" = "#AA0000", "Theory" = "#0080FF")) +
    labs(
        title = "[B] Variance of difference in scaled population size between one realization and the mean",
        x = "Time", y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")

final_plot <- p_top / p_bottom
print(final_plot)
