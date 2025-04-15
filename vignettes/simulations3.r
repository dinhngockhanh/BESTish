setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")

library(dplyr)

plotname <- "Illustration-1a"
w_vec <- c(0.5, 1.08, 1.1)
v_vec <- c(1, 0.1, 0)

# plotname <- "Illustration-1b"
# w_vec <- c(0.5, 0.08, 0.7)
# v_vec <- c(0.45, 0.55, 0)

# plotname <- "Illustration-1c"
# w_vec <- c(0.6, 0.01, 0.7)
# v_vec <- c(0.45, 0.55, 0)

# plotname <- "Illustration-2a"
# w_vec <- c(0.5, 1.1, 1)
# v_vec <- c(1, 0.6, 0)

# plotname <- "Illustration-2b"
# w_vec <- c(1.07, 1.1, 1)
# v_vec <- c(0.1, 0.2, 0)

# plotname <- "Illustration-2c"
# w_vec <- c(0.7, 1.1, 0.5)
# v_vec <- c(0.1, 0.2, 0)

lambda_vec <- w_vec + 2 * v_vec
u_vec <- c(v_vec[1:(length(v_vec) - 1)] / lambda_vec[1:(length(lambda_vec) - 1)], 0)

alpha <- 1.0
r_initial <- 100000
max_time <- 20
n_reps <- 1000
tau <- 0.01

theory_plot_times <- seq(0, max_time, by = 1)

########################################################################
simulations_raw <- run_replicates(
    n_reps = n_reps,
    r = r_initial,
    lambda_vec = lambda_vec,
    u_vec = u_vec,
    alpha = alpha,
    max_time = max_time,
    tau = tau,
    seed = 123
)

simulations_long <- reshape2::melt(
    simulations_raw,
    id.vars = c("time", "replicate"),
    measure.vars = grep("^Nbar", names(simulations_raw), value = TRUE),
    variable.name = "compartment",
    value.name = "value"
)
simulations_long$source <- "replicate"

theory_mean_raw <- mean_ODE_solver(lambda_vec, u_vec, max_time, dt = 0.1)
theory_variance_raw <- variance_ODE_solver(theory_mean_raw, lambda_vec, u_vec, alpha)

## Generalize for arbitrary number of compartments
comp_cols <- grep("^Nbar", names(simulations_raw), value = TRUE)

simulation_mean <- simulations_raw %>%
    group_by(time) %>%
    summarise_at(vars(all_of(comp_cols)), ~ mean(.x, na.rm = TRUE))

simulation_mean <- reshape2::melt(
    simulation_mean,
    id.vars = "time",
    measure.vars = comp_cols,
    variable.name = "compartment",
    value.name = "value"
)
simulation_mean$source <- "mean"

theory_mean <- reshape2::melt(
    theory_mean_raw,
    id.vars = "time",
    measure.vars = comp_cols,
    variable.name = "compartment",
    value.name = "value"
)
theory_mean$source <- "ODE"

theory_mean <- theory_mean %>%
    group_by(compartment) %>%
    filter(time %in% theory_plot_times | abs(time - theory_plot_times[which.min(abs(time - theory_plot_times))]) == min(abs(time - theory_plot_times))) %>%
    ungroup()

color_codes <- c(
    "N0" = "#B22C2C",
    "Nbar0" = "#B22C2C",
    "N1" = "#85B22C",
    "Nbar1" = "#85B22C",
    "N2" = "#2C85B2",
    "Nbar2" = "#2C85B2"
)

p_left <- ggplot()
for (rep in unique(simulations_raw$replicate)) {
    p_left <- p_left +
        geom_line(
            data = simulations_long %>% filter(replicate == rep),
            aes(x = time, y = value, color = compartment),
            alpha = 0.01
        )
}

p_left <- p_left +
    geom_line(
        data = simulation_mean,
        aes(x = time, y = value, color = compartment),
        size = 1
    ) +
    geom_point(
        data = theory_mean,
        aes(x = time, y = value, color = compartment),
        size = 5
    ) +
    scale_color_manual(values = color_codes) +
    labs(
        title = "Scaled population size",
        x = "Time", y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")



f_list <- list()
comp_cols <- grep("^Nbar", names(theory_mean_raw), value = TRUE)
for (cn in comp_cols) {
    f_list[[cn]] <- approxfun(theory_mean_raw$time, theory_mean_raw[[cn]], rule = 2)
}
all_sims_nhat <- simulations_raw %>%
    group_by(replicate) %>%
    mutate(
        across(
            all_of(comp_cols),
            ~ sqrt(r_initial) * (.x - f_list[[cur_column()]](time)),
            .names = "Nhat_{.col}"
        )
    )
nhat_cols <- grep("^Nhat_Nbar", names(all_sims_nhat), value = TRUE)
var_nhat_emp <- all_sims_nhat %>%
    group_by(time) %>%
    summarize(
        across(
            all_of(nhat_cols),
            ~ var(.x, na.rm = TRUE),
            .names = "varEmp_{.col}"
        ),
        .groups = "drop"
    )
simulation_variance <- data.frame()
theory_variance <- data.frame()
n_types <- length(lambda_vec)
for (i in seq_len(n_types)) {
    emp_col <- paste0("varEmp_Nhat_Nbar", i - 1)
    the_col <- paste0("varNhat", i - 1, "_the")
    comp_name <- paste0("N", i - 1)
    simulation_variance <- rbind(
        simulation_variance,
        data.frame(
            time = var_nhat_emp$time,
            compartment = comp_name,
            source = "Empirical",
            value = var_nhat_emp[[emp_col]]
        )
    )
    theory_variance <- rbind(
        theory_variance,
        data.frame(
            time = theory_variance_raw$time,
            compartment = comp_name,
            source = "Theory",
            value = theory_variance_raw[[the_col]]
        )
    )
}
theory_variance <- theory_variance %>%
    group_by(compartment) %>%
    filter(time %in% theory_plot_times | abs(time - theory_plot_times[which.min(abs(time - theory_plot_times))]) == min(abs(time - theory_plot_times))) %>%
    ungroup()

p_right <- ggplot() +
    geom_line(
        data = simulation_variance,
        aes(x = time, y = value, color = compartment),
        size = 1
    ) +
    geom_point(
        data = theory_variance,
        aes(x = time, y = value, color = compartment),
        size = 5
    ) +
    scale_color_manual(values = color_codes) +
    labs(
        title = "Variance",
        x = "Time", y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")

final_plot_new <- p_left | p_right
png(paste0(plotname, ".png"), res = 300, width = 16, height = 8, units = "in")
print(final_plot_new)
dev.off()
