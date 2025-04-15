setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")

w_vec <- c(1, 0.7, 0.8)
v_vec <- c(0.4, 0.3, 0)
lambda_vec <- w_vec + 2 * v_vec
u_vec <- c(v_vec[1:(length(v_vec) - 1)] / lambda_vec[1:(length(lambda_vec) - 1)], 0)

alpha <- 1.0
r_initial <- 100000
max_time <- 100
n_reps <- 100
tau <- 0.01

# 7) RUN REPLICATES, ODE, VARIANCE, THEN PLOT
all_sims_new <- run_replicates(
    n_reps, r_initial, lambda_vec, u_vec, alpha, max_time, tau,
    seed = 123
)

mean_ode_df <- mean_ODE_solver(lambda_vec, u_vec, max_time, dt = 0.1)
the_var_df <- variance_ODE_solver(mean_ode_df, lambda_vec, u_vec, alpha)

rep_melt <- reshape2::melt(
    all_sims_new,
    id.vars = c("time", "replicate"),
    measure.vars = grep("^Nbar", names(all_sims_new), value = TRUE),
    variable.name = "compartment",
    value.name = "value"
)
rep_melt$source <- "replicate"

avg_sims <- all_sims_new %>%
    group_by(time) %>%
    summarize(across(starts_with("Nbar"), mean), .groups = "drop")

mean_melt <- reshape2::melt(
    avg_sims,
    id.vars = "time",
    measure.vars = grep("^Nbar", names(avg_sims), value = TRUE),
    variable.name = "compartment",
    value.name = "value"
)
mean_melt$source <- "mean"

ode_melt <- reshape2::melt(
    mean_ode_df,
    id.vars = "time",
    measure.vars = grep("^Nbar", names(mean_ode_df), value = TRUE),
    variable.name = "compartment",
    value.name = "value"
)
ode_melt$source <- "ODE"

lines_df <- bind_rows(rep_melt, mean_melt, ode_melt)
lines_df$compartment <- sub("Nbar", "N", lines_df$compartment)

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
    labs(title = "[A] Scaled population size", x = "Time", y = NULL) +
    theme_minimal() +
    theme(legend.position = "none")

# Variance of difference
f_list <- list()
comp_cols <- grep("^Nbar", names(mean_ode_df), value = TRUE)
for (cn in comp_cols) {
    f_list[[cn]] <- approxfun(mean_ode_df$time, mean_ode_df[[cn]], rule = 2)
}

all_sims_nhat <- all_sims_new %>%
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

var_long <- data.frame()
n_types <- length(lambda_vec)
for (i in seq_len(n_types)) {
    emp_col <- paste0("varEmp_Nhat_Nbar", i - 1)
    the_col <- paste0("varNhat", i - 1, "_the")
    cptname <- paste0("N", i - 1)

    var_long <- rbind(
        var_long,
        data.frame(
            time = var_nhat_emp$time,
            compartment = cptname,
            source = "Empirical",
            value = var_nhat_emp[[emp_col]]
        ),
        data.frame(
            time = the_var_df$time,
            compartment = cptname,
            source = "Theory",
            value = the_var_df[[the_col]]
        )
    )
}

p_bottom <- ggplot() +
    geom_line(
        data = subset(var_long, source == "Empirical"),
        aes(x = time, y = value, color = source),
        size = 2
    ) +
    geom_line(
        data = subset(var_long, source == "Theory"),
        aes(x = time, y = value, color = source),
        size = 2, linetype = "dashed"
    ) +
    facet_wrap(~compartment, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c("Empirical" = "#AA0000", "Theory" = "#0080FF")) +
    labs(
        title = "[B] Variance of difference in scaled population size",
        x = "Time", y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")

final_plot_new <- p_top / p_bottom
png("Illustration.png", res = 300, width = 16, height = 12, units = "in")
print(final_plot_new)
dev.off()
