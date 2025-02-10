# 1) Load Packages & Define Simulation Function
library(ggplot2)
library(reshape2)
library(dplyr)



setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/R")
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd("/Users/dinhngockhanh/My Drive (knd2127@columbia.edu)/RESEARCH AND EVERYTHING/Projects/GITHUB/DriverSelectionSweep/vignettes")



# 2) Define Multiple Scenarios of lambda_vec
scenario_params <- list(
    list(
        label = "Scenario A: (1.0, 0.8, 1.4)",
        lambda_vec = c(0.5, 0.4, 0.6)
    ),
    list(
        label = "Scenario B: (1.4, 1.2, 1.1)",
        lambda_vec = c(1.4, 1.2, 1.1)
    ),
    list(
        label = "Scenario C: (0.9, 1.2, 1.6)",
        lambda_vec = c(0.9, 1.2, 1.6)
    ),
    list(
        label = "Scenario D: (1.1, 1.1, 1.0)",
        lambda_vec = c(1.1, 1.1, 1.0)
    )
)

r_initial <- 100000
u_vec <- c(0.01, 0.01, 0.0) # two-hit: last type can't mutate further
alpha <- 1.0
max_time <- 100
tau <- 0.01

# 3) Run the Simulation for Each Scenario
all_results <- list()

for (sc in scenario_params) {
    # run the simulation
    sim_df <- simulate_continuous_moran_tau(
        r = r_initial,
        lambda_vec = sc$lambda_vec,
        u_vec = u_vec,
        alpha = alpha,
        max_time = max_time,
        tau = tau,
        # seed = 123
    )

    # store scenario label
    sim_df$scenario <- sc$label

    # accumulate
    all_results[[sc$label]] <- sim_df
}

# combine into one data frame
all_results_df <- dplyr::bind_rows(all_results)


# 4) Plot (All Scenarios in Facets or Overlaid)
long_df <- reshape2::melt(all_results_df,
    id.vars = c("time", "scenario"),
    variable.name = "Type", value.name = "Count"
)

# a) Facet by scenario, each with lines for N0, N1, N2
p_facet <- ggplot(long_df, aes(x = time, y = Count, color = Type)) +
    geom_line() +
    scale_y_continuous(trans = "log1p") + # or log10
    facet_wrap(~scenario, scales = "free_y") + # separate panels
    labs(
        title = "Continuous Moran (tau-leaping)",
        subtitle = "Multiple Lambda-Vec Scenarios",
        x = "Time", y = "Cell Count (log1p scale)"
    ) +
    theme_minimal()

print(p_facet)
