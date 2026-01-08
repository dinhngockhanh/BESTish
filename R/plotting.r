plot_BESTish_marginal <- function(grid,
                                  filename = NULL,
                                  filetype = "png",
                                  width = 30, height = 15, units = "in", res = 150) {
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(viridis))
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(gridExtra))
    #---Find variables to plot
    tmp <- sapply(grid, function(x) length(unique(x)))
    cols_to_plot <- intersect(names(tmp[tmp > 1]), c("w0", "w1", "log10v0", "alpha", "R"))
    cat(bold(yellow("Plotting marginal distributions from BESTish for variables:", paste0(cols_to_plot, collapse = ", "), "\n")))
    #---Create subplots for each marginal distribution
    plot_marginal <- function(df, col) {
        df <- df %>%
            dplyr::group_by_at(col) %>%
            dplyr::summarise(prob = sum(posterior), .groups = "drop") %>%
            dplyr::arrange_at(col) %>%
            dplyr::mutate(prob = prob / sum(prob))
        ggplot(df, aes(x = !!sym(col), y = prob)) +
            geom_area(alpha = 0.3) +
            geom_line(linewidth = 1.05) +
            labs(x = col, y = NULL, title = NULL) +
            theme_minimal(base_size = 30) +
            theme(
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white")
            )
    }
    subplots <- lapply(cols_to_plot, function(col) plot_marginal(grid, col))
    #---Output combined plot
    p <- gridExtra::arrangeGrob(
        grobs = subplots,
        ncol = length(cols_to_plot)
    )
    if (!is.null(filename)) {
        if (filetype == "png") {
            png(paste0(filename, ".png"), res = res, width = width, height = height, units = units)
        } else if (filetype == "jpeg" | filetype == "jpg") {
            jpeg(paste0(filename, ".jpg"), res = res, width = width, height = height, units = units, quality = 95)
        } else if (filetype == "svg") {
            svg(paste0(filename, ".svg"), width = width, height = height)
        } else if (filetype == "tiff" | filetype == "tif") {
            tiff(paste0(filename, ".tiff"), res = res, width = width, height = height, units = units)
        } else if (filetype == "eps") {
            setEPS()
            postscript(paste0(filename, ".eps"), width = width, height = height)
        } else if (filetype == "pdf") {
            pdf(paste0(filename, ".pdf"), width = width, height = height)
        }
        capture.output(
            {
                grid.draw(p)
                dev.off()
            },
            file = NULL
        )
        return()
    } else {
        return(p)
    }
}

plot_BESTish_joint <- function(grid,
                               filename = NULL,
                               filetype = "png",
                               width = 30, height = 15, units = "in", res = 150) {
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(viridis))
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(gridExtra))
    #---Find variables to plot
    tmp <- sapply(grid, function(x) length(unique(x)))
    cols_to_plot <- intersect(names(tmp[tmp > 1]), c("w0", "w1", "log10v0", "alpha", "R"))
    cat(bold(yellow("Plotting joint distributions from BESTish for variables:", paste0(cols_to_plot, collapse = ", "), "\n")))
    #---Create subplots for each joint distribution
    plot_joint <- function(df, xcol, ycol) {
        grid_step <- function(x) {
            ux <- sort(unique(x))
            if (length(ux) < 2) 0.01 else median(diff(ux))
        }
        df <- df %>%
            dplyr::group_by_at(c(xcol, ycol)) %>%
            dplyr::summarise(prob = sum(posterior), .groups = "drop")
        df$v1 <- df[[xcol]]
        df$v2 <- df[[ycol]]
        w1 <- grid_step(df$v1)
        w2 <- grid_step(df$v2)
        ggplot(df, aes(x = v1, y = v2)) +
            geom_tile(aes(fill = prob), width = w1, height = w2) +
            scale_fill_viridis(option = "C", name = "Posterior") +
            labs(x = xcol, y = ycol, title = NULL) +
            theme_minimal(base_size = 30) +
            theme(
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white")
            )
    }
    pairs <- combn(cols_to_plot, 2, simplify = FALSE)
    subplots <- lapply(pairs, function(pair) plot_joint(grid, pair[1], pair[2]))
    #---Output combined plot
    p <- gridExtra::arrangeGrob(
        grobs = subplots,
        ncol = length(pairs)
    )
    if (!is.null(filename)) {
        if (filetype == "png") {
            png(paste0(filename, ".png"), res = res, width = width, height = height, units = units)
        } else if (filetype == "jpeg" | filetype == "jpg") {
            jpeg(paste0(filename, ".jpg"), res = res, width = width, height = height, units = units, quality = 95)
        } else if (filetype == "svg") {
            svg(paste0(filename, ".svg"), width = width, height = height)
        } else if (filetype == "tiff" | filetype == "tif") {
            tiff(paste0(filename, ".tiff"), res = res, width = width, height = height, units = units)
        } else if (filetype == "eps") {
            setEPS()
            postscript(paste0(filename, ".eps"), width = width, height = height)
        } else if (filetype == "pdf") {
            pdf(paste0(filename, ".pdf"), width = width, height = height)
        }
        capture.output(
            {
                grid.draw(p)
                dev.off()
            },
            file = NULL
        )
        return()
    } else {
        return(p)
    }
}

plot_BESTish_vaf <- function(parameters,
                             data = NULL,
                             simulation_count = 0,
                             simulation_mean = FALSE,
                             time_max,
                             time_step = 0.01,
                             filename = NULL,
                             filetype = "png",
                             width = 30, height = 15, units = "in", res = 150) {
    suppressPackageStartupMessages(library(ggplot2))
    cat(bold(yellow("Plotting comparison of VAF trajectories between theory, simulations and data\n")))
    #---Compute theoretical mean/variance trajectories with selected parameters
    map_mean_var_df <- compute_expected_means_and_variances(
        time_grids = seq(0, time_max, by = time_step),
        h = time_step,
        w0 = parameters$w0,
        w1 = parameters$w1,
        v0 = parameters$v0,
        alpha = parameters$alpha
    )
    map_mean_var_df$vaf_var <- NA_real_
    for (i in 1:nrow(map_mean_var_df)) {
        ti <- map_mean_var_df$time[i]
        x <- c(map_mean_var_df$N0_bar[i], map_mean_var_df$N1_bar[i])
        Vt <- matrix(c(
            map_mean_var_df$V11[i], map_mean_var_df$V12[i],
            map_mean_var_df$V12[i], map_mean_var_df$V22[i]
        ), nrow = 2, byrow = TRUE)
        J <- Jg(x[1], x[2])
        map_mean_var_df$vaf_var[i] <- as.numeric(J %*% Vt %*% t(J)) / parameters$R
    }
    map_mean_var_df$vaf_lowerCI <- map_mean_var_df$vaf_mu - 1.96 * sqrt(map_mean_var_df$vaf_var)
    map_mean_var_df$vaf_upperCI <- map_mean_var_df$vaf_mu + 1.96 * sqrt(map_mean_var_df$vaf_var)
    #---Create stochastic simulations with selected parameters
    if (simulation_count > 0) {
        lambda0 <- parameters$w0 + 2 * parameters$v0
        lambda1 <- parameters$w1
        u0 <- parameters$v0 / lambda0
        u1 <- 0
        map_simulations <- run_replicates(
            n_reps = simulation_count,
            r = parameters$R,
            lambda_vec = c(lambda0, lambda1),
            u_vec = c(u0, u1),
            alpha = parameters$alpha,
            max_time = time_max,
            tau = time_step
        )
        map_simulations$vaf <- 0.5 * map_simulations$N1 / (map_simulations$N0 + map_simulations$N1)
        map_simulations_mean <- map_simulations %>%
            group_by(time) %>%
            summarise(across(where(is.numeric) & !matches("replicate"), mean, na.rm = TRUE))
    }
    #---Create plot for the VAF comparison
    title_text <- "VAF trajectories: theory (red)"
    p <- ggplot() +
        geom_ribbon(
            data = map_mean_var_df,
            aes(x = time, ymin = vaf_lowerCI, ymax = vaf_upperCI),
            fill = "#BC3C29", alpha = 0.3
        )
    if (simulation_count > 0) {
        title_text <- paste0(title_text, ", simulation (blue)")
        p <- p +
            geom_line(data = map_simulations, aes(x = time, y = vaf, group = replicate), color = "#0072B5", size = 1, alpha = 0.2)
    }
    p <- p +
        geom_line(data = map_mean_var_df, aes(x = time, y = vaf_mu), color = "#BC3C29", size = 5)
    if (simulation_count > 0 && simulation_mean) {
        p <- p +
            geom_line(data = map_simulations_mean, aes(x = time, y = vaf), color = "#0072B5", linetype = "dotted", size = 5)
    }
    if (!is.null(data)) {
        title_text <- paste0(title_text, ", data (green)")
        p <- p +
            geom_point(data = data, aes(x = Age, y = VAF), color = "white", size = 14) +
            geom_point(data = data, aes(x = Age, y = VAF), color = "#00A087", size = 10)
    }
    p <- p +
        labs(
            x = "Age",
            y = "Variant Allele Frequency (VAF)",
            title = title_text
        ) +
        theme_minimal(base_size = 50) +
        theme(
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "none"
        )
    #---Output combined plot
    if (!is.null(filename)) {
        if (filetype == "png") {
            png(paste0(filename, ".png"), res = res, width = width, height = height, units = units)
        } else if (filetype == "jpeg" | filetype == "jpg") {
            jpeg(paste0(filename, ".jpg"), res = res, width = width, height = height, units = units, quality = 95)
        } else if (filetype == "svg") {
            svg(paste0(filename, ".svg"), width = width, height = height)
        } else if (filetype == "tiff" | filetype == "tif") {
            tiff(paste0(filename, ".tiff"), res = res, width = width, height = height, units = units)
        } else if (filetype == "eps") {
            setEPS()
            postscript(paste0(filename, ".eps"), width = width, height = height)
        } else if (filetype == "pdf") {
            pdf(paste0(filename, ".pdf"), width = width, height = height)
        }
        capture.output(
            {
                print(p)
                dev.off()
            },
            file = NULL
        )
        return()
    } else {
        return(p)
    }
}
