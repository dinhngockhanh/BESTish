suppressPackageStartupMessages({
    library(BESTish)
})
#-------------------------Create grid system for BESTish joint posterior
BESTish_grid <- make_BESTish_grid(
    w1_min = 1,
    w1_max = 1.3,
    w1_nbins = 500,
    log10v0_min = -6,
    log10v0_max = -3,
    log10v0_nbins = 500,
    alpha_min = 0.5,
    alpha_max = 1,
    alpha_nbins = 10,
    w0 = 1,
    R = 20000
)
#-------Infer posterior distributions of (w1,log10v0,alpha) with BESTish
#------------for each (study,mutation) pair in Watson et al. cohort data
watson <- read.csv("watson.csv")
watson$Age <- as.numeric(watson$Age)
watson$VAF <- as.numeric(watson$VAF)
watson_inference_list <- read.csv("watson_inference_list.csv")
for (row in 1:nrow(watson_inference_list)) {
    #---Retrieve (age,VAF) dataset for each (study,mutation) pair
    Study <- watson_inference_list$Study[row]
    Gene_ID <- watson_inference_list$Gene_ID[row]
    Protein_change_ID <- watson_inference_list$Protein_change_ID[row]
    outdir <- paste0("BESTish_cohort_", Gene_ID, "_", Protein_change_ID, "_", Study)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    watson_short <- watson[
        watson$Study == Study &
            watson$Gene_ID == Gene_ID &
            watson$Protein_change_ID == Protein_change_ID,
    ]
    data <- data.frame(
        Age = watson_short$Age,
        VAF = watson_short$VAF
    )
    #---Run BESTish to compute joint posterior for (w1,log10v0,alpha)
    BESTish_grid <- BESTish_inference(
        data = data,
        grid = BESTish_grid,
        mode = "cohort",
        time_step = 0.0005,
        parallel = TRUE
    )
    saveRDS(BESTish_grid, file = file.path(outdir, paste0("posterior_", Gene_ID, "_", Protein_change_ID, "_", Study, ".rds")))
    #---Plot marginal and joint distributions inferred from BESTish
    plot_BESTish_marginal(
        grid = BESTish_grid,
        filename = file.path(outdir, paste0("marginal_", Gene_ID, "_", Protein_change_ID, "_", Study)),
        filetype = "png"
    )
    plot_BESTish_joint(
        grid = BESTish_grid,
        filename = file.path(outdir, paste0("joint_", Gene_ID, "_", Protein_change_ID, "_", Study)),
        filetype = "png"
    )
    #---Plot comparison of VAF trajectories between theory and simulations
    #---with MAP estimates from BESTish and data
    plot_BESTish_vaf(
        parameters = BESTish_grid[which.max(BESTish_grid$posterior), ],
        data = data,
        simulation_count = 100,
        time_max = max(data$Age, na.rm = TRUE) + 1,
        filename = file.path(outdir, paste0("vaf_trajectory_", Gene_ID, "_", Protein_change_ID, "_", Study)),
        filetype = "png"
    )
}
