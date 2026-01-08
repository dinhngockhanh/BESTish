suppressPackageStartupMessages({
    devtools::load_all() ############################################### REPLACE WITH library(BESTish)
})
#-------------------------Create grid system for BESTish joint posterior
BESTish_grid <- make_BESTish_grid(
    w1_min = 1,
    w1_max = 1.3,
    w1_nbins = 5, ###################################################### REPLACE WITH w1_nbins = 500
    log10v0_min = -6,
    log10v0_max = -3,
    log10v0_nbins = 5, ################################################# REPLACE WITH log10v0_nbins = 500
    alpha = 1,
    w0 = 1,
    R = 20000
)
#-------------Infer posterior distributions of (w1,log10v0) with BESTish
#-----for each (patient,mutation) pair in Fabre et al. longitudinal data
fabre <- read.csv("fabre.csv")
fabre$Age <- as.numeric(fabre$Age)
fabre$VAF <- as.numeric(fabre$VAF)
fabre_inference_list <- read.csv("fabre_inference_list.csv")
for (row in 1:nrow(fabre_inference_list)) {
    #---Retrieve (age,VAF) dataset for each (patient,mutation) pair
    Sample_ID <- fabre_inference_list$Sample_ID[row]
    Gene_ID <- fabre_inference_list$Gene_ID[row]
    Protein_change_ID <- fabre_inference_list$Protein_change_ID[row]
    outdir <- paste0("BESTish_longitudinal_", Gene_ID, "_", Protein_change_ID, "_", Sample_ID)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    fabre_short <- fabre[
        fabre$Sample_ID == Sample_ID &
            fabre$Gene_ID == Gene_ID &
            fabre$Protein_change_ID == Protein_change_ID,
    ]
    data <- data.frame(
        Age = fabre_short$Age,
        VAF = fabre_short$VAF
    )
    #---Run BESTish to compute joint posterior for (w1,log10v0)
    BESTish_grid <- BESTish_inference(
        data = data,
        grid = BESTish_grid,
        mode = "longitudinal",
        time_step = 0.001, ############################################# REPLACE WITH time_step = 0.0005
        parallel = TRUE
    )
    saveRDS(BESTish_grid, file = file.path(outdir, paste0("posterior_", Gene_ID, "_", Protein_change_ID, "_", Sample_ID, ".rds")))
    #---Plot marginal and joint distributions inferred from BESTish
    plot_BESTish_marginal(
        grid = BESTish_grid,
        filename = file.path(outdir, paste0("marginal_", Gene_ID, "_", Protein_change_ID, "_", Sample_ID)),
        filetype = "png"
    )
    plot_BESTish_joint(
        grid = BESTish_grid,
        filename = file.path(outdir, paste0("joint_", Gene_ID, "_", Protein_change_ID, "_", Sample_ID)),
        filetype = "png"
    )
    #---Plot comparison of VAF trajectories between theory and simulations
    #---with MAP estimates from BESTish and data
    plot_BESTish_vaf(
        parameters = BESTish_grid[which.max(BESTish_grid$posterior), ],
        data = data,
        simulation_count = 3, ########################################## REPLACE WITH simulation_count = 100
        time_max = max(data$Age, na.rm = TRUE) + 1,
        filename = file.path(outdir, paste0("vaf_trajectory_", Gene_ID, "_", Protein_change_ID, "_", Sample_ID)),
        filetype = "png"
    )
}
