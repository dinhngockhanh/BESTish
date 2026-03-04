#   Inference of mutation and selection rates from time-series and cohort mutant variant allele frequency datasets

##  Installation

The BESTish library (Bayesian Estimate for SelecTion Incorporating Scaling-limit to detect mutant Heterogeneity) [1] can be installed with

```{r}
devtools::install_github("dinhngockhanh/BESTish")
```

which usually takes under 1 minute.

### Software dependencies

BESTish requires the following R packages:
-   `dplyr` - data manipulation and transformation
-   `tidyr` - data tidying utilities
-   `readxl` - reading Excel files 

##  Vignettes

The `vignettes` folder contains examples of using BESTish to infer mutation rates and selection coefficients for individual mutations from longitudinal data in [2] and population-level data in [3,4], as tabulated in [5]:

-   File `data_processing.r` details our filtering process to determine potential drivers of clonal hematopoiesis (CH) from [2] and [5]. For each dataset, we further extract the variant allele frequency (VAF) values and corresponding time points (in years) for each potential CH driver. This file creates processed tables `fabre.csv`, `fabre_inference_list.csv`, `watson.csv` and `watson_inference_list.csv` for BESTish implementations.
-   File `longitudinal_inference.r` implements BESTish to compute the posterior distributions for the mutation rate and selection coefficient for each potential CH driver, given the corresponding (VAF, age) values from time-series data in [2] summarized in `fabre.csv` and `fabre_inference_list.csv`.
-   File `cohort_inference.r` implements BESTish to infer the mutation rates and selection coefficients for CH drivers from population-level data in [5] summarized in `watson.csv` and `watson_inference_list.csv`.

Runtime depends on the complexity (i.e., number of parameters) and granularity (i.e., number of bins per parameter) of the specific configuration, and can take up to a few hours per mutation.
Progress is visualized during the duration of inference; if the runtime is too long, we suggest reducing the grid sizes.

##  References
1.  Wang RY, Dinh KN, Taketomi K, Pang G, King KY, Kimmel M. BESTish: A diffusion-approximation framework for inferring selection and mutation in clonal hematopoiesis. bioRxiv (2026). https://doi.org/10.64898/2026.01.27.702030
2.  Fabre MA, de Almeida JG, Fiorillo E, Mitchell E, Damaskou A, Rak J, Orrù V, Marongiu M, Chapman MS, Vijayabaskar MS, Baxter J. The longitudinal dynamics and natural history of clonal haematopoiesis. Nature 606:335-42 (2022). https://doi.org/10.1038/s41586-022-04785-z
3.  McKerrell T, Park N, Moreno T, Grove CS, Ponstingl H, Stephens J, Crawley C, Craig J, Scott MA, Hodkinson C, Baxter J. Leukemia-associated somatic mutations drive distinct patterns of age-related clonal hemopoiesis. Cell Reports 10:1239-45 (2015). https://doi.org/10.1016/j.celrep.2015.02.005
4.  Coombs CC, Zehir A, Devlin SM, Kishtagari A, Syed A, Jonsson P, Hyman DM, Solit DB, Robson ME, Baselga J, Arcila ME. Therapy-related clonal hematopoiesis in patients with non-hematologic cancers is common and associated with adverse clinical outcomes. Cell Stem Cell 21:374-82 (2017). https://doi.org/10.1016/j.stem.2017.07.010
5.  Watson CJ, Papula AL, Poon GY, Wong WH, Young AL, Druley TE, Fisher DS, Blundell JR. The evolutionary dynamics and fitness landscape of clonal hematopoiesis. Science 367:1449-54 (2020). https://doi.org/10.1126/science.aay9333
