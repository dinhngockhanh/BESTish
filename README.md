#   Inference of mutation and selection rates from time-series and cohort mutant variant allele frequency datasets

##  Installation

The [BESTish](https://doi.org/10.64898/2026.01.27.702030) (**B**ayesian **e**stimate for **s**elec**t**ion **i**ncorporating **s**caling-limit to detect mutant **h**eterogeneity) library can be installed with

```{r}
devtools::install_github("dinhngockhanh/BESTish")
```

##  Vignettes

The `vignettes` folder contains examples of using BESTish to infer mutation rates and selection coefficients for individual mutations from longitudinal data in [Fabre et al.](https://doi.org/10.1038/s41586-022-04785-z) and population-level data in [McKerrell et al.](https://doi.org/10.1016/j.celrep.2015.02.005) and [Coombs et al.](https://doi.org/10.1016/j.stem.2017.07.010), as tabulated in [Watson et al.](https://doi.org/10.1126/science.aay9333):

-   File `data_processing.r` details our filtering process to determine potential drivers of clonal hematopoiesis (CH) from [Fabre et al.](https://doi.org/10.1038/s41586-022-04785-z) and [Watson et al.](https://doi.org/10.1126/science.aay9333). For each dataset, we further extract the variant allele frequency (VAF) values and corresponding time points (in years) for each potential CH driver. This file creates processed tables `fabre.csv`, `fabre_inference_list.csv`, `watson.csv` and `watson_inference_list.csv` for BESTish implementations.
-   File `longitudinal_inference.r` implements BESTish to compute the posterior distributions for the mutation rate and selection coefficient for each potential CH driver, given the corresponding (VAF, age) values from time-series data in [Fabre et al.](https://doi.org/10.1038/s41586-022-04785-z) summarized in `fabre.csv` and `fabre_inference_list.csv`.
-   File `cohort_inference.r` implements BESTish to infer the mutation rates and selection coefficients for CH drivers from population-level data in [Watson et al.](https://doi.org/10.1126/science.aay9333) summarized in `watson.csv` and `watson_inference_list.csv`.

Runtime depends on the complexity (i.e., number of parameters) and granularity (i.e., number of bins per parameter) of the specific configuration, and can take up to a few hours per mutation.
Progress is visualized during the duration of inference; if the runtime is too long, we suggest reducing the grid sizes.

##  References
1.  Wang RY, Dinh KN, Taketomi K, Pang G, King KY, Kimmel M. [BESTish: A diffusion-approximation framework for inferring selection and mutation in clonal hematopoiesis](https://doi.org/10.64898/2026.01.27.702030). bioRxiv (2026).
2.  Fabre MA, de Almeida JG, Fiorillo E, Mitchell E, Damaskou A, Rak J, Orrù V, Marongiu M, Chapman MS, Vijayabaskar MS, Baxter J. [The longitudinal dynamics and natural history of clonal haematopoiesis](https://doi.org/10.1038/s41586-022-04785-z). Nature 606:335-42 (2022).
3.  McKerrell T, Park N, Moreno T, Grove CS, Ponstingl H, Stephens J, Crawley C, Craig J, Scott MA, Hodkinson C, Baxter J. [Leukemia-associated somatic mutations drive distinct patterns of age-related clonal hemopoiesis](https://doi.org/10.1016/j.celrep.2015.02.005). Cell Reports 10:1239-45 (2015).
4.  Coombs CC, Zehir A, Devlin SM, Kishtagari A, Syed A, Jonsson P, Hyman DM, Solit DB, Robson ME, Baselga J, Arcila ME. [Therapy-related clonal hematopoiesis in patients with non-hematologic cancers is common and associated with adverse clinical outcomes](https://doi.org/10.1016/j.stem.2017.07.010). Cell Stem Cell 21:374-82 (2017).
5.  Watson CJ, Papula AL, Poon GY, Wong WH, Young AL, Druley TE, Fisher DS, Blundell JR. [The evolutionary dynamics and fitness landscape of clonal hematopoiesis](https://doi.org/10.1126/science.aay9333). Science 367:1449-54 (2020).
