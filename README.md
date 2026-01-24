#   Inference of mutation and selection rates from time-series and cohort mutant variant allele frequency datasets

##  Installation

This is <u>underlined text</u>.

The BESTish library (**B**ayesian **e**stimate for **s**elec**t**ion **i**ncorporating **s**caling-limit to detect mutant **h**eterogeneity) can be installed with

```{r}
devtools::install_github("dinhngockhanh/BESTish")
```

##  Vignettes

The `vignettes` folder contains examples of using BESTish to infer mutation rates and selection coefficients for individual mutations from longitudinal data in [1] and population-level data in [2,3], as tabulated in [4]:
-   File `data_processing.r` details our filtering process to determine potential drivers of clonal hematopoiesis (CH) from [1] and [4]. For each dataset, we further extract the variant allele frequency (VAF) values and corresponding time points (in years) for each potential CH driver. This file creates processed tables `fabre.csv`, `fabre_inference_list.csv`, `watson.csv` and `watson_inference_list.csv` for BESTish implementations.
-   File `longitudinal_inference.r` implements BESTish to compute the posterior distributions for the mutation rate and selection coefficient for each potential CH driver, given the corresponding (VAF, age) values from time-series data in [1] summarized in `fabre.csv` and `fabre_inference_list.csv`.
-   File `cohort_inference.r` implements BESTish to infer the mutation rates and selection coefficients for CH drivers from population-level data in [4] summarized in `watson.csv` and `watson_inference_list.csv`.

##  References
1.  Fabre MA, de Almeida JG, Fiorillo E, Mitchell E, Damaskou A, Rak J, Orrù V, Marongiu M, Chapman MS, Vijayabaskar MS, Baxter J. The longitudinal dynamics and natural history of clonal haematopoiesis. Nature 606:335-42 (2022). https://doi.org/10.1038/s41586-022-04785-z
2.  McKerrell T, Park N, Moreno T, Grove CS, Ponstingl H, Stephens J, Crawley C, Craig J, Scott MA, Hodkinson C, Baxter J. Leukemia-associated somatic mutations drive distinct patterns of age-related clonal hemopoiesis. Cell Reports 10:1239-45 (2015). https://doi.org/10.1016/j.celrep.2015.02.005
3.  Coombs CC, Zehir A, Devlin SM, Kishtagari A, Syed A, Jonsson P, Hyman DM, Solit DB, Robson ME, Baselga J, Arcila ME. Therapy-related clonal hematopoiesis in patients with non-hematologic cancers is common and associated with adverse clinical outcomes. Cell Stem Cell 21:374-82 (2017). https://doi.org/10.1016/j.stem.2017.07.010
4.  Watson CJ, Papula AL, Poon GY, Wong WH, Young AL, Druley TE, Fisher DS, Blundell JR. The evolutionary dynamics and fitness landscape of clonal hematopoiesis. Science 367:1449-54 (2020). https://doi.org/10.1126/science.aay9333
