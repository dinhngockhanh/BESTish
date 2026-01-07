#   Inference of mutation and selection rates from time-series and cohort mutant VAF datasets

##  Installation

The BESTish library can be installed with

```{r}
devtools::install_github("dinhngockhanh/BESTish")
```

##  Vignettes

The `vignettes` folder contains examples of using BESTish to infer mutation rates and selection coefficients for individual mutations from population-level data in [1,2], as tabulated in [3], and longitudinal data in [4]:
-   sfdsdfs


For each dataset, the input consists of the list of VAF values and corresponding time points (in years).

##  References
1.  McKerrell T, Park N, Moreno T, Grove CS, Ponstingl H, Stephens J, Crawley C, Craig J, Scott MA, Hodkinson C, Baxter J. Leukemia-associated somatic mutations drive distinct patterns of age-related clonal hemopoiesis. Cell Reports 10:1239-45 (2015). https://doi.org/10.1016/j.celrep.2015.02.005
2.  Coombs CC, Zehir A, Devlin SM, Kishtagari A, Syed A, Jonsson P, Hyman DM, Solit DB, Robson ME, Baselga J, Arcila ME. Therapy-related clonal hematopoiesis in patients with non-hematologic cancers is common and associated with adverse clinical outcomes. Cell Stem Cell 21:374-82 (2017). https://doi.org/10.1016/j.stem.2017.07.010
3.  Watson CJ, Papula AL, Poon GY, Wong WH, Young AL, Druley TE, Fisher DS, Blundell JR. The evolutionary dynamics and fitness landscape of clonal hematopoiesis. Science 367:1449-54 (2020). https://doi.org/10.1126/science.aay9333
4.  Fabre MA, de Almeida JG, Fiorillo E, Mitchell E, Damaskou A, Rak J, Orrù V, Marongiu M, Chapman MS, Vijayabaskar MS, Baxter J. The longitudinal dynamics and natural history of clonal haematopoiesis. Nature 606:335-42 (2022). https://doi.org/10.1038/s41586-022-04785-z