#   Inference of mutation and selection rates from time-series and cohort mutant VAF datasets

##  Installation

The BESTish library can be installed with

```{r}
devtools::install_github("dinhngockhanh/BESTish")
```

##  Vignettes

The `vignettes` folder contains examples of using BESTish to infer mutation rates and selection coefficients for specific mutations from cohort data in [1] 


ABC-SMC-(D)RF, ABC-RF, ABC-DRF and traditional ABC methods to infer parameters in different mathematical models. These examples are showcased in the paper "Approximate Bayesian Computation sequential Monte Carlo via random forests" [1].

##  References
1.  Coombs CC, Zehir A, Devlin SM, Kishtagari A, Syed A, Jonsson P, Hyman DM, Solit DB, Robson ME, Baselga J, Arcila ME. Therapy-related clonal hematopoiesis in patients with non-hematologic cancers is common and associated with adverse clinical outcomes. Cell Stem Cell 21(3):374-82 (2017). 
2.  McKerrell T, Park N, Moreno T, Grove CS, Ponstingl H, Stephens J, Crawley C, Craig J, Scott MA, Hodkinson C, Baxter J. Leukemia-associated somatic mutations drive distinct patterns of age-related clonal hemopoiesis. Cell Reports 10(8):1239-45 (2015).
3.  Fabre MA, de Almeida JG, Fiorillo E, Mitchell E, Damaskou A, Rak J, Orrù V, Marongiu M, Chapman MS, Vijayabaskar MS, Baxter J. The longitudinal dynamics and natural history of clonal haematopoiesis. Nature 606(7913):335-42 (2022).