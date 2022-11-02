# L3AssessRPackage
Collection of age and length-based catch curve and per recruit analyses.

"L3Assess" contains a range of age and length-based catch curve and per recruit analyses, with various extensions, intended mainly for assessing stocks in data-limited situations. As these methods have strong (equilibrium) assumptions, they are likely most suited to data-limited situations when an assessment using a dynamic model is not possible, e.g. due to lack of a reliable abundance and/or catch indiex, but when (representative) age and/or length composition data to exist, together with biological information (i.e. growth, maturity, weight-length relationships etc.). The various methods in L3Assess (and similar approaches in other packages) may also have some value in complementing more sophisticated assessment analyses involving dynamic methods, but should not be used instead of these.
For more information on the various methods available in "L3Assess" and associated documentation, refer to help file and vignette. In some situations, there may be some benefit in using the R package "WAFishBioogy" to produce inputs (i.e. estimates of various biological parameters) required for anlaysis in "L3Assess". 

As with other packages in github, to install L3Assess, first ensure you have the 'devtools' package, otherwise install using install.packages("devtools').
Then, use: 
library(devtools)
devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE)
library(L3Assess).


