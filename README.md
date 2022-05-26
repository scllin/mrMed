# mrMed
mrMed (**MR-based mediation analysis**) is a function to perform causal medaiton analysis using summary satistics from genome-wide association studies (GWAS) based on Mendelian Randomization (MR).

mrMed provides three default methods to perform the MR-based mediation analysis and estimates of total effect (TE), direct effect (DE), indirect effect (IE), and mediation proportion (rho):
1. Diff_IVW
2. Prod_IVW
3. Prod_Median

### Reference
Causal Mediation Analysis: A Summary-Data Mendelian Randomization Approach
<https://...>

### 1. Install and load mrMed
To install and load the latest version from GitHub, run the following lines:
```r
if (!require("devtools")) { install.packages("devtools") } else {}
devtools::install_github("scllin/mrMed")
```

```r
library(mrMed)
```

### 2. Examples
```r
# Load the dataset with exposure:WHR, mediator:T2D, and outcome:CAD
data(WHR_T2D_CAD)

# Run mrMed with the default methods: Diff_IVW, Prod_IVW, Prod_Median
mrMed(WHR_T2D_CAD)

# Run mrMed with the the other methods
mrMed(WHR_T2D_CAD, method_list=c("Diff_IVW_0","Prod_IVW_0"))

# Run Diff-Median method (require certain run time due to the bootstrap procedure)
mrMed(WHR_T2D_CAD, method_list=c("Diff_Median"))
```