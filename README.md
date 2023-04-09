# mrMed
**mrMed** is a function to perform MR-based causal medaiton analysis based on summary-data Mendelian randomization.

Three default methods are provided to estimate the total effect (TE), direct effect (DE), indirect effect (IE), and mediation proportion (\rho):
1. Diff_IVW
2. Prod_IVW
3. Prod_Median


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

### Reference & Data source
```r
# 1. Pulit SL, Stoneman C, Morris AP, Wood AR, Glastonbury CA, 785 Tyrrell J, et al. Meta-analysis of genome-wide association 786 studies for body fat distribution in 694 649 individuals of 787 European ancestry. Hum Mol Genet. 2019;28(1):166-74. 788
# 2. Nikpay M, Goel A, Won HH, Hall LM, Willenborg C, Kanoni 789 S, et al. A comprehensive 1,000 Genomes-based genome-wide 790 association meta-analysis of coronary artery disease. Nat 791 Genet. 2015;47(10):1121-30. 792
# 3. Mahajan A, Taliun D, Thurner M, Robertson NR, Torres JM, 793 Rayner NW, et al. Fine-mapping type 2 diabetes loci to single-794 variant resolution using high-density imputation and islet-795 specific epigenome maps. Nat Genet. 2018;50(11):1505-13. 796
# 4. Wootton RE, Richmond RC, Stuijfzand BG, Lawn RB, Sallis 797 HM, Taylor GMJ, et al. Evidence for causal effects of lifetime 798 smoking on risk for depression and schizophrenia: a Mendelian 799 randomisation study. Psychol Med. 2020;50(14):2435-43.
# 5. Scott RA, Scott LJ, Mägi R, Marullo L, Gaulton KJ, Kaakinen M, et al. An expanded genome-wide association study of type 2 diabetes in Europeans. Diabetes. 2017;66(11):2888-902.
```