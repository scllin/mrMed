# mrMed
devtools::install_github("scllin/mrMed")
library(mrMed)

#example 1
data(WHR_SMK_CAD)
mrMed(WHR_SMK_CAD)

#example 2
data(WHR_T2D_CAD)
mrMed(WHR_T2D_CAD)

#example 3
data(WHR_T2Dnoukb_CAD)
mrMed(WHR_T2Dnoukb_CAD)