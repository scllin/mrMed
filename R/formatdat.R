
#add columns to fit required format for TwoSampleMR
form_dat <- function(dat_mrMed){
	if(!"id.X"%in%colnames(dat_mrMed)){id.X <- rep("X",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.X)}
	if(!"X"%in%colnames(dat_mrMed)){X <- rep("X",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,X)}
	if(!"id.M"%in%colnames(dat_mrMed)){id.M <- rep("M",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.M)}
	if(!"M"%in%colnames(dat_mrMed)){M <- rep("M",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,M)}
	if(!"id.Y"%in%colnames(dat_mrMed)){id.Y <- rep("Y",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,id.Y)}
	if(!"Y"%in%colnames(dat_mrMed)){Y <- rep("Y",dim(dat_mrMed)[1]);dat_mrMed<-cbind(dat_mrMed,Y)}
	return(dat_mrMed)
}
