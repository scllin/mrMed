library(mrMed)

setwd("/staging/biology/kk037356/res/mrMed/4_mrMed_realdata")

#output:
#Table3 
#Table4
#TableS1
#TableS5
#TableS2-S4

#===Table3: MR-based mediation X:WHR, M:SMK, Y:CAD
set.seed(100)
data(WHR_SMK_CAD)
res1 <- mrMed(WHR_SMK_CAD,method_list=c("Diff_IVW_0","Diff_IVW","Diff_Egger","Diff_Median","Prod_IVW_0","Prod_IVW","Prod_Egger","Prod_Median"))
tab3 <- array(0,c(8,4))
temp <- apply(cbind(round(exp(res1$TE[,c("CI_lower")]),2),round(exp(res1$TE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab3[,1] <- paste0(round(exp(unlist(res1$TE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res1$DE[,c("CI_lower")]),2),round(exp(res1$DE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab3[,2] <- paste0(round(exp(unlist(res1$DE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res1$IE[,c("CI_lower")]),2),round(exp(res1$IE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab3[,3] <- paste0(round(exp(unlist(res1$IE[,c("b")])),2),temp)
temp <- apply(cbind(round(res1$rho[,c("CI_lower")],3),round(res1$rho[,c("CI_upper")],3)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab3[,4] <- paste0(round(res1$rho[,c("b")],3),temp)
write.table(tab3,file="Table3.csv",sep=",",row.names=FALSE,col.names=FALSE)

#===Table4: MR-based mediation X:WHR, M:T2D, Y:CAD
set.seed(100)
data(WHR_T2D_CAD)
res2 <- mrMed(WHR_T2D_CAD,method_list=c("Diff_IVW_0","Diff_IVW","Diff_Egger","Diff_Median","Prod_IVW_0","Prod_IVW","Prod_Egger","Prod_Median"))
tab4 <- array(0,c(8,4))
temp <- apply(cbind(round(exp(res2$TE[,c("CI_lower")]),2),round(exp(res2$TE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab4[,1] <- paste0(round(exp(unlist(res2$TE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res2$DE[,c("CI_lower")]),2),round(exp(res2$DE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab4[,2] <- paste0(round(exp(unlist(res2$DE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res2$IE[,c("CI_lower")]),2),round(exp(res2$IE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab4[,3] <- paste0(round(exp(unlist(res2$IE[,c("b")])),2),temp)
temp <- apply(cbind(round(res2$rho[,c("CI_lower")],3),round(res2$rho[,c("CI_upper")],3)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tab4[,4] <- paste0(round(res2$rho[,c("b")],3),temp)
write.table(tab4,file="Table4.csv",sep=",",row.names=FALSE,col.names=FALSE)

#===TableS1 MR-based mediation X:WHR, M:T2D(noukb), Y:CAD
set.seed(100)
data(WHR_T2Dnoukb_CAD)
res3 <- mrMed(WHR_T2Dnoukb_CAD,method_list=c("Diff_IVW_0","Diff_IVW","Diff_Egger","Diff_Median","Prod_IVW_0","Prod_IVW","Prod_Egger","Prod_Median"))
tabS1 <- array(0,c(8,4))
temp <- apply(cbind(round(exp(res3$TE[,c("CI_lower")]),2),round(exp(res3$TE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tabS1[,1] <- paste0(round(exp(unlist(res3$TE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res3$DE[,c("CI_lower")]),2),round(exp(res3$DE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tabS1[,2] <- paste0(round(exp(unlist(res3$DE[,c("b")])),2),temp)
temp <- apply(cbind(round(exp(res3$IE[,c("CI_lower")]),2),round(exp(res3$IE[,c("CI_upper")]),2)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tabS1[,3] <- paste0(round(exp(unlist(res3$IE[,c("b")])),2),temp)
temp <- apply(cbind(round(res3$rho[,c("CI_lower")],3),round(res3$rho[,c("CI_upper")],3)),1,paste,collapse=",")
temp <- apply(cbind(rep("(",8),temp),1,paste,collapse="")
temp <- apply(cbind(temp,rep(")",8)),1,paste,collapse="")
tabS1[,4] <- paste0(round(res3$rho[,c("b")],3),temp)
write.table(tabS1,file="TableS1.csv",sep=",",row.names=FALSE,col.names=FALSE)



#====Table S5: F-stat
library(MVMR)
Fstat <- array(0,c(6,2))

Fstat[1,1] <- mean((WHR_SMK_CAD$beta.X[WHR_SMK_CAD$Gx_plum==1]/WHR_SMK_CAD$se.X[WHR_SMK_CAD$Gx_plum==1])^2)
Fstat[2,1] <- mean((WHR_SMK_CAD$beta.M[WHR_SMK_CAD$Gm_plum==1]/WHR_SMK_CAD$se.M[WHR_SMK_CAD$Gm_plum==1])^2)
indx <- which(WHR_SMK_CAD$G_mvmr==1)
r_input <- format_mvmr(
	BXGs = WHR_SMK_CAD[indx,c("beta.X","beta.M")],
	BYG = WHR_SMK_CAD$beta.Y[indx],
	seBXGs = WHR_SMK_CAD[indx,c("se.X","se.M")],
	seBYG = WHR_SMK_CAD$se.Y[indx],RSID = WHR_SMK_CAD$SNP[indx])
Fstat[1:2,2] <- as.numeric(strength_mvmr(r_input, 0))
colnames(Fstat) <- c("meanF","condi_F")

Fstat[3,1] <- mean((WHR_T2D_CAD$beta.X[WHR_T2D_CAD$Gx_plum==1]/WHR_T2D_CAD$se.X[WHR_T2D_CAD$Gx_plum==1])^2)
Fstat[4,1] <- mean((WHR_T2D_CAD$beta.M[WHR_T2D_CAD$Gm_plum==1]/WHR_T2D_CAD$se.M[WHR_T2D_CAD$Gm_plum==1])^2)
indx <- which(WHR_T2D_CAD$G_mvmr==1)
r_input <- format_mvmr(
	BXGs = WHR_T2D_CAD[indx,c("beta.X","beta.M")],
	BYG = WHR_T2D_CAD$beta.Y[indx],
	seBXGs = WHR_T2D_CAD[indx,c("se.X","se.M")],
	seBYG = WHR_T2D_CAD$se.Y[indx],RSID = WHR_T2D_CAD$SNP[indx])
Fstat[3:4,2] <- as.numeric(strength_mvmr(r_input, 0))
colnames(Fstat) <- c("meanF","condi_F")

Fstat[5,1] <- mean((WHR_T2Dnoukb_CAD$beta.X[WHR_T2Dnoukb_CAD$Gx_plum==1]/WHR_T2Dnoukb_CAD$se.X[WHR_T2Dnoukb_CAD$Gx_plum==1])^2)
Fstat[6,1] <- mean((WHR_T2Dnoukb_CAD$beta.M[WHR_T2Dnoukb_CAD$Gm_plum==1]/WHR_T2Dnoukb_CAD$se.M[WHR_T2Dnoukb_CAD$Gm_plum==1])^2)
indx <- which(WHR_T2Dnoukb_CAD$G_mvmr==1)
r_input <- format_mvmr(
	BXGs = WHR_T2Dnoukb_CAD[indx,c("beta.X","beta.M")],
	BYG = WHR_T2Dnoukb_CAD$beta.Y[indx],
	seBXGs = WHR_T2Dnoukb_CAD[indx,c("se.X","se.M")],
	seBYG = WHR_T2Dnoukb_CAD$se.Y[indx],RSID = WHR_T2Dnoukb_CAD$SNP[indx])
Fstat[5:6,2] <- as.numeric(strength_mvmr(r_input, 0))
colnames(Fstat) <- c("meanF","condi_F")

round(Fstat,1)
write.table(Fstat,file="TableS5.csv",sep=",",row.names=FALSE,col.names=FALSE)



#===TableS2 
library(TwoSampleMR)
set.seed(1234567)
gamma=0.05
z <- abs(qnorm(gamma/2))
tabS2 <- array(0,c(8,3))

#mr2: X->Y using Gx
dat_temp <- WHR_SMK_CAD[!is.na(WHR_SMK_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_XtoY <- harmonise_data(datforX,datforY)
mrXtoY_Gx <- mr(dat_XtoY[dat_XtoY$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[1,1] <- mrXtoY_Gx$nsnp[1] 
tabS2[2:4,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx$b,2),rep("(",3),round(mrXtoY_Gx$b-z*mrXtoY_Gx$se,2),",",round(mrXtoY_Gx$b+z*mrXtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->Y using Gx_plum 
mrXtoY_Gx_plum <- mr(dat_XtoY[dat_XtoY$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[5,1] <- mrXtoY_Gx_plum$nsnp[1] 
tabS2[6:8,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx_plum$b,2),rep("(",3),round(mrXtoY_Gx_plum$b-z*mrXtoY_Gx_plum$se,2),",",round(mrXtoY_Gx_plum$b+z*mrXtoY_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx
dat_temp <- WHR_SMK_CAD[!is.na(WHR_SMK_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.M),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]
datforM$M<-"M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","outcome",colnames(datforM))
dat_XtoM <- harmonise_data(datforX,datforM)
mrXtoM_Gx <- mr(dat_XtoM[dat_XtoM$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[1,2] <- mrXtoM_Gx$nsnp[1] 
tabS2[2:4,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx$b,2),rep("(",3),round(mrXtoM_Gx$b-z*mrXtoM_Gx$se,2),",",round(mrXtoM_Gx$b+z*mrXtoM_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx_plum 
mrXtoM_Gx_plum <- mr(dat_XtoM[dat_XtoM$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[5,2] <- mrXtoM_Gx_plum$nsnp[1] 
tabS2[6:8,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx_plum$b,2),rep("(",3),round(mrXtoM_Gx_plum$b-z*mrXtoM_Gx_plum$se,2),",",round(mrXtoM_Gx_plum$b+z*mrXtoM_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm
dat_temp <- WHR_SMK_CAD[!is.na(WHR_SMK_CAD$beta.M),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M","Gx","Gx_plum","Gm","Gm_plum")]
datforM$M <- "M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","exposure",colnames(datforM))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_MtoY <- harmonise_data(datforM,datforY)
mrMtoY_Gx <- mr(dat_MtoY[dat_MtoY$Gm==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[1,3] <- mrMtoY_Gx$nsnp[1] 
tabS2[2:4,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gx$b,2),rep("(",3),round(mrMtoY_Gx$b-z*mrMtoY_Gx$se,2),",",round(mrMtoY_Gx$b+z*mrMtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm_plum 
mrMtoY_Gm_plum <- mr(dat_MtoY[dat_MtoY$Gm_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS2[5,3] <- mrMtoY_Gm_plum$nsnp[1] 
tabS2[6:8,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gm_plum$b,2),rep("(",3),round(mrMtoY_Gm_plum$b-z*mrMtoY_Gm_plum$se,2),",",round(mrMtoY_Gm_plum$b+z*mrMtoY_Gm_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

write.table(tabS2,file="TableS2.csv",sep=",",row.names=FALSE,col.names=FALSE)



#===TableS3 
set.seed(1234567)
gamma=0.05
z <- abs(qnorm(gamma/2))
tabS3 <- array(0,c(8,3))

#mr2: X->Y using Gx
dat_temp <- WHR_T2D_CAD[!is.na(WHR_T2D_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_XtoY <- harmonise_data(datforX,datforY)
mrXtoY_Gx <- mr(dat_XtoY[dat_XtoY$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[1,1] <- mrXtoY_Gx$nsnp[1] 
tabS3[2:4,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx$b,2),rep("(",3),round(mrXtoY_Gx$b-z*mrXtoY_Gx$se,2),",",round(mrXtoY_Gx$b+z*mrXtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->Y using Gx_plum 
mrXtoY_Gx_plum <- mr(dat_XtoY[dat_XtoY$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[5,1] <- mrXtoY_Gx_plum$nsnp[1] 
tabS3[6:8,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx_plum$b,2),rep("(",3),round(mrXtoY_Gx_plum$b-z*mrXtoY_Gx_plum$se,2),",",round(mrXtoY_Gx_plum$b+z*mrXtoY_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx
dat_temp <- WHR_T2D_CAD[!is.na(WHR_T2D_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.M),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]
datforM$M<-"M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","outcome",colnames(datforM))
dat_XtoM <- harmonise_data(datforX,datforM)
mrXtoM_Gx <- mr(dat_XtoM[dat_XtoM$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[1,2] <- mrXtoM_Gx$nsnp[1] 
tabS3[2:4,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx$b,2),rep("(",3),round(mrXtoM_Gx$b-z*mrXtoM_Gx$se,2),",",round(mrXtoM_Gx$b+z*mrXtoM_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx_plum 
mrXtoM_Gx_plum <- mr(dat_XtoM[dat_XtoM$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[5,2] <- mrXtoM_Gx_plum$nsnp[1] 
tabS3[6:8,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx_plum$b,2),rep("(",3),round(mrXtoM_Gx_plum$b-z*mrXtoM_Gx_plum$se,2),",",round(mrXtoM_Gx_plum$b+z*mrXtoM_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm
dat_temp <- WHR_T2D_CAD[!is.na(WHR_T2D_CAD$beta.M),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M","Gx","Gx_plum","Gm","Gm_plum")]
datforM$M <- "M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","exposure",colnames(datforM))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_MtoY <- harmonise_data(datforM,datforY)
mrMtoY_Gx <- mr(dat_MtoY[dat_MtoY$Gm==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[1,3] <- mrMtoY_Gx$nsnp[1] 
tabS3[2:4,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gx$b,2),rep("(",3),round(mrMtoY_Gx$b-z*mrMtoY_Gx$se,2),",",round(mrMtoY_Gx$b+z*mrMtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm_plum 
mrMtoY_Gm_plum <- mr(dat_MtoY[dat_MtoY$Gm_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS3[5,3] <- mrMtoY_Gm_plum$nsnp[1] 
tabS3[6:8,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gm_plum$b,2),rep("(",3),round(mrMtoY_Gm_plum$b-z*mrMtoY_Gm_plum$se,2),",",round(mrMtoY_Gm_plum$b+z*mrMtoY_Gm_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

write.table(tabS3,file="TableS3.csv",sep=",",row.names=FALSE,col.names=FALSE)




#===TableS4 
set.seed(1234567)
gamma=0.05
z <- abs(qnorm(gamma/2))
tabS4 <- array(0,c(8,3))

#mr2: X->Y using Gx
dat_temp <- WHR_T2Dnoukb_CAD[!is.na(WHR_T2Dnoukb_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_XtoY <- harmonise_data(datforX,datforY)
mrXtoY_Gx <- mr(dat_XtoY[dat_XtoY$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[1,1] <- mrXtoY_Gx$nsnp[1] 
tabS4[2:4,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx$b,2),rep("(",3),round(mrXtoY_Gx$b-z*mrXtoY_Gx$se,2),",",round(mrXtoY_Gx$b+z*mrXtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->Y using Gx_plum 
mrXtoY_Gx_plum <- mr(dat_XtoY[dat_XtoY$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[5,1] <- mrXtoY_Gx_plum$nsnp[1] 
tabS4[6:8,1] <- apply(cbind(apply(cbind(round(mrXtoY_Gx_plum$b,2),rep("(",3),round(mrXtoY_Gx_plum$b-z*mrXtoY_Gx_plum$se,2),",",round(mrXtoY_Gx_plum$b+z*mrXtoY_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx
dat_temp <- WHR_T2Dnoukb_CAD[!is.na(WHR_T2Dnoukb_CAD$beta.X),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.M),]
datforX <- dat_temp[,c("SNP","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X","Gx","Gx_plum","Gm","Gm_plum")]
datforX$X <- "X"; datforX$id.X <- "X" 
colnames(datforX) <- gsub("X","exposure",colnames(datforX))
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]
datforM$M<-"M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","outcome",colnames(datforM))
dat_XtoM <- harmonise_data(datforX,datforM)
mrXtoM_Gx <- mr(dat_XtoM[dat_XtoM$Gx==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[1,2] <- mrXtoM_Gx$nsnp[1] 
tabS4[2:4,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx$b,2),rep("(",3),round(mrXtoM_Gx$b-z*mrXtoM_Gx$se,2),",",round(mrXtoM_Gx$b+z*mrXtoM_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: X->M using Gx_plum 
mrXtoM_Gx_plum <- mr(dat_XtoM[dat_XtoM$Gx_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[5,2] <- mrXtoM_Gx_plum$nsnp[1] 
tabS4[6:8,2] <- apply(cbind(apply(cbind(round(mrXtoM_Gx_plum$b,2),rep("(",3),round(mrXtoM_Gx_plum$b-z*mrXtoM_Gx_plum$se,2),",",round(mrXtoM_Gx_plum$b+z*mrXtoM_Gx_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm
dat_temp <- WHR_T2Dnoukb_CAD[!is.na(WHR_T2Dnoukb_CAD$beta.M),]
dat_temp <- dat_temp[!is.na(dat_temp$beta.Y),]
datforM <- dat_temp[,c("SNP","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M","Gx","Gx_plum","Gm","Gm_plum")]
datforM$M <- "M"; datforM$id.M <- "M" 
colnames(datforM) <- gsub("M","exposure",colnames(datforM))
datforY <- dat_temp[,c("SNP","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
datforY$Y <- "Y"; datforY$id.Y <- "Y" 
colnames(datforY) <- gsub("Y","outcome",colnames(datforY))
dat_MtoY <- harmonise_data(datforM,datforY)
mrMtoY_Gx <- mr(dat_MtoY[dat_MtoY$Gm==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[1,3] <- mrMtoY_Gx$nsnp[1] 
tabS4[2:4,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gx$b,2),rep("(",3),round(mrMtoY_Gx$b-z*mrMtoY_Gx$se,2),",",round(mrMtoY_Gx$b+z*mrMtoY_Gx$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")

#mr2: M->Y using Gm_plum 
mrMtoY_Gm_plum <- mr(dat_MtoY[dat_MtoY$Gm_plum==1,],method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))	
tabS4[5,3] <- mrMtoY_Gm_plum$nsnp[1] 
tabS4[6:8,3] <- apply(cbind(apply(cbind(round(mrMtoY_Gm_plum$b,2),rep("(",3),round(mrMtoY_Gm_plum$b-z*mrMtoY_Gm_plum$se,2),",",round(mrMtoY_Gm_plum$b+z*mrMtoY_Gm_plum$se,2)),1,paste,collapse=""),rep(")",3)),1,paste,collapse="")


write.table(tabS4,file="TableS4.csv",sep=",",row.names=FALSE,col.names=FALSE)
