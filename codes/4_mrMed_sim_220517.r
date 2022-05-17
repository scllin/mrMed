library(mrMed)
library(optparse)
library(mvtnorm)
library(MVMR)

#================================================#
#==========functions used in simulation==========#
#================================================#

reproducible_runif <- function(N,a,b,seed = NULL){    
	set.seed(seed)
	return(runif(N,a,b))
    
}

reproducible_sample <- function(G,size,seed = NULL){    
	set.seed(seed)
	return(sample(G,size))
}

reproducible_rnorm <- function(N,a,b,seed = NULL){    
	set.seed(seed)
	return(rnorm(N,a,b))
}

rcorrbinom <- function(n, prob1, prob2,corr){

	P00   <- (1-prob1)*(1-prob2) + corr*sqrt(prob1*prob2*(1-prob1)*(1-prob2))
	P01   <- 1 - prob1 - P00
	P10   <- 1 - prob2 - P00
	P11   <- P00 + prob1 + prob2 - 1
	PROBS <- c(P00, P01, P10, P11)
	if (min(PROBS) < 0)       { stop('Error: corr is not in the allowable range') }

	tri <- list(c(0,0),c(0,1),c(1,0),c(1,1))
	temp_t1 <- tri[sample(x=1:length(tri),size=n,prob=c(P00,P01,P10,P11),replace=TRUE)]
	t1 <- t(matrix(unlist(temp_t1),2,n))
	temp_t2 <- tri[sample(x=1:length(tri),size=n,prob=c(P00,P01,P10,P11),replace=TRUE)]
	t2 <- t(matrix(unlist(temp_t2),2,n))
	t_12 <- t1+t2

	return(t_12)
}

myGWAS <- function(XMY,Gk_XMY){

	XMY_dmean <- XMY - mean(XMY)
	Gk_XMY_dmean <- sapply(1:dim(Gk_XMY)[2], function(j){Gk_XMY[,j]-mean(Gk_XMY[,j])})


	beta_XMY <- cov(Gk_XMY_dmean,XMY_dmean)/apply(Gk_XMY_dmean,2,var)

	se_XMY <- sapply(1:dim(Gk_XMY)[2],function(j){sqrt(sum((XMY_dmean-beta_XMY[j]*Gk_XMY_dmean[,j])^2)/(dim(Gk_XMY)[1]-2))})/sqrt(apply(Gk_XMY_dmean^2,2,sum))

	pval_XMY <- 2*pt(-abs(beta_XMY/se_XMY),dim(Gk_XMY)[1]-2)

	return(data.frame(beta_XMY=beta_XMY,se_XMY=se_XMY,pval_XMY=pval_XMY))
	}


#=====================================#
#==========simulation set-up==========#
#=====================================#

#===optional parameters
option_list <- list(
  make_option("--R", type="integer", default=1000, help="number of simulation runs"),
  make_option("--N", type="integer", default=80000, help="sample size"),
  make_option("--alpha", type="numeric", default=0.3, help="effect of X on M"),
  make_option("--beta", type="numeric", default=0.3, help="effect of M on Y"),
  make_option("--delta", type="numeric", default=0.21, help="effect of X on Y after controlling M"),
  make_option("--pleiotropy", type="integer", default=1, help="S1-S11")
)
parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

#===number of SNPs
N1_Gx <- 60; N2_Gx <- 20; N1_Gm <- 60; N2_Gm <- 20 #here we set symmetric setting for Gx Gm, should be revised with carefully examination
N_Gx <- N1_Gx + N2_Gx
N_Gm <- N1_Gm + N2_Gm
N_Gk <- N_Gx + N_Gm

#===fixed parameters
seedC <- 1 #changed for another parameter setting 

a_xk <-  c(reproducible_runif(N_Gx,0.1,0.3,1000*seedC)) 
a_mk <-  c(reproducible_runif(N_Gm,0.1,0.3,3000*seedC))
b_xk <-  rep(0,N_Gx)
b_mk <-  rep(0,N_Gm)
c_xk <-  rep(0,N_Gx)
c_mk <-  rep(0,N_Gm)
d_x1k <-  rep(0,N_Gx)
d_x2k <-  rep(0,N_Gx)
d_x3k <-  rep(0,N_Gx)
d_m1k <-  rep(0,N_Gm)
d_m2k <-  rep(0,N_Gm)
d_m3k <-  rep(0,N_Gm)

omega <- 5 #random error std
Sig <- matrix(c(3,1,2,1,3,1,2,1,3), ncol=3) #cov matrix for U1, U2, U3 (vars affected by confounders)

#===selected pleiotropy SNPs
leng_1 <- N1_Gx/5; leng_2 <- N2_Gx/5
leng_3 <- N1_Gm/5; leng_4 <- N2_Gm/5
S_xb <- c(reproducible_sample(1:N1_Gx,leng_1,5000*seedC),reproducible_sample((N1_Gx+1):(N1_Gx+N2_Gx),leng_2,6000*seedC))
S_mb <- c(reproducible_sample(1:N1_Gm,leng_3,7000*seedC),reproducible_sample((N1_Gm+1):(N1_Gm+N2_Gm),leng_4,8000*seedC))
S_xc <- c(reproducible_sample(1:N1_Gx,leng_1,9000*seedC),reproducible_sample((N1_Gx+1):(N1_Gx+N2_Gx),leng_2,10000*seedC))
S_mc <- c(reproducible_sample(1:N1_Gm,leng_3,11000*seedC),reproducible_sample((N1_Gm+1):(N1_Gm+N2_Gm),leng_4,12000*seedC))
S_xd  <- c(reproducible_sample(1:N1_Gx,leng_1,13000*seedC),reproducible_sample((N1_Gx+1):(N1_Gx+N2_Gx),leng_2,14000*seedC))
S_md  <- c(reproducible_sample(1:N1_Gm,leng_3,15000*seedC),reproducible_sample((N1_Gm+1):(N1_Gm+N2_Gm),leng_4,16000*seedC))

#===allele freq
p_Gx <- reproducible_runif(N_Gx,0.3,0.7,17000*seedC)
p_Gm <- reproducible_runif(N_Gm,0.3,0.7,18000*seedC)
p_Gk <- c(p_Gx,p_Gm)

#===sample size
n_gwasX <- opt$N
n_gwasM <- opt$N
n_gwasY <- opt$N

#===Scenarios set-up for S1-S11
if(opt$pleiotropy==2|opt$pleiotropy==5|opt$pleiotropy==6){
	b_xk[S_xb] <-  reproducible_runif(leng_1+leng_2,-0.15,0.15,19000*seedC)
	b_mk[S_mb] <-  reproducible_runif(leng_3+leng_4,-0.15,0.15,20000*seedC)
}

if(opt$pleiotropy==3|opt$pleiotropy==5|opt$pleiotropy==6){
	c_xk[S_xc] <- reproducible_runif(leng_1+leng_2,-0.15,0.15,21000*seedC)
	c_mk[S_mc] <- reproducible_runif(leng_3+leng_4,-0.15,0.15,22000*seedC)
}

if(opt$pleiotropy==4|opt$pleiotropy==5|opt$pleiotropy==6){
	d_x1k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.15,0.15,23000*seedC)
	d_x2k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.15,0.15,24000*seedC)
	d_x3k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.15,0.15,25000*seedC)

	d_m1k[S_md] <- reproducible_runif(leng_3+leng_4,-0.15,0.15,26000*seedC)
	d_m2k[S_md] <- reproducible_runif(leng_3+leng_4,-0.15,0.15,27000*seedC)
	d_m3k[S_md] <- reproducible_runif(leng_3+leng_4,-0.15,0.15,28000*seedC)	
}

if(opt$pleiotropy==7|opt$pleiotropy==10|opt$pleiotropy==11){
	b_xk[S_xb] <-  reproducible_runif(leng_1+leng_2,-0.1,0.2,29000*seedC)
	b_mk[S_mb] <-  reproducible_runif(leng_3+leng_4,-0.1,0.2,30000*seedC)
}

if(opt$pleiotropy==8|opt$pleiotropy==10|opt$pleiotropy==11){
	c_xk[S_xc] <- reproducible_runif(leng_1+leng_2,-0.1,0.2,31000*seedC)
	c_mk[S_mc] <- reproducible_runif(leng_3+leng_4,-0.1,0.2,32000*seedC)
}

if(opt$pleiotropy==9|opt$pleiotropy==10|opt$pleiotropy==11){
	d_x1k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.1,0.2,33000*seedC)
	d_x2k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.1,0.2,34000*seedC)
	d_x3k[S_xd] <- reproducible_runif(leng_1+leng_2,-0.1,0.2,35000*seedC)

	d_m1k[S_md] <- reproducible_runif(leng_3+leng_4,-0.1,0.2,36000*seedC)
	d_m2k[S_md] <- reproducible_runif(leng_3+leng_4,-0.1,0.2,37000*seedC)
	d_m3k[S_md] <- reproducible_runif(leng_3+leng_4,-0.1,0.2,38000*seedC)		
}

#===Pseudo GWAS data frame
gwasX <- data.frame(SNP.X=paste("rs",seq(N_Gk),sep=""),X=rep("X",N_Gk),id.X=rep("X",N_Gx),
effect_allele.X=rep("A",N_Gk), other_allele.X=rep("G",N_Gk), 
eaf.X=p_Gk, beta.X=rep(0,N_Gk), se.X=rep(0,N_Gk), pval.X=rep(0,N_Gk))

gwasM <- gwasX; gwasM$X <- "M"; gwasM$id.X <- "M"
gwasY <- gwasX; gwasY$X <- "Y"; gwasY$id.X <- "Y"

names(gwasM) <- gsub("X","M",names(gwasX))
names(gwasY) <- gsub("X","Y",names(gwasX))


#==============================#
#==========start Loop==========#
#==============================#

set.seed(1000)  #fixed seed for simulation
mylist <- list()
runtimelist <- list()
snpStren <- array(0,c(opt$R,6)) #SNPs strength 

for(r in 1:opt$R){

	#===generating X
	#===Gx'
	G_xk_gwasX <- array(0, c(n_gwasX,N_Gx))
	for(i in 1:N1_Gx){G_xk_gwasX[,i] <- rbinom(n_gwasX, 2, p_Gx[i])}
	#===Gm'
	G_mk_gwasX <- array(0, c(n_gwasX,N_Gm))
	for(i in 1:N1_Gm){G_mk_gwasX[,i] <- rbinom(n_gwasX, 2, p_Gm[i])}
	#===paired variables for Gx/Gx' Gm/Gm'
	for(i in (N1_Gx+1):(N1_Gx+N2_Gx)){tmp <- rcorrbinom(n_gwasX,prob1=p_Gx[i],prob2=p_Gm[i],corr=0.3);G_xk_gwasX[,i]<-tmp[,1];G_mk_gwasX[,i]<-tmp[,2]}

	epsilon_gwasX <- rmvnorm(n=n_gwasX, mean=c(0,0,0), sigma=Sig)
	U1_gwasX <- (epsilon_gwasX[,1] + G_xk_gwasX%*%d_x1k + G_mk_gwasX%*%d_m1k)
	U3_gwasX <- (epsilon_gwasX[,3] + G_xk_gwasX%*%d_x3k + G_mk_gwasX%*%d_m3k)
	Zx_gwasX <- U1_gwasX + U3_gwasX 
	X <- rnorm(n_gwasX,0,omega) + Zx_gwasX + G_xk_gwasX%*%a_xk + G_mk_gwasX%*%b_mk   #a_Gk contained b_k

	#===generating M
	#===Gx'
	G_xk_gwasM <- array(0, c(n_gwasM,N_Gx))
	for(i in 1:N1_Gx){G_xk_gwasM[,i] <- rbinom(n_gwasM, 2, p_Gx[i])}
	#===Gm'
	G_mk_gwasM <- array(0, c(n_gwasM,N_Gm))
	for(i in 1:N1_Gm){G_mk_gwasM[,i] <- rbinom(n_gwasM, 2, p_Gm[i])}
	#===paired variables for Gx/Gx' Gm/Gm'
	for(i in (N1_Gx+1):(N1_Gx+N2_Gx)){tmp <- rcorrbinom(n_gwasM,prob1=p_Gx[i],prob2=p_Gm[i],corr=0.3);G_xk_gwasM[,i]<-tmp[,1];G_mk_gwasM[,i]<-tmp[,2]}

	epsilon_gwasM <- rmvnorm(n=n_gwasM, mean=c(0,0,0), sigma=Sig)
	U1_gwasM <- (epsilon_gwasM[,1] + G_xk_gwasM%*%d_x1k + G_mk_gwasM%*%d_m1k)
	U2_gwasM <- (epsilon_gwasM[,2] + G_xk_gwasM%*%d_x2k + G_mk_gwasM%*%d_m2k)
	U3_gwasM <- (epsilon_gwasM[,3] + G_xk_gwasM%*%d_x3k + G_mk_gwasM%*%d_m3k)
	Zx_gwasM <- U1_gwasM + U3_gwasM 
	Zm_gwasM <- U1_gwasM + U2_gwasM 	
	X_gwasM <- rnorm(n_gwasM,0,omega) + Zx_gwasM + G_xk_gwasM%*%a_xk  + G_mk_gwasM%*%b_mk
	M <- rnorm(n_gwasM,0,omega) + Zm_gwasM + opt$alpha*X_gwasM + G_mk_gwasM%*%a_mk + G_xk_gwasM%*%b_xk

	#===generating Y
	#===Gx'
	G_xk_gwasY <- array(0, c(n_gwasY,N_Gx))
	for(i in 1:N1_Gx){G_xk_gwasY[,i] <- rbinom(n_gwasY, 2, p_Gx[i])}
	#===Gm'
	G_mk_gwasY <- array(0, c(n_gwasY,N_Gm))
	for(i in 1:N1_Gx){G_mk_gwasY[,i] <- rbinom(n_gwasY, 2, p_Gm[i])}
	#===paired variables for Gx/Gx' Gm/Gm'
	for(i in (N1_Gx+1):(N1_Gx+N2_Gx)){tmp <- rcorrbinom(n_gwasY,prob1=p_Gx[i],prob2=p_Gm[i],corr=0.3);G_xk_gwasY[,i]<-tmp[,1];G_mk_gwasY[,i]<-tmp[,2]}
	
	epsilon_gwasY <- rmvnorm(n=n_gwasY, mean=c(0,0,0), sigma=Sig)
	U1_gwasY <- (epsilon_gwasY[,1] + G_xk_gwasY%*%d_x1k + G_mk_gwasY%*%d_m1k )
	U2_gwasY <- (epsilon_gwasY[,2] + G_xk_gwasY%*%d_x2k + G_mk_gwasY%*%d_m2k )
	U3_gwasY <- (epsilon_gwasY[,3] + G_xk_gwasY%*%d_x3k + G_mk_gwasY%*%d_m3k )
	Zx_gwasY <- U1_gwasY + U3_gwasY 
	Zm_gwasY <- U1_gwasY + U2_gwasY
	Zy_gwasY <- U2_gwasY + U3_gwasY	
	X_gwasY <- rnorm(n_gwasY,0,omega) + Zx_gwasY + G_xk_gwasY%*%a_xk + G_mk_gwasY%*%b_mk 
	M_gwasY <- rnorm(n_gwasY,0,omega) + Zm_gwasY + opt$alpha*X_gwasY + G_mk_gwasY%*%a_mk + G_xk_gwasY%*%b_xk 
	Y <- rnorm(n_gwasY,0,omega) + Zy_gwasY + opt$delta*X_gwasY + opt$beta*M_gwasY + G_xk_gwasY%*%c_xk + G_mk_gwasY%*%c_mk

	#===use myGWAS() to generate gwas statistics
	Gk_X <- cbind(G_xk_gwasX,G_mk_gwasX)
	Gk_M <- cbind(G_xk_gwasM,G_mk_gwasM)
	Gk_Y <- cbind(G_xk_gwasY,G_mk_gwasY)
	if(opt$pleiotropy==6|opt$pleiotropy==11){
		gwasX[,c("beta.X","se.X","pval.X")] <- myGWAS(X,Gk_X)
	}else{
		gwasX[,c("beta.X","se.X","pval.X")] <- myGWAS(X_gwasM,Gk_M)
	}
	gwasM[,c("beta.M","se.M","pval.M")] <- myGWAS(M,Gk_M)
	gwasY[,c("beta.Y","se.Y","pval.Y")] <- myGWAS(Y,Gk_Y)

	#===categorize the IVs: Gx Gx' Gm Gm'
	ctg_Gk <- data.frame(SNP=gwasX$SNP,Gx=c(rep(1,N_Gx),rep(0,N_Gm)),Gx_plum=c(rep(1,N1_Gx),rep(0,N2_Gx+N_Gm)),Gm_plum=c(rep(0,N_Gx),rep(1,N1_Gm),rep(0,N2_Gm)),G_mvmr=c(rep(1,N1_Gx),rep(0,N2_Gx),rep(1,N1_Gm),rep(0,N2_Gm)))
	dat_sim <- cbind(ctg_Gk,gwasX,gwasM,gwasY)

	#rsq and averaged F-stat of IVs, computing based on sample of Y
	snpStren[r,1] <- var(G_xk_gwasY[,1:N1_Gx]%*%a_xk[1:N1_Gx] + G_xk_gwasY[,1:N1_Gx]%*%d_x1k[1:N1_Gx] + G_xk_gwasY[,1:N1_Gx]%*%d_x3k[1:N1_Gx])/var(X_gwasY)
	snpStren[r,2] <- var(G_mk_gwasY[,1:N1_Gm]%*%a_mk[1:N1_Gm] + G_mk_gwasY[,1:N1_Gm]%*%d_m1k[1:N1_Gm] + G_mk_gwasY[,1:N1_Gm]%*%d_m2k[1:N1_Gm])/var(M_gwasY)
	snpStren[r,3] <- mean((gwasX$beta.X[ctg_Gk$Gx_plum==1]/gwasX$se.X[ctg_Gk$Gx_plum==1])^2)
	snpStren[r,4] <- mean((gwasM$beta.M[ctg_Gk$Gm_plum==1]/gwasM$se.M[ctg_Gk$Gm_plum==1])^2)	
	indx <- which(dat_sim$G_mvmr==1)
	r_input <- format_mvmr(
		BXGs = dat_sim[indx,c("beta.X","beta.M")],
		BYG = dat_sim$beta.Y[indx],
		seBXGs = dat_sim[indx,c("se.X","se.M")],
		seBYG = dat_sim$se.Y[indx],RSID =dat_sim$SNP[indx])
	cF <- strength_mvmr(r_input, 0)
	snpStren[r,5] <- cF[[1]]
	snpStren[r,6] <- cF[[2]]
	colnames(snpStren)<- c("r2.SNP.X","r2.SNP.M","F.X","F.M","cF.X","cF.M")

	method_list=c("Diff_IVW_0","Diff_IVW","Diff_Egger","Diff_Median","Prod_IVW_0","Prod_IVW","Prod_Egger","Prod_Median")
	mylist[[r]] <- mrMed(dat_sim,method_list)

}


#==================#
#===output table===#
#==================#

IE_true <- opt$alpha*opt$beta
DE_true <- opt$delta 
TE_true <- IE_true + opt$delta 
rho_true <- IE_true/TE_true

#table for point estimate 
#we omit the rowname(method_list) and colnames (TE_mean,TE_mse,DE_mean,DE_mse,IE_mean,IE_mse,rho_mean,rho_mse)
table_mse <- array(0,c(length(method_list),8))
table_mse[,1] <- Reduce("+", lapply(lapply(mylist, "[[", c("TE")),function(x) x[,c("b")]))/ length(mylist)
table_mse[,2] <- Reduce("+", lapply(lapply(mylist, "[[", c("TE")),function(x) (x[,c("b")]-TE_true)^2))/length(mylist)
table_mse[,3] <- Reduce("+", lapply(lapply(mylist, "[[", c("DE")),function(x) x[,c("b")]))/ length(mylist)
table_mse[,4] <- Reduce("+", lapply(lapply(mylist, "[[", c("DE")),function(x) (x[,c("b")]-DE_true)^2))/length(mylist)
table_mse[,5] <- Reduce("+", lapply(lapply(mylist, "[[", c("IE")),function(x) x[,c("b")]))/ length(mylist)
table_mse[,6] <- Reduce("+", lapply(lapply(mylist, "[[", c("IE")),function(x) (x[,c("b")]-IE_true)^2))/length(mylist)
table_mse[,7] <- Reduce("+", lapply(lapply(mylist, "[[", c("rho")),function(x) x[,c("b")]))/ length(mylist)
table_mse[,8] <- Reduce("+", lapply(lapply(mylist, "[[", c("rho")),function(x) (x[,c("b")]-rho_true)^2))/length(mylist)


#table for ci with coverage rate (cr) and interval length (il)
#we omit the rowname(method_list) and colnames (TE_cr,TE_il,DE_cr,DE_il,IE_cr,IEil,rho_cr,rho_il)
table_ci <- array(0,c(length(method_list),8))
table_ci[,1] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("TE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) (x[,c("CI_lower")]<=TE_true)*(x[,c("CI_upper")]>=TE_true)))/ length(mylist)
table_ci[,2] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("TE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) x[,c("CI_upper")]-x[,c("CI_lower")]))/ length(mylist)
table_ci[,3] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("DE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) (x[,c("CI_lower")]<=DE_true)*(x[,c("CI_upper")]>=DE_true)))/ length(mylist)
table_ci[,4] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("DE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) x[,c("CI_upper")]-x[,c("CI_lower")]))/ length(mylist)
table_ci[,5] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("IE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) (x[,c("CI_lower")]<=IE_true)*(x[,c("CI_upper")]>=IE_true)))/ length(mylist)
table_ci[,6] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("IE")),function(x) x[,c("CI_lower","CI_upper")]),function(x) x[,c("CI_upper")]-x[,c("CI_lower")]))/ length(mylist)
table_ci[,7] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("rho")),function(x) x[,c("CI_lower","CI_upper")]),function(x) (x[,c("CI_lower")]<=rho_true)*(x[,c("CI_upper")]>=rho_true)))/ length(mylist)
table_ci[,8] <- Reduce("+", lapply(lapply(lapply(mylist, "[[", c("rho")),function(x) x[,c("CI_lower","CI_upper")]),function(x) x[,c("CI_upper")]-x[,c("CI_lower")]))/ length(mylist)

#average r-square, F-stat, conditional F-stat for SNPs
snpStat <- t(apply(snpStren,2,mean))


write.table(table_mse,file=paste("table_mse_",opt$alpha,"alpha_",opt$beta,"beta_",opt$delta,"delta_S",opt$plei,".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(table_ci,file=paste("table_ci_",opt$alpha,"alpha_",opt$beta,"beta_",opt$delta,"delta_S",opt$plei,".csv",sep=""),sep=",",row.names=FALSE,col.names=FALSE)
write.table(snpStat,file=paste("table_snpStren_",opt$alpha,"alpha_",opt$beta,"beta_",opt$delta,"delta_S",opt$plei,".csv",sep=""),sep=",",row.names=FALSE)