\name{mrMed}
\alias{mrMed}
\title{
    a function to perform the MR-based mediation analysis using summary statistics from genome-wide association studies (GWAS)
}
\description{
    mrMed includes several methods to perform MR-based mediation analysis and provides corresponding estimates of the total effect (TE), direct effect (DE), indirect effect (IE), and mediation proportion (rho).
}
\usage{
    mrMed(dat_mrMed, method_list=c("Diff_IVW","Prod_IVW","Prod_Median"))
}
\arguments{
  \item{dat_mrMed}{
    dataframe, the required format are refered to the examples
}
  \item{method_list}{
    vector of characters, name(s) of the MR-based mediation methods
}

}

% \details{}

\value{
    \item{TE}{Results for the total effect including point estimate (b), standard error of the estimate (se), p-value (pval), and the upper bound and lower bound of 95\% confidence interval.}
    \item{DE}{Results for the direct effect including point estimate (b), standard error of the estimate (se), p-value (pval), and the upper bound and lower bound of 95\% confidence interval.}
    \item{IE}{Results for the indirect effect including point estimate (b), standard error of the estimate (se), p-value (pval), and the upper bound and lower bound of 95\% confidence interval.}
    \item{rho}{Results for the mediation proportion including point estimate (b), standard error of the estimate (se), p-value (pval), and the upper bound and lower bound of 95\% confidence interval.}
}
\references{
    Causal Mediation Analysis: A Summary-Data Mendelian Randomization Approach
}
\author{
    Shu-Chin Lin
}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{}

\examples{

data(WHR_T2D_CAD)
mrMed(dat_mrMed=WHR_T2D_CAD, method_list=c("Diff_IVW","Prod_IVW","Prod_Median"))
}
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line