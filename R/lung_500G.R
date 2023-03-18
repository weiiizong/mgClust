#' LGRC lung disease data (Top 500 variant genes)
#'
#' Gene expression data are collected from Gene Expression Omnibus (GEO) GSE47460
#' and clinical information obtained from Lung Genomics Research Consortium
#' (https://ltrcpublic.com/). The data has gene expression of 500 genes for n=259 pateints,
#' three prognostic covariates (age, gender, BMI) and seven outomes ((FEV1, FVC, FEV1/FVC,
#' RV, WBC1, WBC4, BODE). See our paper for details.
#'
#' @docType data
#'
#' @usage data(lung_500G)
#'
#' @format a list containing gene expression matrix `G`, covariate matrix `X` and outcome
#' matrix `Y`
#'
#' @keywords datasets
#'
#' @source Lung Genomics Research Consortium (https://ltrcpublic.com/)
#'         Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
#'
#' @examples
#' data(lung_500G)
"lung_500G"
