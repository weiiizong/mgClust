#' LGRC lung disease data (diagnosis information)
#'
#' The clinical diagnosis information of 259 subjects in lung_1000G and lung_500G.
#' The "Major_Diagnosis_Final_Clinical" contains the broad classification of samples as
#' "COPD" or "ILD" while the "diagnosis" contains detailed subtype information of ILD
#' patients.
#'
#' @docType data
#'
#' @usage data(lung_diagnosis)
#'
#' @format a list containing the broad ("Major_Diagnosis_Final_Clinical") and detailed
#' ("diagnosis") classification of the 259 patients in
#' lung_1000G and lung_500G data.
#'
#' @keywords datasets
#'
#' @source Lung Genomics Research Consortium (https://ltrcpublic.com/)
#'         Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
#'
#' @examples
#' data(lung_diagnosis)
"lung_diagnosis"
