% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lung_diagnosis.R
\docType{data}
\name{lung_diagnosis}
\alias{lung_diagnosis}
\title{LGRC lung disease data (diagnosis information)}
\format{
a list containing the broad ("Major_Diagnosis_Final_Clinical") and detailed
("diagnosis") classification of the 259 patients in
lung_1000G and lung_500G data.
}
\source{
Lung Genomics Research Consortium (https://ltrcpublic.com/)
        Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)
}
\usage{
data(lung_diagnosis)
}
\description{
The clinical diagnosis information of 259 subjects in lung_1000G and lung_500G.
The "Major_Diagnosis_Final_Clinical" contains the broad classification of samples as
"COPD" or "ILD" while the "diagnosis" contains detailed subtype information of ILD
patients.
}
\examples{
data(lung_diagnosis)
}
\keyword{datasets}
