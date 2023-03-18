#' Predict subtype/cluster assignment probabilities, labels and outcome matrix for new data.
#'
#' @param G The matrix for new omics data. Rows are samples and columns are omic features.
#' @param X The matrix for new confounding variables. Rows are samples and columns are confounding
#' variables. If NULL, only predicted subtype/cluster assignment probabilities and labels will
#' be generated.
#' @param Y The matrix for new multivariate outcomes. Rows are samples and columns are outcome
#' variables. If NULL, only predicted subtype/cluster assignment probabilities and labels will
#' be generated.
#' @param result The list output from function \code{fit_mgClust}.
#'
#' @return A list object
#' \itemize{
#'  \item{\code{pred_prob}}{predicted subtype/cluster assignment probabilities.}
#'  \item{\code{pred_lb}}{predicted subtype/cluster assignment labels.}
#'  \item{\code{pred_outcome}}{predicted outcome matrix. Only generated if both X and Y are
#'  provided.}
#' }
#' @export predict_mgClust
predict_mgClust = function(G, X = NULL, Y = NULL,result){
  G=as.matrix(G)
  n = nrow(G)
  NG = ncol(G)

  Sigma = result$Sigma
  B = result$B
  C = result$C
  K = result$K

  #pai: (n*K)
  gamma = result$gamma
  gamma=cbind(gamma,0)
  pai=sapply(1:K, function(k)
    exp(G %*% gamma[,k,drop=F])/rowSums(exp(G %*% gamma)))
  predict_lb = factor(apply(pai,1,which.max))

  res = list(pred_prob = pai,
             pred_lb = predict_lb)

  if(!(is.null(X) | is.null(Y))){
    X=as.matrix(X)
    Y=as.matrix(Y)
    if(nrow(G) != nrow(X) | nrow(G) != nrow(Y)){
      stop("Number of samples of G, X, Y should be matched!")
    }
    np = ncol(X)
    nq = ncol(Y)
    Yhat = pai%*%t(B) + X%*%C

    res$pred_outcome = Yhat
  }

  return(res)
}
