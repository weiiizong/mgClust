#' Fit multivariate outcome guided clustering model with EM algorithm
#'
#' @param G The matrix for omics data. Rows are samples and columns are omic features.
#' @param X The matrix for confounding variables. Rows are samples and columns are confounding
#' variables.
#' @param Y The matrix for multivariate outcomes. Rows are samples and columns are outcome
#' variables.
#' @param gamma Initial Gamma matrix (i.e., coefficients in the disease subtyping model).
#' @param Sigma Initial Sigma matrix (i.e., covariance matrix).
#' @param B Initial B matrix (i.e., intercept matrix in the outcome association model).
#' @param C Initial C matrix (i.e., common covariate effect in the outcome association model).
#' @param K Number of subtypes/clusters.
#' @param lambdaB A user speficied tunning parameter for outcome sparsity.
#' @param lambdaG.seq The tunning parameter for omics feature sparsity. If NULL, it will
#' automatically select the best lambdaG in each EM step by cross validation (i.e., using
#' cv.glmnet to fit the disease subtyping model). If a vector is provided, cross validation will
#' be restricted to these values. If a numeric value is provided, the disease subtyping model
#' will be fitted with this lambdaG without internal cross validation. Default is NULL.
#' @param runs The maximum number of EM updates, Default is 200.
#' @param quite Whether report the process of EM step. Default is TRUE.
#' @param alphaG The elasticnet mixing parameter. 1 is Lasso penalty and 0 is Ridge penalty.
#' @param seed "fixed" or a numeric value. "fixed" set seed for each run = r to be r. "random"
#' will result in random seeds.
#'
#' @return A list object
#' \itemize{
#'  \item{\code{gamma}}{fitted Gamme matrix}
#'  \item{\code{Sigma}}{fitted Sigma matrix}
#'  \item{\code{B}}{fitted B matrix}
#'  \item{\code{C}}{fitted C matrix}
#'  \item{\code{K}}{K}
#'  \item{\code{lambdaB}}{lambdaB}
#'  \item{\code{lambdaG_path}}{selected lambdaG along the updating steps}
#'  \item{\code{BIC}}{BIC}
#'  \item{\code{R2_RMSE_table}}{a table of R^2 and RMSE for each outcome}
#'  \item{\code{pred_prob}}{predicted probability for belonging to each subgroup}
#'  \item{\code{pred_lb}}{predicted subtype/cluster assignement}
#'  \item{\code{pred_outcome}} {prediced outcome matrix}
#' }
#' @export fit_mgClust
#'
#' @examples
#' \dontrun{
#' X = lung_1000G$X
#' G = lung_1000G$G
#' Y = lung_1000G$Y
#' K = 4
#' n = nrow(G)
#' np = ncol(X)
#' nq = ncol(Y)
#' NG = ncol(G)
#'
#' set.seed(12345)
#' library(sparcl)
#' fit = lm(Y~X)
#' C = coef(fit)[-1,]
#' Y1 = Y-X%*%C
#' km.perm = KMeansSparseCluster.permute(Y1,K=K,nperms=20)
#' km.out = KMeansSparseCluster(Y1,K=K,wbounds=km.perm$bestw,nstart = 150)
#' Cs = km.out[[1]]$Cs
#' w_old = t(sapply(1:length(Cs), function(x){
#' a = rep(0,K)
#' a[Cs[x]] = 1
#' return(a)
#' }))
#' library(glmnet)
#' cvfit=cv.glmnet(x=G, y=w_old, family="multinomial", type.multinomial = "grouped", intercept=T,
#'  alpha = 0.5)
#' coef = coef(cvfit,s="lambda.1se")
#' gamma_new_matrix=sapply(1:K, function(x) as.numeric(coef[[x]]))[-1,]
#' gamma=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])[,-K]
#' center.ls = lapply(1:K, function(k) apply(Y1[which(km.out[[1]]$Cs==k),],2,mean))
#' B = do.call(cbind,center.ls)
#' Sigma = diag(abs(rnorm(ncol(Y))),nrow = nq)
#' initials_list = list(B=B, C=C, Sigma=Sigma, gamma=gamma)
#' full_fit = fit_mgClust(G, X, Y, gamma = initials_list$gamma, Sigma = initials_list$Sigma,
#' B = initials_list$B, C = initials_list$C, =4, lambdaB = 0.002436381, lambdaG.seq=0.1507389,
#' runs = 200, quite=FALSE, alphaG = 0.5, seed = "fixed")
#' }
fit_mgClust = function(G, X, Y, gamma, Sigma, B, C, K, lambdaB, lambdaG.seq=NULL,
                        runs = 200, quite=TRUE, alphaG = 0.5, seed = "fixed"){
  X=as.matrix(X)#fit seperate intercept as B
  G=as.matrix(G)
  Y=as.matrix(Y)

  if(nrow(G) != nrow(X) | nrow(G) != nrow(Y)){
    stop("Number of samples of G, X, Y should be matched!")
  }
  n = nrow(G)
  np = ncol(X)
  nq = ncol(Y)
  NG = ncol(G)

  criterion = c()
  criterion.pai = c()
  #criterion.w = c()
  criterion.gamma = c()
  criterion.Sigma = c()
  criterion.B = c()
  criterion.C = c()
  lambdaG = c()
  l=1
  dis=500
  while (dis>10^(-2) & l<=runs) {
    if(seed == "fixed"){
      set.seed(l)
    }else{
      print("random seed per run")
    }
    if(quite==FALSE){
      print(paste("Run EM",l,"Times",sep = " "))
    }
    t1 = Sys.time()
    gamma_old = gamma#dim = ((NG),(K-1))
    Sigma_old = Sigma#dim = (nq,nq)
    B_old = B#dim = (nq,K)
    C_old = C#dim = (np,np)


    #==E-STEP==#
    gamma_old=cbind(gamma_old,0)
    pai_old=exp(G %*% gamma_old)/rowSums(exp(G %*% gamma_old))

    #likelihood for sample i belongs to subtype k
    func_matrix_old = sapply(1:K, function(k) sapply(1:n, function(x) Rfast::dmvnorm(Y[x, ,drop=F ],as.numeric(B_old[,k]+X[x,,drop=F] %*% C_old),Sigma_old)))#dim = (n,K)
    func_matrix_old[which(func_matrix_old==0)]=1e-300

    #calculate the expected value of Z
    w_old=pai_old*func_matrix_old/rowSums(pai_old*func_matrix_old)
    w_old[which(w_old==0)]=1e-300


    #==M-STEP==#
    if(length(lambdaG.seq) == 1){
      lambdaG_new = lambdaG.seq
      fit=glmnet::glmnet(x=G, y=w_old, lambda=lambdaG.seq, family = "multinomial", type.multinomial = "grouped", intercept = F, alpha = alphaG)
      gamma_new_matrix=sapply(1:K, function(x) as.numeric(fit$beta[[x]]))
      gamma_new_matrix=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])
    }else{
      lambdaG_new = "NA"
      gamma_new_matrix = tryCatch({
        cvfit=glmnet::cv.glmnet(x=G, y=w_old, family="multinomial", type.multinomial = "grouped", intercept=F, alpha = alphaG, lambda = lambdaG.seq)
        lambdaG_new = cvfit$lambda.1se
        coef = coef(cvfit,s="lambda.1se")
        gamma_new_matrix=sapply(1:K, function(x) as.numeric(coef[[x]]))[-1,]#remove intercept
        gamma_new_matrix=sapply(1:K, function(x) gamma_new_matrix[,x]-gamma_new_matrix[,K])
      }, error =function(e) {
        return(gamma_old)
      })
    }
    pai_new=exp(G %*% gamma_new_matrix)/rowSums(exp(G %*% gamma_new_matrix))

    L_T_inv = Matrix::solve(Matrix::chol(Sigma_old))#Sigma=L%*%L_T
    Y1 = Y-X%*%C_old
    Y11 = Y1 %*% L_T_inv
    IMat_k = matrix(1,nrow = K,ncol = 1)
    Y1Vec_w = as.numeric(do.call(rbind,lapply(1:nq, function(j) sqrt(as.numeric(w_old))*(IMat_k %x% Y11[,j]))))
    Int_Mat_ls = list()
    for (j in 1:nq) {
      block_ls = list()
      for (k in 1:K) {
        block_ls[[k]] = as.matrix(sqrt(w_old[,k])) %x% t(as.matrix(L_T_inv[,j]))
      }
      Int_Mat_ls[[j]] = as.matrix(Matrix::bdiag(block_ls))
    }
    Int_Mat = do.call(rbind,Int_Mat_ls)

    group = rep(seq(1:nq),K)
    permute.group = group[order(group)]
    permute.Int_Mat = Int_Mat[,order(group)]
    data.SGL = list(x=permute.Int_Mat,y=Y1Vec_w)
    fit = SGL::SGL(data = data.SGL,index=permute.group,alpha=0,standardize=FALSE,lambdas = lambdaB)
    B_new = matrix(fit$beta,nrow = nq, ncol = K, byrow=T)

    #updata C
    Y2 = matrix(0, nrow = n, ncol = nq)
    for (i in 1:n) {
      w.res = lapply(1:K, function(k) w_old[i,k]*(Y[i,]-B_new[,k]))
      res.mat = do.call(rbind,w.res)
      Y2[i,] = apply(res.mat,2,sum)
    }
    C_new = Matrix::solve(t(X)%*%X) %*% (t(X)%*%Y2)

    #updata Sigma
    Sigma_new = Reduce("+",lapply(1:K, function(k) t(Y-matrix(1,nrow=n)%*%matrix(B_new[,k],nrow = 1)-X%*%C_new) %*% diag(w_old[,k]) %*% (Y-matrix(1,nrow=n)%*%matrix(B_new[,k],nrow = 1)-X%*%C_new)))/n

    #func_matrix_new = sapply(1:K, function(k) sapply(1:n, function(x) Rfast::dmvnorm(Y[x, ,drop=F ],as.numeric(B_new[,k]+X[x,,drop=F] %*% C_new),Sigma_new)))#dim = (n,K)
    #w_new=pai_new*func_matrix_new/rowSums(pai_new*func_matrix_new)

    dis_gamma = sum(as.numeric((gamma_new_matrix - gamma_old)^2),na.rm = T)
    dis_pai = sum(as.numeric((pai_new - pai_old)^2),na.rm = T)
    #dis_w = sum(as.numeric((w_new - w_old)^2),na.rm = T)

    dis_Sigma = sum(as.numeric((Sigma_new-Sigma_old)^2),na.rm = T)
    dis_B= sum(as.numeric((B_new-B_old)^2),na.rm = T)
    dis_C = sum(as.numeric((C_new-C_old)^2),na.rm = T)

    dis = sqrt(sum(dis_pai+dis_Sigma+dis_B+dis_C,na.rm = T))

    gamma = as.matrix(gamma_new_matrix[,-K])
    Sigma = Sigma_new
    B  = B_new
    C = C_new

    lambdaG[l] = lambdaG_new
    criterion[l] = dis
    criterion.gamma[l] = dis_gamma
    criterion.pai[l] = dis_pai
    #criterion.w[l] = dis_w
    criterion.Sigma[l] = dis_Sigma
    criterion.B[l] = dis_B
    criterion.C[l] = dis_C
    l=l+1
    t2 = Sys.time()

  }

  ##R2_RMSE table--------------
  gamma0=cbind(gamma,0)
  pai=sapply(1:K, function(k)
    exp(G %*% gamma0[,k,drop=F])/rowSums(exp(G %*% gamma0)))#predicted pai
  pred_lb = factor(apply(pai,1,which.max))
  Yhat = pai%*%t(B) + X%*%C

  RMSE.ind=apply(Y-Yhat,2,function(x) sqrt(sum(x^2)/length(x)))
  RSS.ind=apply(Y-Yhat,2,function(x) sum(x^2))
  Ybar = apply(Y,2,mean)
  Ybar_mat = matrix(rep(Ybar,nrow(Y)),nrow = nrow(Y), ncol = length(Ybar), byrow = T)
  TSS.ind=apply(Y-Ybar_mat,2,function(x) sum(x^2))
  R2.ind=1-RSS.ind/TSS.ind
  R2_RMSE_table = data.frame(RMSE.ind=RMSE.ind, R2.ind=R2.ind)

  ##BIC--------------
  kfunc = function(k) {
    sapply(1:n,function(x) (2*pi)^(-nq/2)*(det(Sigma))^(-1/2)*exp(-0.5*(as.matrix(Y[x,]- B[,k]-X[x,] %*% C))%*%Matrix::solve(Sigma)%*% t(as.matrix(Y[x,]- B[,k] -X[x,] %*% C))))
  }
  func_matrix = sapply(1:K, kfunc)
  l1 = sum(log(diag(pai %*% t(func_matrix))))
  par.num = length(which(as.numeric(B) !=0))+length(as.numeric(C))+nq*(nq-1)/2+length(which(as.numeric(gamma) !=0))
  BIC = log(n)*par.num-2*(c(l1))

  return(list(gamma=gamma, Sigma=Sigma, B=B, C=C,
              lambdaG_path = lambdaG,
              lambdaB = lambdaB,
              R2_RMSE_table = R2_RMSE_table,
              BIC = BIC,
              K=K,
              pred_prob = pai,
              pred_lb = pred_lb,
              pred_outcome = Yhat
              # criterion=criterion, criterion.gamma=criterion.gamma,
              # criterion.pai=criterion.pai,
              # criterion.Sigma=criterion.Sigma,criterion.B=criterion.B,criterion.C=criterion.C,
  ))
}
