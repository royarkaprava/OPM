#' @title The function to fit OPM algorithm for variable selection under two data-types Adaptive lasso in Step I and OPM in step 2
#' @description Takes the response Y (a binary vector of class labels) and predictor matrix X
#' @references Diego et. al. (2020)
#'     "Combining Phenotypic and Genomic Data]{Combining phenotypic and genomic data to improve prediction of binary traits"
#'
#' @param Y is the binary vector of class labels
#' @param X is the predictor matrix structures as [Macrovariable, Microvariables]
#' @param LAMBDA is the penalty term
#' @param ncores is to signify number of cores for parallelization
#' @param test is the index vector testing set
#' @param train is the index vector training set
#' @param IndMa is the number of Macrovariables

#' @return fit.ALASSOOPM returns a list of errors (prediction and estimation in that order) and the coefficients \cr



fit.ALASSOOPM <- function(Y, X, LAMBDA, ncores, test, train, IndMa){
  
  
  
  library(parallel)
  library(matrixcalc)
  library(Hmisc)
  library(Matrix)
  library(compiler)
  library(foreach)
  library(doSNOW)
  library(matrixcalc)
  library(mvtnorm)
  library(sn)
  library(glmnet)
  
  
  matrix_prod <- function(x, beta, m){
    r <- matrix(0, nrow(x), nrow(beta))
    t <- 0
    if(is.null(m))
      r <- r 
    else{
      m <- rbind(m)
      while(t < nrow(m)){
        t <- t + 1
        col.no <- m[t, 1]
        multi <- beta[col.no, m[t, 2]]
        r[, col.no] <- r[, col.no] + x[, m[t, 2]] * multi}}
    return(r)
  }
  
  
  loglike <- function(y, x, beta, lambda, m){
    alpha <- exp(matrix_prod(x, beta, m))
    pc <- as.matrix(cbind(1:length(y), y))
    y.zero <- (which(y == 0))
    if(length(y.zero) != 0)
    {pc <- pc[ - y.zero, ]}
    alpha1 <- matrix(0, dim(x)[1], dim(beta)[1])
    alpha1[pc] <- alpha[pc]
    dv <- 1 + rowSums(alpha)
    d <- matrix(rep(dv, dim(beta)[1]), dim(x)[1], dim(beta)[1])
    temp <- alpha1 / d
    l <- sum(log(temp[pc])) - sum(log(dv[y.zero])) - lambda*(sum(abs(beta[,-c(1:IndMa)])))
    if(is.nan(l)){l <- -Inf}
    if(is.na(l)){l <- -Inf}
    return(l)
  }
  logl <- cmpfun(loglike)
  prb <- function(x, beta){
    alpha <- exp(tcrossprod(x, beta))
    dv <- rowSums(alpha)
    d <- matrix(rep(dv, dim(beta)[1]), dim(x)[1], dim(beta)[1])
    probs <- alpha / d
    return(probs)
  }
  
  
  newton_compo <- function(l, j, y, x, beta, m, b){
    beta[l, j] = b
    alpha <- exp(matrix_prod(x, beta, m))
    c <- which(y == l)
    dv <- 1 + rowSums(alpha)
    f <- alpha[, l] / dv
    if(l == 0){f <- 1 / dv}
    t2 <- x[, j] * f
    f.fd <- sum(x[c, j]) - sum(t2)
    f.sd <- (x[, j] ^ 2) * f * (1 - f)
    f.sd <- - sum(f.sd)
    return(list(f.fd = f.fd, f.sd = f.sd))
  }
  
  nc <- cmpfun(newton_compo)
  
  newton <- function(l, j, X, Y, beta, lambda, m, tolfun = 1e-4)
  {
    A <- c(l, j)
    b <- beta[l, j]
    itrlogl <- logl(Y, X, beta, lambda, m)
    m1 <- as.matrix(rbind(m, A))
    oldlogl <- 0
    while(abs(itrlogl - oldlogl) > tolfun * (abs(oldlogl) + 1))
    {
      c <- b
      oldlogl <- itrlogl
      if(oldlogl == Inf){break}
      beta[l, j] = b
      alpha <- exp(matrix_prod(X, beta, m))
      cl <- which(Y == l)
      dv <- 1 + rowSums(alpha)
      f <- alpha[, l] / dv
      if(l == 0){f <- 1 / dv}
      t2 <- X[, j] * f
      f.fd <- sum(X[cl, j]) - sum(t2)
      f.sd <- (X[, j] ^ 2) * f * (1 - f)
      f.sd <- - sum(f.sd)
      NC <- list(f.fd = f.fd, f.sd = f.sd)
      itrlogl1 <- -Inf
      b1 <- ifelse(b - ((NC$f.fd - lambda) / NC$f.sd) > 0, b - ((NC$f.fd - lambda)) / NC$f.sd, ifelse(b - ((NC$f.fd + lambda) / NC$f.sd) < 0, b - ((NC$f.fd + lambda)) / NC$f.sd, break))
      
      beta[l, j] = b1
      itrlogl <- logl(Y, X, beta, lambda, m1)
      t <- 1
      if(is.nan(itrlogl)){break}
      if(is.na(itrlogl)){break}
      if(itrlogl == - Inf){break}
      while(itrlogl < oldlogl){
        b1 <- (b - (2^(-t))*((NC$f.fd - sign(b1) * lambda) / NC$f.sd))
        t <- t + 1
        beta[l, j] = b1
        itrlogl <- logl(Y, X, beta, lambda, m1)
        if(t > 11)
        {break}
        if(is.nan(itrlogl)){break}
        if(is.na(itrlogl)){break}
        if(itrlogl == - Inf){break}
      }
      b <- b1
      if(itrlogl < oldlogl){break}
    }
    if(abs(c) < .0001){c <- 0}
    return(list(est = c, ll = oldlogl))
  }
  
  ne <- cmpfun(newton)
  
  sol <- function(X, Y, lambda = 0, k)
  { 
    X <- xTRN
    y <- yTRN
    lambda <- LAMBDA
    p <- ncol(X)
    beta <- matrix(0, k, p)
    m <- cbind(1, 1:p)
    s <- k * p
    h <- rep(0, s)
    h1 <- rep(0, s)
    for(i in 2:k)
    {
      t <- cbind(i, 1:p)
      m <- rbind(m, t)
    }
    m <- as.matrix(m)
    m1 <- matrix(0, s, 2)
    m2 <- matrix(0, s, 1)
    g <- list(X, Y, beta, lambda, NULL)
    fun1 <- function(i, j, G)
    {
      X  <- G[[1]]
      Y <- G[[2]]
      beta <- G[[3]]
      lambda <- G[[4]]
      m <- G[[5]]
      r <- ne(i, j, X, Y, beta, lambda, m)
      return(r)
    }  
    
    for(d in 1:s)
    {
      if(d == s){break}
      
      if(d<=IndMa)
      {
        g[[4]] <- 0
        h <- mcmapply(fun1, m[d:IndMa, 1], m[d:IndMa, 2], MoreArgs = list(G = g), mc.cores=ncores)
        
        #print(c(lambda,d))
        h <- unlist(h)
        ll <- h[(2 * rep(1:(s - d + 1)))]
        est <- h[ - (2 * rep(1:(s - d + 1)))]
        a <- which.max(ll)
        h1 <- est[a]
        a <- a + d - 1
        m1[d, ] <- m[a, ]
        i <- m[a, 1]
        j <- m[a, 2]
        #print(c(i,j))
        
      }else{
        g[[4]] <- lambda
        h <- mcmapply(fun1, m[d:s, 1], m[d:s, 2], MoreArgs = list(G = g), mc.cores=ncores)
        
        #print(c(lambda,d))
        h <- unlist(h)
        ll <- h[(2 * rep(1:(s - d + 1)))]
        est <- h[ - (2 * rep(1:(s - d + 1)))]
        a <- which.max(ll)
        h1 <- est[a]
        a <- a + d - 1
        m1[d, ] <- m[a, ]
        i <- m[a, 1]
        j <- m[a, 2]
        #print(c(i,j))
        if(h1 == 0){break}
        
      }
      
      
      beta[i, j] = h1
      #print(h1)
      m <- m[ - a, ]
      m1[(d+1):s, ] <- m[d:(s - 1), ]
      m <- m1
      g <- list(X, Y, beta, lambda, m[1:d,])
    }
    return(beta)
  }
  
  
  
  
  n <- nrow(X)
  #X <- data[,-(1:6)]
  X2 <- X
  XM <- as.matrix(X[,1:IndMa])        
  Xm <- as.matrix(X[,-c(1:IndMa)])    
  p <- dim(Xm)[2]
  
  
  i <- 1
  tf <- matrix(NA,nrow=dim(XM)[1],ncol=dim(XM)[2])
  
  for(i in 1:IndMa)      
  {
    y <- XM[,i]  
    ridge1_cv <- cv.glmnet(x = Xm, y = y,
                           ## type.measure: loss to use for cross-validation.
                           type.measure = "mse",
                           ## K = 10 is the default.
                           nfold = 10,
                           ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                           alpha = 0)
    
    ## The intercept estimate should be dropped.
    best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
    ## Perform adaptive LASSO
    alasso1 <- glmnet(x = Xm, y = y,
                      ## type.measure: loss to use for cross-validation.
                      ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                      alpha = 1,
                      ##
                      ## penalty.factor: Separate penalty factors can be applied to each
                      ##           coefficient. This is a number that multiplies ‘lambda’ to
                      ##           allow differential shrinkage. Can be 0 for some variables,
                      ##           which implies no shrinkage, and that variable is always
                      ##           included in the model. Default is 1 for all variables (and
                      ##           implicitly infinity for variables listed in ‘exclude’). Note:
                      ##           the penalty factors are internally rescaled to sum to nvars,
                      ##           and the lambda sequence will reflect this change.
                      penalty.factor = 1 / abs(best_ridge_coef))
    ## Perform adaptive LASSO with 10-fold CV
    alasso1_cv <- cv.glmnet(x = Xm, y = y,
                            ## type.measure: loss to use for cross-validation.
                            type.measure = "mse",
                            ## K = 10 is the default.
                            nfold = 10,
                            ## ‘alpha = 1’ is the lasso penalty, and ‘alpha = 0’ the ridge penalty.
                            alpha = 1,
                            ##
                            ## penalty.factor: Separate penalty factors can be applied to each
                            ##           coefficient. This is a number that multiplies ‘lambda’ to
                            ##           allow differential shrinkage. Can be 0 for some variables,
                            ##           which implies no shrinkage, and that variable is always
                            ##           included in the model. Default is 1 for all variables (and
                            ##           implicitly infinity for variables listed in ‘exclude’). Note:
                            ##           the penalty factors are internally rescaled to sum to nvars,
                            ##           and the lambda sequence will reflect this change.
                            penalty.factor = 1 / abs(best_ridge_coef),
                            ## prevalidated array is returned
                            keep = TRUE)
    
    ## Extract coefficients at the error-minimizing lambda
    alasso1_cv$lambda.min
    
    best_alasso_coef1 <- coef(alasso1_cv, s = alasso1_cv$lambda.min)
    best_alasso_coef1
    
    ## https://en.wikipedia.org/wiki/Coefficient_of_determination
    r_squared <- function(y, yhat) 
    { 
      ybar <- mean(y)
      ## Total SS
      ss_tot <- sum((y - ybar)^2)
      ## Residual SS
      ss_res <- sum((y - yhat)^2)
      ## R^2 = 1 - ss_res/ ss_tot
      1 - (ss_res / ss_tot)
    }
    
    ## Function for Adjusted R^2
    ## n sample size, p number of prameters
    adj_r_squared <- function(r_squared, n, p) 
    {
      1 - (1 - r_squared) * (n - 1) / (n - p - 1)
    }
    
    ## Obtain R^2
    r_squared_alasso1 <- r_squared(as.vector(y), as.vector(predict(alasso1, newx = Xm, s = alasso1_cv$lambda.min)))
    r_squared_alasso1 
    yH <- as.vector(predict(alasso1, newx = Xm, s = alasso1_cv$lambda.min))
    
    res <- y-yH  
    cor(y,yH)
    tf[,i] <- res
    cat("Step 1: Macro-regression for ",i)
  }
  
  
  tf.sd <- scale(tf) 
  XM2 <- tf.sd
  X[,c(1:IndMa)] <- XM2   
  
  dim(X)
  x <- X[test, ]  
  y <- Y[test]          
  y <- as.vector(y)
  xte <- X[train,]  
  yte <- Y[train]   
  k=1
  
  xTRN <- xte 
  yTRN <- yte 
  
  xTST <- x 
  yTST <- y 
  
  s <- sol(xTRN,yTRN,LAMBDA/2,k)
  
  s <- as.matrix(s)
  s <- rbind(0, s)
  
  tf <- rep(NA,2)
  
  
  predprob <- prb(xTST, s)
  predyex <- as.vector(apply(predprob, 1, which.max))
  prederror <- sum((predyex - 1) != yTST)/length(yTST)
  
  tf[1] <- prederror
  
  prob <- prb(xte, s)
  yex <- as.vector(apply(prob, 1, which.max))
  error <- sum((yex - 1) != yte)/length(yte)
  #print(error)
  
  tf[2] <- error
  all.out <- matrix(tf,ncol=2)
  colnames(all.out) <- c('P_Error_B','P_Error_B_L')
  
  out <- list(errors = all.out, coefficients = s)
  
  return(out)
  
}

##################################################################################
########################  USE OF THE FUNCTION ####################################
##################################################################################



load("data.rda")
phenotypes2 <- Macro.micro_2_test
data <- phenotypes2   
Ymat <- data[, 1:6]
X    <- data[, -c(1:6)]
test <- 121:277
train <- 1:120
IndMa <- 20

out <- fit.ALASSOOPM(Ymat[, 2], X, LAMBDA = 0.5, ncores = 1, test, train, IndMa)

  

#################################################################################
#################################################################################



