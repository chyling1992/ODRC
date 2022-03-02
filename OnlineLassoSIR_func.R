require(MASS)
require(dr)
require(LassoSIR)
#require(progress)
#require(matlib)
require(glmnet)

# Description: generate covariance matrix
#
# Input parameters:
# - rho: parameter to control the covariance 
# - p: dimension of covariance matrix with param rho
#
# Output: \Sigma = (sigma_{i,j})_{1<= i,j <= p}, sigma_{i,j} = rho^{|i-j|}
#
S <- function(rho, p){
  a <- rho^(1:(p-1))
  b <- diag(1, p)
  for(i in 1:(p-1)){
    b[i,(i+1):p] <- a[1:(p-i)]
  }
  b <- b + t(b) - diag(1,p)
  b
}

# Description: data generated progress
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - beta: the coefficients of X
# - rho: parameter to control the covariance \Sigma, X ~ N(0,\Sigma)
# - type: choose which model is to be generated
#
# Output: (X_{i},Y_{i}),i=1,...,n
#
generate_data <- function(n, p, k, beta, rho, type){
  X = mvrnorm(n, mu = rep(0,p), Sigma = S(rho,p))
  epsilon = rnorm(n)
  Y = switch(type, 
             model1 = (X%*%beta) + epsilon,
             model2 = sin(X%*%beta)*exp(X%*%beta) + 0.1*epsilon,
             model3 = sign(X%*%beta[,1,drop=FALSE]) * (abs((X%*%beta[,2,drop=FALSE])/4 + 2))^3 + epsilon,
             model4 = 3*(X%*%beta[,1,drop=FALSE])*(2 + (X%*%beta[,2,drop=FALSE])/5)^2 + epsilon
             )
  
  list(X = X, Y = Y)
}

# Description: soft threshold for the W
soft_threshold <- function(w, theta, g, gamma){
  if(w > 0 && w <= theta){
    w_new = max(w - g*gamma, 0)
  }else if(w < 0 && w >= -theta){
    w_new = min(w + g*gamma, 0)
  }else{
    w_new = w
  }
  return(w_new)
}

# Description: metric in Cai2020
metric <- function(A,B){
  return( 1 - abs(det(t(A)%*%B)) )
}

# Description: \|Proj_{A} - Proj_{B}\|_F
metric_Frobenius <- function(A,B){
  return( norm(A%*%ginv(t(A)%*%A)%*%t(A)-B%*%ginv(t(B)%*%B)%*%t(B),"F") )
}
Proj <- function (y, X, list = FALSE) 
{
  if (is.vector(y)) 
    y <- matrix(y, ncol = 1)
  if (is.vector(X)) 
    X <- matrix(X, ncol = 1)
  XPX <- crossprod(X)
  P <- X %*% MASS::ginv(XPX) %*% t(X)
  if (!list) 
    return(c(P %*% y))
  else return(list(y = c(P %*% y), P = P))
}

# Description: schmidt orthogonal for x
schmidt_orthogonal <- function(x){
  n <- dim(x)[1]
  p <- dim(x)[2]
  z <- matrix(0,n,p)
  for(i in 1:p){
    if(norm(x[,i,drop=FALSE],"F") == 0){
      z[,i] = x[,i]
    }else{
      if(i == 1){ 
        z[,i] = x[,i]/norm(x[,i,drop=FALSE],"F")
      }else{
        z[,i] = x[,i]
        for(j in 1:(i-1)){
          z[,i] = z[,i] - Proj(x[,i], z[,j])
        }
        z[,i] = z[,i]/norm(x[,i,drop=FALSE],"F")
      }
    }

  }
  return(z)
}

schmidt_orthogonal_otp_norm <- function(x){
  n <- dim(x)[1]
  p <- dim(x)[2]
  z <- matrix(0,n,p)
  v <- array(p)
  for(i in 1:p){
    if(norm(x[,i,drop=FALSE],"F") == 0){
      z[,i] = x[,i]
    }else{
      v[i] = norm(x[,i,drop=FALSE],"F")
      if(i == 1){ 
        z[,i] = x[,i]/v[i]
      }else{
        z[,i] = x[,i]
        for(j in 1:(i-1)){
          z[,i] = z[,i] - Proj(x[,i], z[,j])
        }
        z[,i] = z[,i]/v[i]
      }
    }
    
  }
  list(z = z, v = v)
}

LassoSIR <-
  function( X, Y, H=0, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=10, screening=TRUE, no.dim=0)
  {
    if( no.dim!= 0)
      choosing.d="given"
    
    if( (categorical==FALSE)&(H==0))
    {
      H <- ( function()
      {
        H <- readline("For the continuous response, please choose the number of slices:   ")
        H <- as.numeric(unlist(strsplit(H, ",")))
        return( dim )
      }
      )()
    }
    
    p <- dim(X)[2]
    n <- dim(X)[1]
    
    
    if( categorical==FALSE)
    {
      ORD <- order( Y )
      X <- X[ORD, ]
      Y <- Y[ORD]
      
      ## Construct the  matrix M
      ms <- array(0, n)
      m <- floor( n/H )
      c <- n%%H
      M <- matrix(0, nrow=H, ncol=n )
      if( c==0 )
      {
        M <- diag( H ) %x% matrix( 1, nrow=1, ncol= m )/m
        ms <- m+ms
      }else{
        for(i in 1:c){
          M[i, ( (m+1)*(i-1)+1):( (m+1)*i )] <- 1/(m+1)
          ms[ ( (m+1)*(i-1)+1):( (m+1)*i) ] <- m
        }
        for( i in (c+1): H ){
          M[i, ( (m+1)*c + (i-c-1)*m +1):( (m+1)*c+(i-c)*m)] <- 1/m
          ms[ ( (m+1)*c +(i-c-1)*m+1):( (m+1)*c+(i-c)*m) ] <- m-1
        }
      }
      
      ## Calculate the screening statistics
      if( screening==TRUE){          
        x.sliced.mean <- M%*%X
        sliced.variance <- apply( x.sliced.mean, 2, var )
        keep.ind <- sort( order( sliced.variance, decreasing=TRUE)[1:n] )
      }else{
        keep.ind <- c(1:p)
      }
      
      X <- X[, keep.ind]
      
      X.H <- matrix(0, nrow=H, ncol= dim(X)[2] )
      grand.mean <- matrix( apply(X, 2, mean ), nrow=1, ncol=dim(X)[2] )
      X.stand.ord <- X - grand.mean %x% matrix(1, nrow=dim(X)[1], ncol=1)   
      X.H <- M%*% X.stand.ord
      
    }else{
      ms <- array(0, n)
      
      Y.unique <- unique(Y)
      H <- length( Y.unique )
      ORD <- which( Y==Y.unique[1] )      
      nH <- sum( Y == Y.unique[1] )
      ms[1:nH] <- nH
      
      for( i in 2:H )
      {
        ORD <- c( ORD, which(Y==Y.unique[i]) )
        nH <- c(nH, sum( Y==Y.unique[i] ) )
        ms[ (sum( nH[1:(i-1)])+1): sum(nH[1:i]) ] <- nH[i]
      }
      X <- X[ORD, ]
      
      ## Construct the matrix M
      M <- matrix( 0,  nrow=H, ncol=n )
      M[ 1, 1:nH[1] ] <- 1/nH[1]
      for(i in 2:H)
        M[ i, (sum(nH[1:(i-1)])+1):sum(nH[1:i]) ] <- 1/nH[i]
      
      ## Calculate the screening statistics
      if (screening == TRUE ){
        x.sliced.mean <- M%*%X
        sliced.variance <- apply( x.sliced.mean, 2, var )
        keep.ind <- sort( order( sliced.variance, decreasing=TRUE)[1:n] )
      }else{
        keep.ind <- c(1:p)
      }
      
      X <- X[, keep.ind]
      
      X.H <- matrix(0, nrow=H, ncol= dim(X)[2] )
      grand.mean <- matrix( apply(X, 2, mean ), nrow=1, ncol=dim(X)[2] )
      X.stand.ord <- X - grand.mean %x% matrix(1, nrow=dim(X)[1], ncol=1)
      X.H <- M%*% X.stand.ord
      
    }
    
    svd.XH <- svd(X.H, nv=p)
    res.eigen.value <- array( 0, p )
    res.eigen.value[ 1:dim(X.H)[1] ] <- (svd.XH$d)^2/H
    
    ## LambdaHat <- t(X.H) %*% X.H/H
    ## temp <- eigen( LambdaHat )
    ## res.eigen.value <- temp$values
    
    if( choosing.d=="manual")
    {
      plot( c(1:p), res.eigen.value, ylab="eigen values" )
      no.dim <- ( function()
      {
        dim <- readline("Choose the number of directions:   ")
        dim <- as.numeric(unlist(strsplit(dim, ",")))
        return( dim )
      }
      )()
      
    }
    if ( choosing.d=="automatic")
    {
      ## eigen.ratio <- res.eigen.value[1:(p-1)]/ res.eigen.value[2:p]
      ## dim.cand <- which( eigen.ratio >2 )
      ## eigen.diff <- res.eigen.value[1:(p-1)] - res.eigen.value[2:p]
      ## no.dim <- dim.cand[ max( which( eigen.diff[ dim.cand ]> 0.5) ) ]
      
      ## To search for unknown d, we start from the maximum value H
      beta.hat <- array(0, c(p, min(p, H) ) )
      Y.tilde <- array(0, c(n, min(p, H ) ) )
      
      for( ii in 1: min(p, H) )
      {
        eii <- matrix( 0, nrow= dim(svd.XH$v)[2], ncol=1 )
        eii[ii] <- 1
        eigen.vec <- solve( t(svd.XH$v), eii )
        
        Y.tilde[,ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.vec/( res.eigen.value[ii] )*matrix( 1/ms, nrow=n, ncol=1 )
      }
      
      mus <- array(0, min(p, H) )
      
      for(ii in 1: min(p, H ) )
      {
        lars.fit.cv <- cv.glmnet( X.stand.ord, Y.tilde[,ii], nfolds=nfolds )
        ## choose the one with the smallest cvm
        
        ind <- max( which( lars.fit.cv$cvm==min(lars.fit.cv$cvm) ) )
        if(ind==1)
          ind <- 2
        
        lambda <- lars.fit.cv$lambda[ind]
        mus[ ii ] <- lambda
        lars.fit <- glmnet(X.stand.ord, Y.tilde[,ii], lambda=lambda)
        beta.hat[ keep.ind, ii ] <- as.double( lars.fit$beta )
      }
      
      ## The statistic for determining d is ||beta_i|| * lambda_i
      temp.2 <- sqrt( apply( beta.hat^2, 2, sum) )* res.eigen.value[1:H] 
      temp <- temp.2/temp.2[1]
      
      
      res.kmeans <- kmeans( temp, centers=2 )
      no.dim <- min( sum( res.kmeans$cluster==1), sum( res.kmeans$cluster==2 ) )
    }
    
    beta.hat <- array(0, c(p, no.dim) )
    Y.tilde <- array(0, c(n, no.dim) )
    
    
    for( ii in 1:no.dim)
    {
      eii <- matrix( 0, nrow= dim( t(svd.XH$v) )[2], ncol=1 )
      eii[ii] <- 1
      eigen.vec <- solve( t(svd.XH$v), eii )
      
      Y.tilde[,ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.vec/( res.eigen.value[ii] )*matrix( 1/ms, nrow=n, ncol=1 )
    }
    
    
    if( solution.path==FALSE )
    {
      mus <- array( 0, no.dim )
      
      if(no.dim==1){
        lars.fit.cv <- cv.glmnet( X.stand.ord, Y.tilde, nfolds=nfolds )
        
        ## choose the one with the smallest cvm
        ind <- max( which( lars.fit.cv$cvm==min(lars.fit.cv$cvm) ) )
        if(ind==1)
          ind <- 2
        
        lambda <- lars.fit.cv$lambda[ind]
        mus = lambda
        lars.fit <- glmnet(X.stand.ord, Y.tilde, lambda=lambda)
        beta.hat[ keep.ind ] <- as.double( lars.fit$beta )
      }else{
        for(ii in 1:no.dim)
        {
          lars.fit.cv <- cv.glmnet( X.stand.ord, Y.tilde[,ii], nfolds=nfolds )
          ## choose the one with the smallest cvm
          
          ind <- max( which( lars.fit.cv$cvm==min(lars.fit.cv$cvm) ) )
          if(ind==1)
            ind <- 2
          
          lambda <- lars.fit.cv$lambda[ind]
          mus[ ii ] <- lambda
          lars.fit <- glmnet(X.stand.ord, Y.tilde[,ii], lambda=lambda)
          beta.hat[ keep.ind, ii ] <- as.double( lars.fit$beta )
        }
      }
      ## list( beta= beta.hat, eigen.value=res.eigen.value, no.dim=no.dim, keep.ind=keep.ind, H=H, choosing.d=choosing.d, categorical=categorical, nfolds=nfolds, screening=screening, mu.lasso=mus )
      list( beta= beta.hat, eigen.value=res.eigen.value, no.dim=no.dim, H=H, categorical=categorical, mus = mus, Y.tilde = Y.tilde[sort(ORD),], X.stand.ord = X.stand.ord[sort(ORD),], X = X[sort(ORD),] )
      
    }else{
      lars.fit.all <- list()
      for(ii in 1:no.dim)
      {
        lars.fit.all[[ii]] <- glmnet( X.stand.ord, Y.tilde[,ii] )
      }
      lars.fit.all
    }
  }

# Description: Online Lasso-SIR with Perturbation
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
# - L: do trucated gradient every L steps
# - gamma: learning rate
#
# Output: error
#
#X = realdata_X; Y = realdata_Y
#X = data_simulation$X;Y = data_simulation$Y;n=1000;p=1000;k=1; batch_num=200;L = 10;
OnlineLassoSIR_Perturbation <- function(X, Y, n, p, k, H, batch_num = 200, L = 10, gamma, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Mt = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Mt[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Mt = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Mt[,i] = as.numeric(Yt == i)
    }  
  }
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  Mx = (1/batch_num)*t(Mt) %*% Xt_tilde
  Dt_hat = (1/H)*t(Mx) %*% Mx
  
  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values
  Vt = Dt_eig$vectors
  
  ## Initial the beta by LassoSIR with a small batch samples
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Dt_sum = matrix(0, nrow = p, ncol = p)
  Dt_sum = Dt_sum + Dt_hat
  Gamma_hat_t = Dt_sum
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]
  
  if(p >= 100){
    theta = sir.lasso$mus
  }else{
    theta = rep(1,k)
  }
  g_small = 0.01#0.01(model 1 model 2)
  #gamma = rep(0.1/n,k)#0.0001(model 1 model 2)
  beta_hat = sir.lasso$beta
  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  i = 1
  while (i <= n-batch_num){
    t = batch_num+i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new_class] = 1
    }else{
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new] = 1
    }
    ## update D_{t+1} i.e Dt_hat_update
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    Mx_new = ((t-1)*Mx + t(Mt_new) %*% X_new_central)*(1/t)
    Dt_hat_update = (1/H)*t(Mx_new)%*%Mx_new
    
    ## Use Perturbation to update the eigenvalues and eigen vectors of cov(E(X|Y))
    temp_mat_1 = Gamma_hat_t - Dt_hat_update
    for(j in 1:k){
      Lamdat_update[j] = Lamdat[j] - (i+1)^(-1)*t(Vt[,j,drop=FALSE])%*%temp_mat_1%*%Vt[,j,drop=FALSE]
      Vt_update[,j] = try(Vt[,j,drop=FALSE] - (i+1)^(-1)*ginv(Lamdat[j]*diag(p) - Gamma_hat_t)%*%temp_mat_1%*%Vt[,j,drop=FALSE])
      #if( !is.null(attr(tmp_vt, "class")) && attr(tmp_vt, "class") == "try-error"){
      #  return(tmp_vt)
      #}else{
      #  Vt_update[,j] = tmp_vt
      #}
      
      Vt_update[,j] = Vt_update[,j,drop=FALSE]/norm(Vt_update[,j,drop=FALSE], "F")
    }
    
    ## Online SIR via Truncated Gradient
    if(k == 1){
      Yt_new = (1/(Lamdat_update*t*H))*Mt_new%*%Mx_new%*%Vt_update
    }
    else{
      Yt_new = (1/(t*H))*Mt_new%*%Mx_new%*%Vt_update%*%diag(1/Lamdat_update)
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }else{
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }
      if( i == n-batch_num){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
      }
      #cat("i:",i,"\n")
    }
    
    
    ## calculate the distance between B and B_hat
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius( beta, beta_hat )
    
    ## iteration
    Xt_mean = Xt_mean_new
    Mx = Mx_new
    Dt_sum = Dt_sum + Dt_hat_update
    Gamma_hat_t = Dt_sum/(i+1)
    Lamdat = Lamdat_update
    Vt = Vt_update
    i = i+1
    
  }
  error_lasso_sir_initial = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )
  error_F_lasso_sir_initial = metric_Frobenius( beta, sir.lasso$beta )
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_lasso_sir_initial = error_lasso_sir_initial,
       error_F_lasso_sir_initial = error_F_lasso_sir_initial,
       error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}

# Description: Online Lasso-SIR with GD
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
# - L: do trucated gradient every L steps
# - gamma: learning rate
#
# Output: error
#
#X = data_simulation$X; Y = data_simulation$Y; n = 1000; p = 20; k = 1; batch_num = 100; L = 10; gamma = rep(0.01, k); categorical = FALSE
OnlineLassoSIR_GD <- function(X, Y, n, p, k, H, batch_num = 100, L = 10, gamma, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Mt = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Mt[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Mt = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Mt[,i] = as.numeric(Yt == i)
    }  
  }
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  Mx = (1/batch_num)*t(Mt) %*% Xt_tilde
  Dt_hat = (1/H)*t(Mx) %*% Mx
  
  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values[1:k]
  Vt = Dt_eig$vectors[,1:k,drop=FALSE]
  
  ## Initial the beta by LassoSIR with a small batch samples
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]
  
  if(p >= 100){
    theta = sir.lasso$mus
  }else{
    theta = rep(1,k)
  }
  
  g_small = 0.01#0.01(model 1 model 2)
  #gamma = rep(0.1/n,k)#0.0001(model 1 model 2)
  beta_hat = sir.lasso$beta
  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  i = 1
  while (i <= n-batch_num){
    t = batch_num+i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new_class] = 1
    }else{
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new] = 1
    }
    ## update D_{t+1} i.e Dt_hat_update
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    Mx_new = ((t-1)*Mx + t(Mt_new) %*% X_new_central)*(1/t)
    Dt_hat_update = (1/H)*t(Mx_new)%*%Mx_new
    
    ## Use GD to update the eigenvalues and eigen vectors of cov(E(X|Y))
    Vt_update = Vt + (i+1)^(-1/2)*(Dt_hat_update%*%Vt)
    Vt_update = schmidt_orthogonal(Vt_update)
    Lamdat_update = Lamdat + (i+1)^(-1/2)*(diag(t(Vt)%*%Dt_hat_update%*%Vt) - Lamdat)

    ## Online SIR via Truncated Gradient
    if(k == 1){
      Yt_new = (1/(Lamdat_update*t*H))*Mt_new%*%Mx_new%*%Vt_update
    }
    else{
      Yt_new = (1/(t*H))*Mt_new%*%Mx_new%*%Vt_update%*%diag(1/Lamdat_update)
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }else{
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }
      if( i == n-batch_num){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
      }
    }
    
    #beta_hat = schmidt_orthogonal(beta_hat)
    ## calculate the distance between B and B_hat
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius( beta, beta_hat )
    
    ## iteration
    Xt_mean = Xt_mean_new
    Mx = Mx_new
    Lamdat = Lamdat_update
    Vt = Vt_update    
    i = i+1
    
  }
  error_lasso_sir_initial = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )
  error_F_lasso_sir_initial = metric_Frobenius( beta, sir.lasso$beta )
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_lasso_sir_initial = error_lasso_sir_initial,
       error_F_lasso_sir_initial = error_F_lasso_sir_initial,
       error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}


# Description: Online Lasso-SIR with IPCA
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
# - L: do trucated gradient every L steps
# - gamma: learning rate
#
# Output: error
#
OnlineLassoSIR_IPCA <- function(X, Y, n, p, k, H, batch_num = 200, L = 10, gamma, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Mt = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Mt[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Mt = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Mt[,i] = as.numeric(Yt == i)
    }  
  }
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  Mx = (1/batch_num)*t(Mt) %*% Xt_tilde
  Dt_hat = (1/H)*t(Mx) %*% Mx
  
  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values[1:k]
  Vt = Dt_eig$vectors[,1:k,drop=FALSE]
  
  ## Initial the beta by LassoSIR with a small batch samples
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]
  
  if(p >= 100){
    theta = sir.lasso$mus
  }else{
    theta = rep(1,k)
  }
  g_small = 0.01#0.01(model 1 model 2)
  #gamma = rep(0.1/n,k)#0.0001(model 1 model 2)
  beta_hat = sir.lasso$beta
  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  i = 1
  while (i <= n-batch_num){
    t = batch_num+i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new_class] = 1
      id_update = Y_new_class
    }else{
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new] = 1
      id_update = Y_new
    }
    ## update D_{t+1} i.e Dt_hat_update
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    Mx_new = ((t-1)*Mx + t(Mt_new) %*% X_new_central)*(1/t)
    Dt_hat_update = (1/H)*t(Mx_new)%*%Mx_new
    
    ## Use Incremental PCA to update the eigenvalues and eigen vectors of cov(E(X|Y))
    Mx_new_k = t(Mx_new[id_update,,drop=FALSE])
    ut = Mx_new_k - Vt[,1:k,drop=FALSE]%*%t(Vt[,1:k,drop=FALSE])%*%Mx_new_k
    eta = cbind(Vt[,1:k,drop=FALSE],ut/norm(ut,"F"))
    Qt_update = t(eta)%*%Dt_hat_update%*%eta
    Qt_update_eig = eigen(Qt_update)
    Vt_update[,1:k] = eta%*%Qt_update_eig$vectors[,1:k]
    Lamdat_update[1:k] = Qt_update_eig$values[1:k]
    
    ## Online SIR via Truncated Gradient
    if(k == 1){
      Yt_new = (1/(Lamdat_update*t*H))*Mt_new%*%Mx_new%*%Vt_update
    }
    else{
      Yt_new = (1/(t*H))*Mt_new%*%Mx_new%*%Vt_update%*%diag(1/Lamdat_update)
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }else{
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*sqrt(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }
      if( i == n-batch_num){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
      }
    }
    
    
    ## calculate the distance between B and B_hat
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius( beta, beta_hat )
    
    ## iteration
    Xt_mean = Xt_mean_new
    Mx = Mx_new
    Lamdat = Lamdat_update
    Vt = Vt_update    
    i = i+1
    
  }
  error_lasso_sir_initial = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )
  error_F_lasso_sir_initial = metric_Frobenius( beta, sir.lasso$beta )
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_lasso_sir_initial = error_lasso_sir_initial,
       error_F_lasso_sir_initial = error_F_lasso_sir_initial,
       error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}


# Description: Online Lasso-SIR with CCIPCA
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
# - L: do trucated gradient every L steps
# - gamma: learning rate
#
# Output: error
#
OnlineLassoSIR_CCIPCA <- function(X, Y, n, p, k, H, batch_num = 200, L = 10, gamma, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Mt = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Mt[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Mt = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Mt[,i] = as.numeric(Yt == i)
    }  
  }
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  Mx = (1/batch_num)*t(Mt) %*% Xt_tilde
  Dt_hat = (1/H)*t(Mx) %*% Mx
  
  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values[1:k]
  Vt = Dt_eig$vectors
  ut = matrix(0,p,k)
  for(j in 1:k){
    ut[,j] = Lamdat[j]*Vt[,j,drop=FALSE]
  }  
  ## Initial the beta by LassoSIR with a small batch samples
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]
  ut_update = ut[,1:k,drop=FALSE]
  
  if(p >= 100){
    theta = sir.lasso$mus
  }else{
    theta = rep(1,k)
  }
  g_small = 0.01#0.01(model 1 model 2)
  #gamma = rep(0.1/n,k)#0.0001(model 1 model 2)
  beta_hat = sir.lasso$beta
  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  i = 1
  while (i <= n-batch_num){
    t = batch_num+i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new_class] = 1
    }else{
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new] = 1
    }
    ## update D_{t+1} i.e Dt_hat_update
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    Mx_new = ((t-1)*Mx + t(Mt_new) %*% X_new_central)*(1/t)
    Dt_hat_update = (1/H)*t(Mx_new)%*%Mx_new
    
    ## Use CCIPCA to update the eigenvalues and eigen vectors of cov(E(X|Y))
    #psi_t = Mx_new%*%ut
    #ut_update = (i/(i+1))*ut + (1/(H*(i+1)))*t(Mx_new)%*%psi_t%*%diag(1/Lamdat_update, k, k)
    ut_update = (i/(i+1))*ut + (1/(i+1))*Dt_hat_update%*%ut%*%diag(1/Lamdat_update, k, k)
    ut_update_schmidt = schmidt_orthogonal_otp_norm(ut_update)
    Vt_update = ut_update_schmidt$z
    Lamdat_update = ut_update_schmidt$v
    
    ## Online SIR via Truncated Gradient
    if(k == 1){
      Yt_new = (1/(as.numeric(Lamdat_update)*t*H))*Mt_new%*%Mx_new%*%Vt_update
    }
    else{
      Yt_new = (1/(t*H))*Mt_new%*%Mx_new%*%Vt_update%*%diag(1/Lamdat_update)
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }else{
        beta_hat[,jj] = beta_hat[,jj] + sqrt(log(p))*gamma[jj]*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      }
      if( i == n-batch_num){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta[jj], g_small*L, gamma[jj]*(1/t)^(1/2))
        if (sum(beta_vec != 0) <= 1) next
        beta_hat[,jj] = matrix(beta_vec,p,1)
      }
    }
    
    #beta_hat = schmidt_orthogonal(beta_hat)
    ## calculate the distance between B and B_hat
    #beta_hat = schmidt_orthogonal(beta_hat)
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius( schmidt_orthogonal(beta), beta_hat )
    
    ## iteration
    Xt_mean = Xt_mean_new
    Mx = Mx_new
    ut = ut_update
    Lamdat = Lamdat_update
    Vt = Vt_update    
    i = i+1
    
  }
  error_lasso_sir_initial = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )
  error_F_lasso_sir_initial = metric_Frobenius( beta, sir.lasso$beta )
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_lasso_sir_initial = error_lasso_sir_initial,
       error_F_lasso_sir_initial = error_F_lasso_sir_initial,
       error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}


# Description: Online SIR with Perturbation
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
#
# Output: error
#
OnlineSIR_Perturbation <- function(X, Y, n, p, k, H, batch_num = 200, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Yt_tilde = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Yt_tilde[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Yt_tilde = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Yt_tilde[,i] = as.numeric(Yt == i)
    }
  }
  Xt_tilde = cbind(1, Xt)
  At = t(Xt_tilde)%*%Xt_tilde
  At_inverse = try(ginv(At), TRUE)
  if(!is.null(attr(At_inverse,"class")) && attr(At_inverse,"class") == "try-error"){
    next
  }
  XYt = (t(Xt_tilde)%*%Yt_tilde)
  I_tilde = cbind(as.matrix(rep(0,p)),diag(p))
  mt_hat = I_tilde%*%At_inverse%*%XYt
  Mt_hat = mt_hat%*%t(mt_hat)
  Mt_eig = eigen(Mt_hat)
  Lamdat = Mt_eig$values
  Vt = Mt_eig$vectors
  
  Mt_sum = matrix(0, nrow = p, ncol = p)
  Mt_sum = Mt_sum + Mt_hat
  Gamma_hat_t = Mt_sum
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]

  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  i = 1
  while (i <= n-batch_num){
    t = batch_num + i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    xt1_tilde = cbind(1, X_new)
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      yt1_tilde = matrix(0, nrow = 1, ncol = H)
      yt1_tilde[1, Y_new_class] = 1
    }else{
      yt1_tilde = matrix(0, nrow = 1, ncol = H)
      yt1_tilde[1, Y_new] = 1
    }
    At1_inverse = At_inverse - (At_inverse%*%t(xt1_tilde)%*%xt1_tilde%*%At_inverse)/as.numeric(1 + xt1_tilde%*%At_inverse%*%t(xt1_tilde))
    Ct1 = I_tilde%*%At1_inverse
    XYt1 = XYt + t(xt1_tilde)%*%yt1_tilde
    mt1_hat = Ct1%*%XYt1
    Mt1_hat = mt1_hat%*%t(mt1_hat)
    
    # online SVD
    temp_mat_1 = Gamma_hat_t - Mt1_hat
    for(j in 1:k){
      Lamdat_update[j] = Lamdat[j] - (i+1)^(-1)*t(Vt[,j,drop=FALSE])%*%temp_mat_1%*%Vt[,j,drop=FALSE]
      Vt_update[,j] = Vt[,j,drop=FALSE] - (i+1)^(-1)*ginv(Lamdat[j]*diag(p) - Gamma_hat_t)%*%temp_mat_1%*%Vt[,j,drop=FALSE]
      Vt_update[,j] = Vt_update[,j,drop=FALSE]/norm(Vt_update[,j,drop=FALSE], "F")
    }
    beta_hat = Vt_update
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius(beta, beta_hat)
    # iteration
    At_inverse = At1_inverse
    XYt = XYt1
    Mt_sum = Mt_sum + Mt1_hat
    Gamma_hat_t = Mt_sum/(i+1)
    Lamdat = Lamdat_update
    Vt = Vt_update
    i = i+1
  }
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}


# Description: Online SIR with GD
#
# Input parameters:
# - n: the size of sample
# - p: the dimension of explanatory varaibles X
# - k: the number of directions
# - batch_num: a small batch for warm start
#
# Output: error
#
# X = train_set_X;Y = as.numeric(train_set_Y);n = n_train_use;
OnlineSIR_GD <- function(X, Y, n, p, k, H, batch_num = 200, categorical = FALSE, beta){
  ## pre-specify the cutting points
  ## and use a small batch to calculate covariance Cov(E(x|y)) for the warm start
  batch_ind = 1:batch_num # using a small batch
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    #H = 5
    Yt_tilde = matrix(0, nrow = batch_num, ncol = H)
    sliced_class = findInterval(Yt,slices_interval, left.open = TRUE)
    for (i in 1:H){
      Yt_tilde[,i] = as.numeric(sliced_class == i)
    }
  }else{
    class_num = length(unique(Y))
    H = class_num
    Yt_tilde = matrix(0, nrow = batch_num, ncol = H)
    for (i in 1:H){
      Yt_tilde[,i] = as.numeric(Yt == i)
    }
  }
  #sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical, screening=FALSE,no.dim=k)
  
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = cbind(1, Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean)
  At = t(Xt_tilde)%*%Xt_tilde
  At_inverse = try(ginv(At), TRUE)
  if(!is.null(attr(At_inverse,"class")) && attr(At_inverse,"class") == "try-error"){
    next
  }
  XYt = (t(Xt_tilde)%*%Yt_tilde)
  I_tilde = cbind(as.matrix(rep(0,p)),diag(p))
  mt_hat = I_tilde%*%At_inverse%*%XYt
  Mt_hat = mt_hat%*%t(mt_hat)
  Mt_eig = eigen(Mt_hat)
  Lamdat = Mt_eig$values
  Vt = Mt_eig$vectors
  Bt_hat = Vt[,1:k,drop=FALSE]
  #Bt_hat = matrix(runif(p*k),p,k)#sir.lasso$beta
  #Vt_update = Vt[,1:k,drop=FALSE]

  error_vec = rep(NA,n-batch_num)
  error_vec_F = rep(NA,n-batch_num)
  
  i = 1
  while (i <= n-batch_num){
    
    t = batch_num + i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    xt1_tilde = cbind(1, X_new - Xt_mean_new)
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      yt1_tilde = matrix(0, nrow = 1, ncol = H)
      yt1_tilde[1, Y_new_class] = 1
    }else{
      yt1_tilde = matrix(0, nrow = 1, ncol = H)
      yt1_tilde[1, Y_new] = 1
    }
    At1_inverse = At_inverse - (At_inverse%*%t(xt1_tilde)%*%xt1_tilde%*%At_inverse)/as.numeric(1 + xt1_tilde%*%At_inverse%*%t(xt1_tilde))
    Ct1 = I_tilde%*%At1_inverse
    XYt1 = XYt + t(xt1_tilde)%*%yt1_tilde
    mt1_hat = Ct1%*%XYt1
    Mt1_hat = mt1_hat%*%t(mt1_hat)
    gamma = sqrt(1/i) 
    Bt1_hat = schmidt_orthogonal(Bt_hat + gamma*Mt_hat%*%Bt_hat)
    beta_hat = Bt1_hat
    
    error_vec[i] = metric( schmidt_orthogonal(beta), schmidt_orthogonal(beta_hat) )
    error_vec_F[i] = metric_Frobenius(beta, beta_hat)
    # iteration
    Mt_hat = Mt1_hat
    Bt_hat = Bt1_hat
    At_inverse = At1_inverse
    XYt = XYt1
    i = i+1
  }
  
  #error_lasso_sir_initial = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )
  #error_F_lasso_sir_initial = metric_Frobenius( beta, sir.lasso$beta )
  error_min = min(error_vec,na.rm = TRUE)
  error_tail = tail(na.omit(error_vec),1)
  error_F_min = min(error_vec_F,na.rm = TRUE)
  error_F_tail = tail(na.omit(error_vec_F),1)
  list(error_min = error_min,
       error_tail = error_tail,
       error_F_min = error_F_min,
       error_F_tail = error_F_tail,
       beta_hat = beta_hat)
}