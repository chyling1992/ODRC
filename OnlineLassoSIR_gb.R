## online SIR via Truncated Gradient
## Cheng Haoyang
## 2020.6.20

set.seed(1234567890)
N = 100
error_mat = array(0, N)
error_tail_mat = array(0, N)
pb <- progress_bar$new(total = N)  #n-batch_num
time_start = proc.time()
for(s in 1:N){
  ## input data
  n = 1000
  p = 20
  k = 1 # the number of directions
  X = mvrnorm(n, mu = rep(0,p), Sigma = S(0.3,p))
  beta = matrix(0, nrow = p, ncol = k)
  beta[1:2,1] = 1 # model 1
  #beta[3,1] = 1 # model 2
  #beta[1,1] = 1 # model 3
  #beta[2,2] = 1 # model 3
  epsilon = rnorm(n)
  Y = (X%*%beta) + epsilon   # model1
  #Y = sin(X%*%beta)*exp(X%*%beta) + 0.1*epsilon # model 2
  #Y = (X%*%beta[,1,drop=FALSE])/(1+(X%*%beta[,2,drop=FALSE]+1)^2) + 0.2*epsilon # model 3
  categorical = FALSE
  batch_num = 50
  batch_ind = 1:batch_num #using a small batch
  
  ## pre-specify the cutting points and calculate covariance Cov(E(x|y))
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    H = 5
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
      Mt[,i] = as.numeric(Y == i)
    }  
  }
  nh = colSums(Mt)
  ms = t(as.matrix(1/nh))
  Xt_mean = t(as.matrix(colMeans(Xt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  
  Mx = t(Mt) %*% Xt_tilde
  ms_mat = t(ms)%*%matrix(1,1,p)
  Mx_mean = ms_mat*Mx
  Dt_hat = (1/H)*t(Mx_mean) %*% Mx_mean

  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values[1:k]
  Vt = Dt_eig$vectors[,1:k,drop=FALSE]
  
  #initial
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]

  error_vec = rep(NA,n-batch_num)
  Lamdat_mat = matrix(0,n-batch_num,k)
  
  theta = 1
  g_small = 1#0.01(model 1 model 2)
  gamma = c(sqrt(1/n),sqrt(1/n))#0.0001(model 1 model 2)
  beta_hat = schmidt_orthogonal(sir.lasso$beta)
  L = 15
  tol = 1E-6
  
  
  #error_metric = metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) )  #schmidt_orthogonal(beta_hat)
  ## update D_{t+1} i.e Dt_hat_update
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
      id_update = Y_new_class
    }
    nh[id_update] = nh[id_update] + 1
    ms[1,id_update] = 1/nh[id_update]
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    
    Mx_new = Mx + t(Mt_new) %*% X_new_central
    ms_mat = t(ms)%*%matrix(1,1,p)
    Mx_mean_new = ms_mat*Mx_new
    Dt_hat_update = (1/H)*t(Mx_mean_new)%*%Mx_mean_new
    ## gb
    Vt_update = Vt + (i+1)^(-1/2)*(Dt_hat_update%*%Vt)
    Vt_update = schmidt_orthogonal(Vt_update)
    Lamdat_update = Lamdat + (i+1)^(-1/2)*(diag(t(Vt)%*%Dt_hat_update%*%Vt) - Lamdat)
    #for(j in 1:k){
      #ut_update[,j] = i/(i+1)*ut[,j,drop=FALSE] + (1/((i+1)*Lamdat[j]))*Dt_hat_update%*%ut[,j,drop=FALSE]
      #ut_update[,j] = ut[,j,drop=FALSE] + (i+1)^(-1)*(Dt_hat_update/Lamdat[j] - diag(p))%*%ut[,j,drop=FALSE]
    #  phit = Mx_mean_new%*%Vt[,j,drop=FALSE]
    #  Vt_update[,j] = Vt[,j,drop=FALSE] + (i+1)^(-1)*(Dt_hat_update%*%Vt[,j,drop=FALSE])
    #  Vt_update[,j] = Vt_update[,j,drop=FALSE]/norm(Vt_update[,j,drop=FALSE], "F")
    #  Lamdat_update[j] = Lamdat[j] + (i+1)^(-1)*(t(phit)%*%phit - Lamdat[j])
    #}    
    #Vt_update = schmidt_orthogonal(Vt_update)
    ## online SIR via Truncated Gradient
    if(k == 1){
      Yt_new = (1/(Lamdat_update))*Mt_new%*%Mx_mean_new%*%Vt_update*ms[1,Y_new_class]
    }
    else{
      Yt_new = Mt_new%*%Mx_mean_new%*%Vt_update%*%diag(1/Lamdat_update)*ms[1,Y_new_class]
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta, g_small*L, gamma[jj])
        beta_hat[,jj] = matrix(beta_vec,p,1)
        #zeros_ind = which(beta_vec == 0)
        beta_hat[,jj] = beta_hat[,jj] + gamma[jj]*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
        #beta_hat[zeros_ind,jj] = 0
      }else{
        beta_hat[,jj] = beta_hat[,jj] + gamma[jj]*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
        #if(length(zeros_ind) != 0){
        #  beta_hat[zeros_ind,jj] = 0
        #}
      }
      if( i == n-batch_num){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta, g_small*L, gamma[jj])
        beta_hat[,jj] = matrix(beta_vec,p,1)
      } 
    }
   
    
    # calculate the distance between B and B_hat
    if(k == 1){
      error_vec[i] = metric( beta/norm(beta, "F"), schmidt_orthogonal(beta_hat) )
    }else{
      error_vec[i] = metric( beta, schmidt_orthogonal(beta_hat) )
    }
    # iteration
    Xt_mean = Xt_mean_new
    Mx = Mx_new
    Lamdat = Lamdat_update
    Vt = Vt_update    
    i = i+1
    
  }
  error_mat[s] = min(error_vec,na.rm = TRUE)
  error_tail_mat[s] = tail(na.omit(error_vec),1) 
  #plot(error_vec)
  #min(error_vec,na.rm = TRUE) 
  #tail(na.omit(error_vec),1) 
  pb$tick()
}

# date 2020.7.26
## model 1
# model 1 p = 20 n = 1000 0.014312, 0.01531882 (0.1529s)
# model 1 p = 20 n = 2000 0.009036844, 0.009693519 (1.2110s)

# model 1 p = 70 n = 1000 0.03498363, 0.0361736 (0.6423s)
# model 1 p = 70 n = 2000 0.03036755, 0.03321354 (2.3962s)

# model 1 p = 150 n = 1000 0.1113956, 0.1314095 (0.9487s)
# model 1 p = 150 n = 2000 0.1091682, 0.147147 (3.6992s)

## model 2
# model 2 p = 20 n = 1000 
# model 2 p = 20 n = 2000

# model 2 p = 70 n = 1000 
# model 2 p = 70 n = 2000

# model 2 p = 120 n = 1000 
# model 2 p = 120 n = 2000

## model 3
# model 3 p = 20 n = 1000 
# model 3 p = 20 n = 2000

# model 3 p = 70 n = 1000
# model 3 p = 70 n = 2000

# model 3 p = 120 n = 1000
# model 3 p = 120 n = 2000

time_end = proc.time()
time_start - time_end

mean(error_mat)
mean(error_tail_mat)
mean(error_mat[error_mat!=1])
mean(error_tail_mat[error_tail_mat!=1])