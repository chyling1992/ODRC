## online SIR via Truncated Gradient
## Cheng Haoyang
## 2020.6.20

set.seed(12345)
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
  X = mvrnorm(n, mu = rep(0,p), Sigma = diag(p))
  beta = matrix(0, nrow = p, ncol = 1)
  #beta[1:2,1] = 1 # model 1
  beta[c(2,4,6),1] = 1 # model 2
  #beta[1,1] = 1 # model 3
  #beta[2,2] = 1 # model 3
  epsilon = rnorm(n)
  #Y = (X%*%beta) + epsilon   # model 1
  Y = sin(X%*%beta)*exp(X%*%beta) + 0.1*epsilon # model 2
  #Y = (X%*%beta[,1,drop=FALSE])/(1+(X%*%beta[,2,drop=FALSE]+1)^2) + 0.2*epsilon # model 3
  categorical = FALSE
  batch_num = 100
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
  #ms = t(as.matrix(1/colSums(Mt)))
  Xt_tilde = Xt - as.matrix(rep(1,batch_num)) %*% Xt_mean
  #Sigma_xx = t(Xt_tilde)%*%Xt_tilde
  Mx = t(Mt) %*% Xt_tilde
  ms_mat = t(ms)%*%matrix(1,1,p)
  Mx_mean = ms_mat*Mx
  Dt_hat = (1/H)*t(Mx_mean) %*% Mx_mean
  #Dt_hat_1 = (1/H)*t(Xt_tilde) %*% (Mt*(matrix(1,batch_num,1)%*%ms)) %*% t(Mt*(matrix(1,batch_num,1)%*%ms)) %*% Xt_tilde
  #Dt_hat = (1/H)*t(Xt_tilde) %*% (Mt) %*% t(Mt) %*% Xt_tilde
  
  Dt_eig = eigen(Dt_hat)
  Lamdat = Dt_eig$values
  Vt = Dt_eig$vectors
  
  #initial
  sir.lasso <- LassoSIR( Xt, Yt, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  Dt_sum = matrix(0, nrow = p, ncol = p)
  Dt_sum = Dt_sum + Dt_hat
  Gamma_hat_t = Dt_sum
  Lamdat_update = Lamdat[1:k]
  Vt_update = Vt[,1:k,drop=FALSE]
  
  error_vec = rep(NA,n-batch_num)
  Lamdat_mat = matrix(0,n-batch_num,k)
  
  theta = Inf
  g = 0
  g_small = 0.04#0.01(model 1 model 2)
  gamma = c(0.0001,0.0001)#0.0001(model 1 model 2)
  beta_hat = sir.lasso$beta
  L = 10
  tol = 1E-6
  
  # metric( schmidt_orthogonal(beta), schmidt_orthogonal(sir.lasso$beta) ) 
  #error_metric = metric( beta, schmidt_orthogonal(sir.lasso$beta) )  #schmidt_orthogonal(beta_hat)
  ## update D_{t+1} i.e Dt_hat_update
  i = 1
  while (i <= n-batch_num){
    t = batch_num+i
    X_new = X[t,,drop=FALSE]
    Y_new = Y[t]
    #Xt = rbind(Xt, X_new)
    #Yt = rbind(Yt, Y_new)
    if (categorical == FALSE){
      Y_new_class = findInterval(Y_new,slices_interval, left.open = TRUE)
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new_class] = 1
      id_update = Y_new_class
      #Mt = rbind(Mt, Mt_new)
    }else{
      Mt_new = matrix(0, nrow = 1, ncol = H)
      Mt_new[1, Y_new] = 1
      id_update = Y_new_class
      #Mt = rbind(Mt, Mt_new)
    }
    nh[id_update] = nh[id_update] + 1
    ms[1,id_update] = 1/nh[id_update]
    #ms = t(as.matrix(1/colSums(Mt)))
    #Xt_mean = t(as.matrix(colMeans(Xt)))
    Xt_mean_new = ((t-1)*Xt_mean + X_new)*(1/t)
    X_new_central = X_new - Xt_mean_new
    #Xt_tilde = Xt - as.matrix(rep(1,t)) %*% Xt_mean
    
    
    Mx_new = Mx + t(Mt_new) %*% X_new_central
    ms_mat = t(ms)%*%matrix(1,1,p)
    Mx_mean_new = ms_mat*Mx_new
    Dt_hat_update = (1/H)*t(Mx_mean_new)%*%Mx_mean_new
    #Sigma_xx_update = Sigma_xx + t(X_new_central)%*%X_new_central
    #Dt_hat_update_1 = (1/(H))*t(Xt_tilde) %*% (Mt*(matrix(1,t,1)%*%ms)) %*% t(Mt*(matrix(1,t,1)%*%ms)) %*% Xt_tilde
    #Dt_hat_update = (1/(H))*t(Xt_tilde) %*% (Mt) %*% t(Mt) %*% Xt_tilde

    ## online SVD
    #Dt_eig_update = eigen(Dt_hat_update)
    #Lamdat_update = Dt_eig_update$values[1:k]
    #Vt_update = Dt_eig_update$vectors[,1:k,drop = FALSE]
    
    temp_mat_1 = Gamma_hat_t - Dt_hat_update
    for(j in 1:k){
      Lamdat_update[j] = Lamdat[j] - (i+1)^(-1)*t(Vt[,j,drop=FALSE])%*%temp_mat_1%*%Vt[,j,drop=FALSE]
      Vt_update[,j] = Vt[,j,drop=FALSE] - (i+1)^(-1)*ginv(Lamdat[j]*diag(p) - Gamma_hat_t)%*%temp_mat_1%*%Vt[,j,drop=FALSE]
      Vt_update[,j] = Vt_update[,j,drop=FALSE]/norm(Vt_update[,j,drop=FALSE], "F")
      #Lamdat_mat[i,j] = Lamdat_update[j]
    }
    
    ## online SIR via Truncated Gradient
    if(k == 1){
      #Yt_new = (1/(Lamdat_update))*Mt_new%*%t(Mt*(matrix(1,t,1)%*%ms))%*%Xt_tilde%*%Vt_update*ms[1,Y_new_class]
      #Yt_new = (1/(t*H))*(1/(Lamdat_update))*Mt_new%*%t(Mt)%*%Xt_tilde%*%Vt_update
      Yt_new = (1/(Lamdat_update))*Mt_new%*%Mx_mean_new%*%Vt_update*ms[1,Y_new_class]
    }
    else{
      #Yt_new = (ms[1,Y_new_class])*Mt_new%*%t(Mt*(matrix(1,t,1)%*%ms))%*%Xt_tilde%*%Vt_update%*%diag(1/Lamdat_update)
      #Yt_new = (1/(t*H))*Mt_new%*%t(Mt)%*%Xt_tilde%*%Vt_update%*%diag(1/Lamdat_update)
      Yt_new = Mt_new%*%Mx_mean_new%*%Vt_update%*%diag(1/Lamdat_update)*ms[1,Y_new_class]
    }
    X_stand = (X_new - Xt_mean_new)
    Yt_hat = X_stand%*%beta_hat
    for(jj in 1:k){
      if ( i %% L == 0){
        beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta, g_small*L, gamma[jj])
        beta_hat[,jj] = matrix(beta_vec,p,1)
      #  zeros_ind = which(beta_vec == 0)
        beta_hat[,jj] = beta_hat[,jj] + 2*0.0001*(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      #  beta_hat[zeros_ind,jj] = 0
      }else{
        beta_hat[,jj] = beta_hat[,jj] + 2*0.0001*(1/i)*t(X_stand)%*%(Yt_new[,jj] - Yt_hat[,jj])
      #  if(length(zeros_ind) != 0){
      #    beta_hat[zeros_ind,jj] = 0
      #  }
      }
    }
    if( i == n-batch_num){
      beta_vec = sapply(beta_hat[,jj,drop=FALSE], soft_threshold, theta, g_small*L, gamma[jj])
      beta_hat[,jj] = matrix(beta_vec,p,1)
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
    Dt_sum = Dt_sum + Dt_hat_update
    Gamma_hat_t = Dt_sum/(i+1)
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
# date 2020.7.27
## model 1
# model 1 p = 20 n = 1000 0.01443034, 0.01542314 (0.6283s)
# model 1 p = 20 n = 2000 0.01269831(1.7183s)

# model 1 p = 70 n = 1000 0.03520114, 0.03636093 (2.4122s)
# model 1 p = 70 n = 2000 0.0280802,0.05837671 (6.2289s)

# model 1 p = 120 n = 1000 0.1123718,0.1321764(6.9438s)
# model 1 p = 120 n = 2000 0.09823827, 0.2461622(16.5928s)

## model 2
# model 2 p = 20 n = 1000 0.02888666,0.03238016 (0.4114s)
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
mean(error_mat[error_mat < 0.5])
mean(error_tail_mat[error_tail_mat < 0.5])