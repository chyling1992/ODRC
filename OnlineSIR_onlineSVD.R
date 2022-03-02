## online SIR via online SVD
## Cheng Haoyang
## 2020.6.25
set.seed(1234567890)
N = 100
error_mat = array(0, N)
error_tail_mat = array(0, N)
pb <- progress_bar$new(total = N)  #n-batch_num
time_start = proc.time()
for(s in 1:N){
  n = 1000
  p = 20
  X = mvrnorm(n, mu = rep(0,p), Sigma = diag(p))
  beta = matrix(0, nrow = p, ncol = 1)
  beta[1:2,1] = 1 # model 1
  #beta[c(2,4,6),1] = 1 # model 2
  #beta[1,1] = 1 # model 3
  #beta[2,2] = 1 # model 3
  epsilon = rnorm(n)
  Y = (X%*%beta) + epsilon   # model 1
  #Y = sin(X%*%beta)*exp(X%*%beta) + 0.1*epsilon # model 2
  #Y = (X%*%beta[,1,drop=FALSE])/(1+(X%*%beta[,2,drop=FALSE]+1)^2) + 0.2*epsilon # model 3
  categorical = FALSE
  k = 1 # the number of directions
  batch_num = 100
  batch_ind = 1:batch_num
  ## pre-specify the cutting points and calculate covariance Cov(E(x|y))
  Xt = X[batch_ind,]
  Yt = Y[batch_ind]
  if ( categorical == FALSE ){
    cutting_points = as.numeric( quantile(Yt, prob = c(0.2,0.4,0.6,0.8)) )
    slices_interval = c(-Inf,cutting_points,Inf)
    H = 5
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
      Yt_tilde[,i] = as.numeric(Y == i)
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
  
  ##
  error_vec = rep(NA,n-batch_num)
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
      #Lamdat_mat[i,j] = Lamdat_update[j]
    }
    beta_hat = Vt_update
    if(k == 1){
      error_vec[i] = metric( beta/norm(beta, "F"), schmidt_orthogonal(beta_hat) )
    }else{
      error_vec[i] = metric( beta, schmidt_orthogonal(beta_hat) )
    }
    
    # iteration
    At_inverse = At1_inverse
    XYt = XYt1
    Mt_sum = Mt_sum + Mt1_hat
    Gamma_hat_t = Mt_sum/(i+1)
    Lamdat = Lamdat_update
    Vt = Vt_update
    i = i+1
  }
  #plot(error_vec)
  #min(error_vec)
  #tail(error_vec,1)
  error_mat[s] = min(error_vec,na.rm = TRUE)
  error_tail_mat[s] = tail(na.omit(error_vec),1) 
  pb$tick()
}

# date 2020.7.26
## model 1
# model 1 p = 20 n = 1000 0.0215139,0.1120043(1000,0.3464s)
# model 1 p = 20 n = 2000 0.01307186, 0.09156421 (0.7082s)

# model 1 p = 70 n = 1000 0.2663033, 0.4799197(2.1506s) 
# model 1 p = 70 n = 2000 0.3011651, 0.6157155(4.6005s)

# model 1 p = 120 n = 1000 0.4959589, 0.8547244(7.3324s)
# model 1 p = 120 n = 2000 0.4991992, 0.8729821(15.25.59s)

## model 2
# model 2 p = 20 n = 1000 0.02565309,0.1050349 (0.3516s)
# model 2 p = 20 n = 2000 0.02548509,0.1455135 (0.7212s)

# model 2 p = 70 n = 1000 0.3232207, 0.623524 (2.1570s)
# model 2 p = 70 n = 2000 0.3076504, 0.5875797 (4.5685s)

# model 2 p = 120 n = 1000 0.5023125, 0.8210043 (7.0683s)
# model 2 p = 120 n = 2000 0.567575, 0.9028283 (15.0725)

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

#model 1 p=20 0.0215139,0.1120043(1000,0.3464s)
#model 2 p=20 0.02271201,0.09956394(1000,0.3661s)
#model 2 p=80 0.3963056,0.8017357(1000, 3.0633s)
#model 3 p=10 0.1527752,0.3718105(1000, 0.3711s) 0.1322651,0.3430155(10000,34.811s)