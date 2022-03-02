## online SIR via gradient descent
## Cheng Haoyang
## 2020.6.25
set.seed(1234567890)
N = 100
error_mat = array(0, N)
error_tail_mat = array(0, N)
error_F_mat = array(0, N)
error_F_tail_mat = array(0, N)
pb <- progress_bar$new(total = N)  #n-batch_num
time_start = proc.time()
for(s in 1:N){
  n = 1000
  p = 20
  k = 1 # the number of directions
  X = mvrnorm(n, mu = rep(0,p), Sigma = S(0.3,p))
  beta = matrix(0, nrow = p, ncol = k)
  beta[1:2,1] = 1 # model 1
  #beta[c(2,4,6),1] = 1 # model 2
  #beta[1,1] = 1 # model 3
  #beta[2,2] = 1 # model 3
  epsilon = rnorm(n)
  Y = (X%*%beta) + epsilon   # model1
  #Y = sin(X%*%beta)*exp(X%*%beta) + 0.1*epsilon # model 2
  #Y = (X%*%beta[,1,drop=FALSE])/(1+(X%*%beta[,2,drop=FALSE]+1)^2) + 0.2*epsilon
  categorical = FALSE
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
  #sir.lasso <- LassoSIR( X, Y, 10, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, screening=FALSE,no.dim=k)
  
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
  #Bt_hat = sir.lasso$beta
  #Vt_update = Vt[,1:k,drop=FALSE]

  ##
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
  #plot(error_vec)
  #min(error_vec)
  #tail(error_vec,1)
  error_mat[s] = min(error_vec,na.rm = TRUE)
  error_tail_mat[s] = tail(na.omit(error_vec),1) 
  error_F_mat[s] = min(error_vec_F,na.rm = TRUE)
  error_F_tail_mat[s] = tail(na.omit(error_vec_F),1)
  pb$tick()
}
# date 2020.7.26
## model 1
# model 1 p = 20 n = 1000 0.1070209, 0.01535596 (0.0814s)
# model 1 p = 20 n = 2000 0.091319, 0.091319 (0.1625s)
# model 1 p = 20 n = 3000 0.00734914, 0.02474031 (3.8022s)
# model 1 p = 20 n = 5000 0.00588265 (8.9597s)

# model 1 p = 70 n = 1000 0.8756856,0.8756863 (0.4126s)
# model 1 p = 70 n = 2000 0.8759037,0.8759046 (0.8553s)

# model 1 p = 120 n = 1000 0.8963301, 0.896333 (0.5238s)
# model 1 p = 120 n = 2000 0.9227305, 0.9227328 (3.0790s)

## model 2
# model 2 p = 20 n = 1000 0.1345215, 0.1345215 (0.0803s)
# model 2 p = 20 n = 2000 0.1253261, 0.1253261 (0.1635s)

# model 2 p = 70 n = 1000 0.8465268,0.8465282 (0.3992s)
# model 2 p = 70 n = 2000 0.8486009,0.8486019 (0.8343s)

# model 2 p = 120 n = 1000 0.9085826,0.9085855 (1.4406s)
# model 2 p = 120 n = 2000 0.9100247,0.9100274 (3.0697s)

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
mean(error_F_mat)
mean(error_F_tail_mat)