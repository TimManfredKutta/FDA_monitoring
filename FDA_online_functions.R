#Libraries
library(expm)
library(MASS)


#mean vector after change


#error generation (Brownian motion)
error_process_B <- function(n,nu,lambda,sigma,x0){
  # Generate equidistant points between 0 and 1
  points <- seq(0, 1, length.out = n)
  # Generate standard Brownian motion evaluations on the points
  x <- c(0, cumsum(sqrt(diff(points)) * rnorm(n )))
  return(x)
}

#error generation (ornstein uhlenbeck)
error_process_O <- function(n,nu,lambda,sigma,x0){
  dw  <- rnorm(n, 0, sqrt(1/n))
  dt  <- 1/n
  x <- c(x0)
  for (i in 2:(n+1)) {
    x[i]  <-  x[i-1] + lambda*(nu-x[i-1])*dt + sigma*dw[i-1]
  }
  return(x);
}

#Generation of discretized functions before and after the change
data_m<- function( mn, cp, sigma) {
  error_matrix =  matrix(data =NA,  nrow = d, ncol =(mn+1))
  if(error_selection==1){
  for(i in 1:(mn+1)){
    error_matrix[,i]= error_process_O((d-1),0,1,sigma,0)
  }
  }
  if(error_selection==2){
    for(i in 1:(mn+1)){
      error_matrix[,i]= error_process_B((d-1),0,1,sigma,0)
    }
  }
  if(dependence ==1){
  error_matrix=error_matrix[,1:mn]
  }
  if(dependence ==2){
  error_matrix=(error_matrix[,1:mn]+error_matrix[,2:(mn+1)])/sqrt(2)
  #MA-1 process of errors
  }
  data_matrix =  matrix(data =NA,  nrow = d, ncol = mn)
  #data before cp
  mean_mat  = matrix(data=rep(mean_vector,(mn-cp)),  nrow = d)
  data_matrix[,1:cp] = error_matrix[,1:cp] 
  data_matrix[,(cp+1):mn] = error_matrix[,(cp+1):mn]+mean_mat
  return(data_matrix) 
}


#Autocovariance estimator
Lag <- function(h, data_mat){
  m = length(data_mat[,1])
  data_mat = t(data_mat[1:(m-h),])%*%data_mat[(h+1):(m),]/m
  return(data_mat)
}

#Quantile simulator, with inputs estimated covariance matrix, zeta, gamma
quantile_simulator<- function(lrmatrix, z, g){
  K = 50
  Q = 1000
  stat=1:Q*0
  for(q in 1:Q){
    n_data = mvrnorm(n = K, mu=rep(0,length(lrmatrix[,1])), lrmatrix, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    for(k in 2:K){
      n_data[k,] = n_data[k,] + n_data[(k-1),]
    }
    for(k in 1:K){
      n_data[k,] = n_data[k,]/max((k/K),z)^g
    }
    stat[q] = max(abs(n_data))/sqrt(K)
  }
  quant = quantile(stat,  probs =c(0.95))
  return(quant)
}

#Test decision function; output is 0 (accept) or 1 (reject)
#Inputs m (M), n (N; time horizon after M), cp (change point location)
#The covariance matrix estimate is stabilized by a small diagonal matrix, to 
#avoid instabilities
test_decision<- function(m,n, cp)
{
  Data = t(data_m((m+n), cp, s))
  mean_before = colSums(Data[1:m,])
  stat=1:n*0
  stat[1] = max(abs(1/m*mean_before-Data[(m+1),])/(sqrt(m)*(1+1/m)*(1/(m+1))^gamma))
  for(k in 2:n){
    mean_after = colSums(Data[(m+1):(m+k),])
    stat[k] = max(abs((k/m)*mean_before-colSums(Data[m:(m+k),]))/(sqrt(m)*(1+k/m)*max((k/(k+m)), zeta)^gamma))
  }
  #Cov = t(Data[1:m,]-colmeans(Data[1:m,]))%*%(Data[1:m,]-colmeans(Data[1:m,]))/m
  Data[1:m,]= Data[1:m,]-colmeans(Data[1:m,])
  Cov = sqrtm(t(Lag(0,Data[1:m,])+2*Lag(1,Data[1:m,])+Lag(2,Data[1:m,]))%*%(Lag(0,Data[1:m,])+2*Lag(1,Data[1:m,])+Lag(2,Data[1:m,])))#$B
  quant = quantile_simulator(Cov, zeta, gamma)
  decision = 0
  if(max(stat)>quant){decision=1}
  #print(c(max(stat),quant))
  #print(l)
  return(decision)
}



quantile_simulator<- function(lrmatrix, z, g){
  K = 50
  Q = 1000
  stat=1:Q*0
  for(q in 1:Q){
    n_data = mvrnorm(n = K, mu=rep(0,length(lrmatrix[,1])), lrmatrix, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
    for(k in 2:K){
      n_data[k,] = n_data[k,] + n_data[(k-1),]
    }
    for(k in 1:K){
      n_data[k,] = n_data[k,]/max((k/K),z)^g
    }
    stat[q] = max(abs(n_data))/sqrt(K)
  }
  quant = quantile(stat,  probs =c(0.95))
  return(quant)
}



#Parameters

N = 400       #Total time horizon
d  = 50       #Number of grid points
s = 1         #volatility of error process
zeta = 0.05   #boundary parameter
gamma = 0   #paramter in statistic
error_selection = 1 #Error type, with 1 OU and 2 BM
CP = 100
dependence = 1 
hyp=0 #parameter that quantifies how deep we are in the alternative; hyp=0 is H_0
#hyp is called a in main paper



#One output 
decision = test_decision(M,N,CP)
print(decision)