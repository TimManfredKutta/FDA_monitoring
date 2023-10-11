#Libraries
library(expm)
library(MASS)


#Generation of discretized denstiy estimates before and after the change
data_m<- function( mn, cp) {
  alpha_array  = runif(mn+1, 0, 1)
  beta_array  = runif(mn+1, min = 0.01, max = 1)
  if(dependence ==1){
    alpha_array = alpha_array[1:mn]
    beta_array = beta_array[1:mn]
  }
  if(dependence ==2){
    alpha_array = (alpha_array[1:mn]+alpha_array[2:(mn+1)])/sqrt(2)
    beta_array = beta_array[1:mn]
  }
  data_matrix =  matrix(data =NA,  nrow = d, ncol = mn)
  L=30
  for(i in 1:cp){
    dens_data =  rnorm(L, alpha_array[i],  beta_array[i]  )
    #dens_data =  rnorm(L, alpha_array[i],  beta_array[i])
    dens_f = density(dens_data, bw = "SJ", adjust = 1, kernel = "gaussian",from = -5, to = 5)
    dens_f_v = approx(dens_f$x, dens_f$y, xout = grid_d)$y
    dens_f_v[is.na(dens_f_v)] <- 0
    data_matrix[,i] = dens_f_v
  }
  for(i in (cp+1):mn){
    dens_data =  rnorm(L, alpha_array[i]+hyp,  beta_array[i])
    dens_f = density(dens_data, bw = "SJ", adjust = 1, kernel = "gaussian",from = -5, to = 5)
    dens_f_v = approx(dens_f$x, dens_f$y, xout = grid_d)$y
    dens_f_v[is.na(dens_f_v)] <- 0
    data_matrix[,i] = dens_f_v
  }
  return(data_matrix)
}

#Function for the G-norm
psi <- function(x, w) {
  result = (x< -w)*exp(x+w)+(x> -w)*(x< w) + (x> w)*exp(-(x-w))
  return(result)
}

#G-norm
gnorm<- function(x,b,w){
  x = abs(x)
  result = max(pmin(x,b)*psi(grid_d,w))
  return(result)
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
    stat[q] = gnorm(abs(n_data)/sqrt(K),30,3)
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
  grid_d =  seq(-5,5, length.out = d)
  Data = t(data_m((m+n), cp))
  mean_before = colSums(Data[1:m,])
  stat=1:n*0
  stat[1] = gnorm(abs(1/m*mean_before-Data[(m+1),])/(sqrt(m)*(1+1/m)*max((1/(m+1)),zeta)^gamma),30,3)
  for(k in 2:n){
    mean_after = colSums(Data[(m+1):(m+k),])
    stat[k] = gnorm(abs((k/m)*mean_before-mean_after)/(sqrt(m)*(1+k/m)*max((k/(k+m)), zeta)^gamma),30,3)
    
  }
   Data[1:m,]= Data[1:m,]-colMeans(Data[1:m,])
  Cov = sqrtm(t(Lag(0,Data[1:m,])+2*Lag(1,Data[1:m,])+Lag(2,Data[1:m,]))%*%(Lag(0,Data[1:m,])+2*Lag(1,Data[1:m,])+Lag(2,Data[1:m,]))+diag(0.00000000001,d))#$B
  quant = quantile_simulator(Cov, zeta, gamma)
  decision = 0
  if(max(stat)>quant){decision=1}
  return(decision)
}



#Parameters
N = 500       #Total time horizon
d  = 50       #Number of grid points
grid_d =  seq(-5,5, length.out = d)
zeta = 0.05      #boundary parameter
gamma = 0.3   #paramter in statistic
M = 75
CP = 100
dependence = 1 
hyp=0 #parameter that quantifies how deep we are in the alternative; hyp=0 is H_0
      #hyp is called a in main paper

#One output 
decision = test_decision(M,N,CP)
print(decision)
