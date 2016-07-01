// http://ito-hi.blog.so-net.ne.jp/2014-08-28
data {
  int <lower = 0> N;
  matrix [1, N] y;
}
transformed data {
  matrix [4, 1] F;
  vector [4] m0;
  cov_matrix [4] C0;
  
  F [1, 1] <- 1;
  F [2, 1] <- 0;
  F [3, 1] <- 1;
  F [4, 1] <- 0;

  m0 [1] <- 0;
  m0 [2] <- 0;
  m0 [3] <- 0;
  m0 [4] <- 0;
  
  C0 <- diag_matrix(rep_vector(1.0e+7,4));
}
parameters {
  // need to impose stationarity constraints
 real <lower=-1,upper=1> phi2; 
 real <upper=(1 - fabs(phi2))> phi1; 
 
 
  real <lower = 0> sigma ; // for V
  vector<lower = 0>[4] W_diag;
  
}
transformed parameters {
  vector [1] V;
  cov_matrix [4] W;
  matrix [4, 4] G;

  G [1, 1] <- 1;
  G [1, 2] <- 1;
  G [1, 3] <- 0;
  G [1, 4] <- 0;
  G [2, 1] <- 0;
  G [2, 2] <- 1;
  G [2, 3] <- 0;
  G [2, 4] <- 0;  
  G [3, 1] <- 0;
  G [3, 2] <- 0;
  G [3, 3] <- phi1;
  G [3, 4] <- 1;
  G [4, 1] <- 0;
  G [4, 2] <- 0;  
  G [4, 3] <- phi2;
  G [4, 4] <- 0;
  
  
  V [1] <- sigma  * sigma ;
  W <- diag_matrix(W_diag);
  W [4, 4] <- 0;

}
model {
  sigma ~ uniform (0, 1.0e-1);
  phi2 ~ normal(0,2.0/3);
  phi1 ~ normal(0,1.0/3);
  W[2,1] ~ inv_gamma(1, 0.001);
  W[3,1] ~ inv_gamma(1, 0.001);
  W[4,1] ~ inv_gamma(1, 0.001);
  y ~ gaussian_dlm_obs (F, G, V, W, m0, C0);
  
}