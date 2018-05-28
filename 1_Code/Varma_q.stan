data {
  int T;  
  int NS;
  int p;
  vector[NS] S_obs[T]; // 
  vector[NS] H[T]; //
}
parameters {
  vector<lower= -10,upper= 10>[NS] S_0[p]; // 
  vector<lower= -10,upper= 10>[NS] S_hat_post[T-p]; // 
  vector<lower= -5,upper= 5>[3] beta[NS]; //
  vector<lower= -5,upper= 5>[3] beta_mu; //
  vector<lower=0.00001>[NS] L_sigma;  // 
  cholesky_factor_corr[NS] L_Omega; // 
  vector<lower=0.00001>[3] L_sigma_beta;  // 
  cholesky_factor_corr[3] L_Omega_beta; // 
  real<lower=0,upper= 1> z; //
  real q;
}
transformed parameters { 
  vector[NS] S_hat[T]; // 
  vector[NS] S_hat_pre[T-p];
  vector[NS] S_mu[T-p]; //
  matrix[NS,NS] L_Sigma;
  matrix[3,3] L_Sigma_beta;
  vector[NS] A1;
  vector[NS] A2;
  vector[NS] R;
  for (i in 1:p){
      S_hat[i] =  S_0[i];
  }
  for(i in 1:NS){
    A1[i] = exp(beta[i,1])./(exp(beta[i,1])+1);
    A2[i] = exp(beta[i,2])./(exp(beta[i,2])+1);
    R[i] =beta[i,3];
  }
  for (t in 1:(T-p)){
    S_hat[t + p] = S_hat_post[t];
    S_hat_pre[t] = log(exp(S_hat_post[t])+H[t+p]);
  }
  // process equation
  S_mu[1]  = log(diag_matrix(A1)*exp(S_hat[p])+
      exp(diag_matrix(A2)*S_hat[1]+R)); //
  for (t in 2:(T-p)) {
    S_mu[t]  = log(diag_matrix(A1)*exp(S_hat[t+p-1])+
      exp(diag_matrix(A2)*S_hat[t+p-3]+R)); //
  }
  L_Sigma = diag_pre_multiply(L_sigma,L_Omega); // 
  L_Sigma_beta = diag_pre_multiply(L_sigma_beta,L_Omega_beta); // 
}
model {
  for(t in 1:T){
    S_obs[t]~lognormal(S_hat[t]+q,z); //
  }
   S_hat_pre~multi_normal_cholesky(S_mu,L_Sigma); //
  beta_mu ~normal(0,1.5);//
  L_Omega ~ lkj_corr_cholesky(2.0); //
  L_sigma ~ cauchy(0, 2.5); //
  for (i in 1:p){
    S_0[p]~normal(0,1);//
  }
  
  // adjust probability accumulator using log-determinant of transform
  for(t in 1:(T-p)){
    target+=S_hat_post[t]-log(exp(S_hat_post[t])+H[t+p]);
  }
  L_sigma_beta ~ cauchy(0,2.5); // vague priors for the ranef SD intercept
  L_Omega_beta ~ lkj_corr_cholesky(2); // vague priors for the ranef corr matrix
  beta ~ multi_normal_cholesky(beta_mu,L_Sigma_beta); // site ranef
  z~cauchy(0,2.5); // error variance parameters
  q ~ normal(0,0.05);
}