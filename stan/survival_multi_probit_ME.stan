data{
  int<lower=1> N;
  int<lower=1> J;
  int<lower=2> K;
  int<lower=1> N_nest;
  int<lower=1> N_x1; // number of observations in repeat measurements
  int<lower=1> N_ind; // number of individuals

  int<lower=0, upper=1> Y[N];
  vector[N_x1] x1; // trait with measurement error

  matrix[N,J] X; // design matrix for other predictors
  matrix[N,K] X1; // design matrix to link x1 at different ages to y

  int<lower=1, upper=N_nest> nest_id[N];
  int<lower=1, upper=K> age[N];

  int<lower=1, upper=N_ind> x1_ind_id[N_x1];
  int<lower=1, upper=N_ind> ind_id[N];

}
transformed data {
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);
}
parameters {
  vector[J] beta_tilde;

  vector<lower=1e-16>[K] sigma_nest;
  cholesky_factor_corr[K] L_corr_nest;
  matrix[N_nest,K] nest_scaled;
  
  real beta0_x1; 
  vector[N_ind] x1_hat_scaled;
  real<lower=0> sigma_E_x1;
  real<lower=0> sigma_ind_x1;

  // slopes of x1_hat and x1_hat^2 on survival
  vector[K] beta_x1_hat;
  vector[K] beta_x1_hat2 ;
}
model {
  vector[N] mu;
  vector[N] nest_effect_vec;
  matrix[N, K] X1_hat;
  
  matrix[N_nest, K] nest_effect = (nest_scaled * (diag_pre_multiply(sigma_nest, L_corr_nest))');

  // repeatability model for x1 
  vector[N_ind] x1_hat = sigma_ind_x1 * x1_hat_scaled;
  
  vector[N_x1] mu_x1  =  beta0_x1 + x1_hat[x1_ind_id];
  
  x1 ~ normal(mu_x1, sigma_E_x1);

  
  for (i in 1:N ){
    nest_effect_vec[i] = nest_effect[nest_id[i],age[i]];
    X1_hat[i,] = x1_hat[ind_id[i]] * X1[i,];  
  }

  

  mu = Phi(Q * beta_tilde 
    + X1_hat  * beta_x1_hat
    + (X1_hat .* X1_hat) * beta_x1_hat2 
    + nest_effect_vec);

  // print("mu =", mu);

  beta_tilde ~ normal(0,100);
  
  beta_x1_hat ~ normal(0,100);
  beta_x1_hat2  ~ normal(0,100);
  beta0_x1 ~ normal(0,100); 

  to_vector(nest_scaled) ~ normal(0, 1);
  x1_hat_scaled ~ normal(0, 1);

  L_corr_nest ~ lkj_corr_cholesky(2);
  sigma_nest ~ cauchy(0, 2.5);
  sigma_ind_x1 ~ cauchy(0, 2.5);
  sigma_E_x1 ~ cauchy(0, 2.5);

  // print("log density before =", target());

  Y ~ bernoulli(mu);
  // print("log density after =", target());

}
generated quantities {
  vector[J] beta = R_inv * beta_tilde;
  matrix [K,K] Sigma_nest = quad_form_diag(multiply_lower_tri_self_transpose(L_corr_nest), sigma_nest);
}

