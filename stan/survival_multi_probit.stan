data{
  int<lower=1> N;
  int<lower=1> J;
  int<lower=2> K;
  int<lower=1> N_nest;

  int<lower=0, upper=1> Y[N];
  matrix[N,J] X;
  int<lower=1, upper=N_nest> nest_id[N];
  int<lower=1, upper=K> age[N];
}
transformed data {
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);
}
parameters {
  vector[J] beta_tilde;

  vector<lower=0>[K] sigma_nest;
  cholesky_factor_corr[K] L_corr_nest;
  matrix[N_nest,K] nest_scaled;

}
model {
  matrix[N_nest, K] nest_effect = (nest_scaled * (diag_pre_multiply(sigma_nest, L_corr_nest))');
  vector[N] mu;
  vector[N] nest_effect_vec;
 
  for (i in 1:N )
    nest_effect_vec[i] = nest_effect[nest_id[i],age[i]];
  
  mu = Phi(Q * beta_tilde + nest_effect_vec);
  // print("mu =", mu);

  beta_tilde ~ normal(0,10);
  to_vector(nest_scaled) ~ normal(0, 1);
  L_corr_nest ~ lkj_corr_cholesky(2);
  sigma_nest ~ cauchy(0, 2.5);
  
  // print("log density before =", target());
  Y ~ bernoulli(mu);
  // print("log density after =", target());

}
generated quantities {
  vector[J] beta = R_inv * beta_tilde;
  matrix [K,K] Sigma_nest = quad_form_diag(multiply_lower_tri_self_transpose(L_corr_nest), sigma_nest);
}

