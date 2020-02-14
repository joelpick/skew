
data {
	int<lower=0> N;                 
	int J;
  int<lower=0> N_nest;

	vector[N] y;                    
  matrix [N, J] X;

  int<lower=0, upper=N_nest> nest_id[N];

}
transformed data{
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);
}

parameters {

  vector[J] beta_tilde;              // fixed effects
  //real xi;

  vector [N_nest] nest_effects;
	real<lower=0> sigma_nest;
  real alpha_nest;
  real<lower=4> nu_nest;

  real<lower=0> sigma_E;
}
transformed parameters{

  real b_nu_delta_nest = (alpha_nest / sqrt(1 + alpha_nest^2)) * sqrt(nu_nest/pi()) * tgamma((nu_nest-1)/2)/tgamma(nu_nest/2) ;
  real omega_nest = sigma_nest/sqrt(nu_nest/(nu_nest-2) - (b_nu_delta_nest)^2);
  real xi_nest = 0 - omega_nest*b_nu_delta_nest;
  vector [N_nest] nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;
 
  vector [N] mu = Q * beta_tilde + nest_effects[nest_id];
}
model {
  beta_tilde ~ normal(0, 100);

  alpha_nest ~ normal(0,10);
	sigma_nest ~ cauchy(0, 10);
  nu_nest ~ uniform(4,40);
  
  sigma_E ~ cauchy(0, 10);
  
  target += log(2) + student_t_lpdf(nest_effects | nu_nest, xi_nest, omega_nest) + student_t_lcdf(nest_scaled | nu_nest,0, 1);
  y ~ normal(mu,sigma_E);

}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x
}

