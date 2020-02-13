data {
	int<lower=0> N;                 
	vector[N] y;                    
	int K;
	matrix [N, K] X;
  int<lower=0> N_nest;
  int<lower=0, upper=N_nest> nest_id[N];         
}
parameters {
	
	vector [K]beta;
  //real xi;

  vector [N_nest] nest_effects;
	real<lower=0> sigma_nest;
  real alpha_nest;
  real<lower=3> nu_nest;

  real<lower=0> sigma_E;
  real alpha_E;
  real<lower=3> nu_E;
}
transformed parameters{
  real delta_nest = alpha_nest / sqrt(1 + alpha_nest^2); 
  real b_nu_nest = sqrt(nu_nest/pi()) * tgamma((nu_nest-1)/2)/tgamma(nu_nest/2);
  real sigma_z_nest = sqrt(nu_nest/(nu_nest-2) - (b_nu_nest*delta_nest)^2);
  real omega_nest = sigma_nest/sigma_z_nest;
  real xi_nest = 0 - omega_nest*b_nu_nest*delta_nest;
  vector [N_nest] nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;

  vector [N] mu = X*beta + nest_effects[nest_id];

  real delta_E = alpha_E / sqrt(1 + alpha_E^2); 
  real b_nu_E = sqrt(nu_E/pi()) * tgamma((nu_E-1)/2)/tgamma(nu_E/2);
  real sigma_z_E = sqrt(nu_E/(nu_E-2) - (b_nu_E*delta_E)^2);
  real omega_E = sigma_E/sigma_z_E;
  real xi_E = 0 - omega_E*b_nu_E*delta_E;

  vector [N] y_scaled = (y - mu - xi_E)/omega_E * alpha_E;

}
model {

  beta ~ normal(0, 10);

  alpha_nest ~ normal(0,10);
  sigma_nest ~ cauchy(0, 10);
  nu_nest ~ uniform(3,30);

  alpha_E ~ normal(0,10);
  sigma_E ~ cauchy(0, 10);
  nu_E ~ uniform(3,30);

  target += log(2) + student_t_lpdf(nest_effects | nu_nest, xi_nest, omega_nest) + student_t_lcdf(nest_scaled | nu_nest,0, 1);

  target += log(2) + student_t_lpdf(y | nu_E, mu+xi_E, omega_E) + student_t_lcdf(y_scaled| nu_E,0, 1);

}
