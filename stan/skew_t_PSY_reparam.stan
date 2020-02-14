data {
	int<lower=0> N;                 
	vector[N] y;                    
	int K;
	matrix [N, K] X ;              
}
parameters {
	
	vector [K]beta;
  real<lower=0> sigma_E;
	real alpha_E;
  real<lower=3> nu_E;
}
transformed parameters{

  real delta_E = alpha_E / sqrt(1 + alpha_E^2); 
  real b_nu_E = sqrt(nu_E/pi()) * tgamma((nu_E-1)/2)/tgamma(nu_E/2);
  real sigma_z_E = sqrt(nu_E/(nu_E-2) - (b_nu_E*delta_E)^2);
  real omega_E = sigma_E/sigma_z_E;
  real xi_E = 0 - omega_E*b_nu_E*delta_E;

  vector [N] mu = xi_E + X*beta;
  vector [N] y_scaled = (y - mu )/omega_E * alpha_E;

}
model {
  beta ~ normal(0, 100);
  sigma_E ~ cauchy(0, 10);
	alpha_E ~ normal(0,10);
  nu_E ~ uniform(3,30);
 target += log(2) + student_t_lpdf(y | nu_E, mu, omega_E) + student_t_lcdf(y_scaled | nu_E, 0, 1);
}
