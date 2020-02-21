functions{
  real get_gamma(real alpha, real nu){
    real delta = alpha / sqrt(1 + alpha^2); 
    real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
    real sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
    real gamma = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
    return(gamma);
  }
}
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
  real<lower=4> nu_E;
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
  nu_E ~ uniform(4,40);
 target += log(2) + student_t_lpdf(y | nu_E, mu, omega_E) + student_t_lcdf(y_scaled | nu_E, 0, 1);
}
generated quantities {
  real gamma_E = get_gamma(alpha_E, nu_E);
}
