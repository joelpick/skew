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
	int J;
  int<lower=0> N_nest;
  int<lower=0> N_ped;             // number of individuals in pedigree
  int<lower=0> N_NoParents;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_ParentsNoOffspring;

	vector[N] y;                    
  matrix [N, J] X;
  vector[N_ped] MSV;               // Mendelian sampling variance

  int<lower=0, upper=N_nest> nest_id[N];
  int<lower=0, upper=N_ped> animal_id[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=N_ped> dam[N_ped];       // dam id in order of animal id
  int<lower=0, upper=N_ped> sire[N_ped];      // sire id  
  int<lower=0, upper=N_ped> NoParents[N_NoParents];
  int<lower=0, upper=N_ped> ParentsOffspring[N_ParentsOffspring];
  int<lower=0, upper=N_ped> ParentsNoOffspring[N_ParentsNoOffspring];
}
transformed data{
  vector[N_ped] MSsd = sqrt(MSV);
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

  vector [N_NoParents] A_NoParents;  // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring;
  real <lower=0> sigma_A;            // sd of A

  real<lower=0> sigma_E;
  real alpha_E;
  real<lower=4> nu_E;
}
transformed parameters{
  real b_nu_delta_nest = (alpha_nest / sqrt(1 + alpha_nest^2)) * sqrt(nu_nest/pi()) * tgamma((nu_nest-1)/2)/tgamma(nu_nest/2) ;
  real omega_nest = sigma_nest/sqrt(nu_nest/(nu_nest-2) - (b_nu_delta_nest)^2);
  real xi_nest = 0 - omega_nest*b_nu_delta_nest;
  vector [N_nest] nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;

  real b_nu_delta_E = (alpha_E / sqrt(1 + alpha_E^2)) * sqrt(nu_E/pi()) * tgamma((nu_E-1)/2)/tgamma(nu_E/2) ;
  real omega_E = sigma_E/sqrt(nu_E/(nu_E-2) - (b_nu_delta_E)^2);
  real xi_E = 0 - omega_E*b_nu_delta_E;

  vector [N_ped] A_scaled;          // scaled breeding value
  vector [N] mu;
  vector [N] y_scaled;
  
  A_scaled[1] = 0;
  A_scaled[NoParents] = A_NoParents;
  A_scaled[ParentsNoOffspring] = A_ParentsNoOffspring ;
  A_scaled[ParentsOffspring] = A_ParentsOffspring;

  mu = Q * beta_tilde + A_scaled[animal_id] * sigma_A + nest_effects[nest_id] + xi_E;

  y_scaled = (y - mu)/omega_E * alpha_E;
}
model {

  real A_mu;
  int PO_i;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring_mean;

  A_NoParents ~ normal( 0 , 1 );
  
  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
    A_mu = (A_scaled[dam[PO_i]] + A_scaled[sire[PO_i]])*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);
  }
  
  A_ParentsNoOffspring_mean = (A_scaled[dam[ParentsNoOffspring]] + A_scaled[sire[ParentsNoOffspring]])*0.5;
  A_ParentsNoOffspring ~ normal( A_ParentsNoOffspring_mean , MSsd[ParentsNoOffspring]);
  
  beta_tilde ~ normal(0,10);

  alpha_nest ~ normal(0,10);
  sigma_nest ~ cauchy(0,10);
  nu_nest ~ uniform(4,40);

  sigma_A ~ cauchy(0,10);

  alpha_E ~ normal(0,10);
  sigma_E ~ cauchy(0,10);
  nu_E ~ uniform(4,40);

  target += log(2) + student_t_lpdf(nest_effects | nu_nest, xi_nest, omega_nest) + student_t_lcdf(nest_scaled | nu_nest,0, 1);
  target += log(2) + student_t_lpdf(y | nu_E, mu, omega_E) + student_t_lcdf(y_scaled| nu_E,0, 1);
}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x
  real gamma_nest = get_gamma(alpha_nest, nu_nest);
  real gamma_E = get_gamma(alpha_E, nu_E);
}
