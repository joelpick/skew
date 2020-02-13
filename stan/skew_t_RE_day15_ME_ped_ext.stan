// functions{
//   real get_mu(real omega, real alpha, real nu){
//     real delta = alpha / sqrt(1 + alpha^2); 
//     real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
//     real mu = 0 + omega*b_nu*delta;
//     return(mu);
//   }
//   real get_sigma2(real omega, real alpha, real nu){
//     real delta = alpha / sqrt(1 + alpha^2); 
//     real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
//     real sigma2_z = nu/(nu-2) - (b_nu*delta)^2;
//     real sigma2 = omega^2 * sigma2_z;
//     return(sigma2);
//   }
//   real get_gamma(real alpha, real nu){
//     real delta = alpha / sqrt(1 + alpha^2); 
//     real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
//     real sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
//     real gamma = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
//     return(gamma);
//   }
// }
data {
	int<lower=0> N;                 
	int J;
  int<lower=0> N_nest;
  int<lower=0> N_ind;
  int<lower=0> N_ped;             // number of individuals in pedigree
  int<lower=0> N_NoParents;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_ParentsNoOffspring;

	vector[N] y;                    
  matrix [N, J] X;
  vector[N_ped] MSV;               // Mendelian sampling variance

  int<lower=0, upper=N_nest> nest_id[N];
  int<lower=0, upper=N_ind> ind_id[N];
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
  // fixed effects
  vector[J] beta_tilde;              

  // nest of rearing effects
  vector [N_nest] nest_effects;
	real<lower=0> sigma_nest;
  real alpha_nest;
  real<lower=4> nu_nest;

  // breeding values
  vector [N_NoParents] A_NoParents;  // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring;
  real<lower=0> sigma_A;
  real alpha_A;
  real<lower=4> nu_A;

  // Residual
  vector [N_ind] ind_effects;
  real<lower=0> sigma_ind;
  real alpha_ind;
  real<lower=4> nu_ind;

  // measurement error
  real<lower=0> sigma_E;
}
transformed parameters{

  real delta_nest = alpha_nest / sqrt(1 + alpha_nest^2); 
  real b_nu_delta_nest = delta_nest * sqrt(nu_nest/pi()) * tgamma((nu_nest-1)/2)/tgamma(nu_nest/2) ;
  real omega_nest = sigma_nest/sqrt(nu_nest/(nu_nest-2) - (b_nu_delta_nest)^2);
  real xi_nest = 0 - omega_nest*b_nu_delta_nest;
  vector [N_nest] nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;

  real delta_A = alpha_A / sqrt(1 + alpha_A^2); 
  real b_nu_A = sqrt(nu_A/pi()) * tgamma((nu_A-1)/2)/tgamma(nu_A/2);
  real sigma_z_A = sqrt(nu_A/(nu_A-2) - (b_nu_A*delta_A)^2);
  real omega_A = sigma_A/sigma_z_A;
  real xi_A = 0 - omega_A*b_nu_A*delta_A;

  real delta_ind = alpha_ind / sqrt(1 + alpha_ind^2); 
  real b_nu_delta_ind = delta_ind * sqrt(nu_ind/pi()) * tgamma((nu_ind-1)/2)/tgamma(nu_ind/2) ;
  real omega_ind = sigma_ind/sqrt(nu_ind/(nu_ind-2) - (b_nu_delta_ind)^2);
  real xi_ind = 0 - omega_ind*b_nu_delta_ind;
  vector [N_ind] ind_scaled = (ind_effects-xi_ind)/omega_ind * alpha_ind;


  // vector [N] mu = X*beta + nest_effects[nest_id];

}
model {
	vector [N] mu;
  vector [N_ped] A_scaled;          // scaled breeding value
  vector [N_ped] A;                    // breeding value
  real A_mu;
  int PO_i;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring_mean;
 
  // base population
  target += log(2) + student_t_lpdf(A_NoParents | nu_A, 0, 1) + student_t_lcdf(A_NoParents * alpha_A | nu_A,0, 1);
  A_scaled[NoParents] = A_NoParents;

  // others 
  // those with parents and offspring need to be in a loop as their parents values needs to be looked up and subsequently their value will need to be looked up for their offsprings breeding value
  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
    A_mu = (A_scaled[dam[PO_i]] + A_scaled[sire[PO_i]])*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);
    A_scaled[PO_i] = A_ParentsOffspring[i];
  }
  
  // those with parents but no offspring - can look their parents breeding values up and future individuals values do not depend on them
  A_ParentsNoOffspring_mean = (A_scaled[dam[ParentsNoOffspring]] + A_scaled[sire[ParentsNoOffspring]])*0.5;
  A_ParentsNoOffspring ~ normal( A_ParentsNoOffspring_mean , MSsd[ParentsNoOffspring]);
  A_scaled[ParentsNoOffspring] = A_ParentsNoOffspring ;
  
  A = A_scaled * omega_A + xi_A;

  mu = Q * beta_tilde + A[animal_id] + nest_effects[nest_id] + ind_effects[ind_id];



  beta_tilde ~ normal(0, 10);

  alpha_nest ~ normal(0,10);
	sigma_nest ~ cauchy(0, 10);
  nu_nest ~ uniform(4,40);
  
  alpha_A ~ normal(0,10);
  sigma_A ~ cauchy(0, 10);
  nu_A ~ uniform(4,40);
  
  alpha_ind ~ normal(0,10);
  sigma_ind ~ cauchy(0, 10);
  nu_ind ~ uniform(4,40);
  
  sigma_E ~ cauchy(0, 10);
  

  target += log(2) + student_t_lpdf(nest_effects | nu_nest, 0, omega_nest) + student_t_lcdf((nest_effects)/omega_nest * alpha_nest | nu_nest,0, 1);

  target += log(2) + student_t_lpdf(ind_effects | nu_ind, 0, omega_ind) + student_t_lcdf((ind_effects)/omega_ind * alpha_ind | nu_ind,0, 1);

  y ~ normal(mu,sigma_E);

}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x

  // real mu_nest = get_mu(omega_nest, alpha_nest, nu_nest);
  // real sigma2_nest = get_sigma2(omega_nest, alpha_nest, nu_nest);
  // real gamma_nest = get_gamma(alpha_nest, nu_nest);

  // real mu_A = get_mu(omega_A, alpha_A, nu_A);
  // real sigma2_A = get_sigma2(omega_A, alpha_A, nu_A);
  // real gamma_A = get_gamma(alpha_A, nu_A);
  
  // real mu_ind = get_mu(omega_ind, alpha_ind, nu_ind);
  // real sigma2_ind = get_sigma2(omega_ind, alpha_ind, nu_ind);
  // real gamma_ind = get_gamma(alpha_ind, nu_ind);

  // real sigma2_E = sigma_E^2;
}
