functions{
  real get_gamma(real delta, real nu){
    real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
    real sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
    real gamma = (b_nu*delta)/sigma_z^3 * ( (nu*(3-delta^2))/(nu-3) - (3*nu)/(nu-2)  + 2*(b_nu*delta)^2 );
    return(gamma);
  }
  real get_omega(real sigma, real delta, real nu){
    real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
    real sigma_z = sqrt(nu/(nu-2) - (b_nu*delta)^2);
    real omega = sigma/sigma_z;
    return(omega);
  }  
  real get_xi(real omega, real delta, real nu){
    real b_nu = sqrt(nu/pi()) * tgamma((nu-1)/2)/tgamma(nu/2);
    real xi = 0 - omega * b_nu * delta;
    return(xi);
  }
  real get_alpha(real delta){
    real alpha = delta/(sqrt(1-delta^2));
    return(alpha);
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

  vector [N_nest] nest_effects;
	real<lower=0> sigma_nest;
  real<lower=-1,upper=1> delta_nest;
  real<lower=4> nu_nest;

  vector [N_NoParents] A_NoParents;  // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring;
  real <lower=0> sigma_A;            // sd of A

  real<lower=0> sigma_E;
  real<lower=-1,upper=1> delta_E;
  real<lower=4> nu_E;
}
transformed parameters{

  real alpha_nest = get_alpha(delta_nest);
  real omega_nest = get_omega(sigma_nest,delta_nest,nu_nest);
  real xi_nest = get_xi(omega_nest,delta_nest,nu_nest);

  real alpha_E = get_alpha(delta_E);
  real omega_E = get_omega(sigma_E,delta_E,nu_E);
  real xi_E = get_xi(omega_E,delta_E,nu_E);

}
model {

  real A_mu;
  int PO_i;
  vector [N_ParentsNoOffspring] A_ParentsNoOffspring_mean;

  vector [N_ped] A_scaled;          // scaled breeding value
  vector [N] mu;
  vector [N] y_scaled;
  vector [N_nest] nest_scaled;

  A_scaled[1] = 0;
  A_NoParents ~ normal( 0 , 1 );
  A_scaled[NoParents] = A_NoParents;
  
  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
    A_mu = (A_scaled[dam[PO_i]] + A_scaled[sire[PO_i]])*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);
    A_scaled[PO_i] = A_ParentsOffspring[i];
  }

  A_ParentsNoOffspring_mean = (A_scaled[dam[ParentsNoOffspring]] + A_scaled[sire[ParentsNoOffspring]])*0.5;
  A_ParentsNoOffspring ~ normal( A_ParentsNoOffspring_mean , MSsd[ParentsNoOffspring]);
    A_scaled[ParentsNoOffspring] = A_ParentsNoOffspring ;


  nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;

  mu = Q * beta_tilde + A_scaled[animal_id] * sigma_A + nest_effects[nest_id] + xi_E;

  y_scaled = (y - mu)/omega_E * alpha_E;

  beta_tilde ~ normal(0,100);

  delta_nest ~ uniform(-1,1);
  sigma_nest ~ cauchy(0,10); 
  nu_nest ~ uniform(4,40);

  sigma_A ~ cauchy(0,10);

  delta_E ~ uniform(-1,1);
  sigma_E ~ cauchy(0,10);
  nu_E ~ uniform(4,40);

  target += log(2) + student_t_lpdf(nest_effects | nu_nest, xi_nest, omega_nest) + student_t_lcdf(nest_scaled | nu_nest,0, 1);
  target += log(2) + student_t_lpdf(y | nu_E, mu, omega_E) + student_t_lcdf(y_scaled| nu_E,0, 1);
}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x
  real gamma_nest = get_gamma(delta_nest, nu_nest);
  real gamma_E = get_gamma(delta_E, nu_E);
}
