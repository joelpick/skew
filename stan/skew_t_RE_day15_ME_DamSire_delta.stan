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
  int<lower=0> N_dam_sire;
  int<lower=0> N_ind;

	vector[N] y;                    
  matrix [N, J] X;

  int<lower=0, upper=N_nest> nest_id[N];
  int<lower=0, upper=N_dam_sire> dam_id[N]; 
  int<lower=0, upper=N_dam_sire> sire_id[N];
  int<lower=0, upper=N_ind> ind_id[N];

}
transformed data{
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

  vector [N_dam_sire] dam_sire_effects;
  real<lower=0> sigma_dam_sire;
  real<lower=-1,upper=1> delta_dam_sire;
  real<lower=4> nu_dam_sire;

  vector [N_ind] ind_effects;
  real<lower=0> sigma_ind;
  real<lower=-1,upper=1> delta_ind;
  real<lower=4> nu_ind;

  real<lower=0> sigma_E;

}
transformed parameters{

  real alpha_nest = get_alpha(delta_nest);
  real omega_nest = get_omega(sigma_nest,delta_nest,nu_nest);
  real xi_nest = get_xi(omega_nest,delta_nest,nu_nest);

  real alpha_dam_sire = get_alpha(delta_dam_sire);
  real omega_dam_sire = get_omega(sigma_dam_sire,delta_dam_sire,nu_dam_sire);
  real xi_dam_sire = get_xi(omega_dam_sire,delta_dam_sire,nu_dam_sire);

  real alpha_ind = get_alpha(delta_ind);
  real omega_ind = get_omega(sigma_ind,delta_ind,nu_ind);
  real xi_ind = get_xi(omega_ind,delta_ind,nu_ind);

}
model {
  vector [N_nest] nest_scaled = (nest_effects-xi_nest)/omega_nest * alpha_nest;

  vector [N_dam_sire] dam_sire_scaled = (dam_sire_effects-xi_dam_sire)/omega_dam_sire * alpha_dam_sire;

  vector [N_ind] ind_scaled = (ind_effects-xi_ind)/omega_ind * alpha_ind;

  vector [N] mu = Q * beta_tilde + nest_effects[nest_id] + dam_sire_effects[dam_id] + dam_sire_effects[sire_id] + ind_effects[ind_id];

  beta_tilde ~ normal(0,100);

  delta_nest ~ uniform(-1,1);
  sigma_nest ~ cauchy(0,10); 
  nu_nest ~ uniform(4,40);

  delta_dam_sire ~ uniform(-1,1);
  sigma_dam_sire ~ cauchy(0, 10);
  nu_dam_sire ~ uniform(4,40);
  
  delta_ind ~ uniform(-1,1);
  sigma_ind ~ cauchy(0, 10);
  nu_ind ~ uniform(4,40);

  sigma_E ~ cauchy(0, 10);

    
  target += log(2) + student_t_lpdf(dam_sire_effects | nu_dam_sire, xi_dam_sire, omega_dam_sire) + student_t_lcdf(dam_sire_scaled| nu_dam_sire,0, 1);
  target += log(2) + student_t_lpdf(nest_effects | nu_nest, xi_nest, omega_nest) + student_t_lcdf(nest_scaled | nu_nest,0, 1);
  target += log(2) + student_t_lpdf(ind_effects | nu_ind, xi_ind, omega_ind) + student_t_lcdf(ind_scaled | nu_ind,0, 1);
  y ~ normal(mu,sigma_E);
}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x

  real gamma_nest = get_gamma(delta_nest, nu_nest);
  real gamma_dam_sire = get_gamma(delta_dam_sire, nu_dam_sire);
  real gamma_ind = get_gamma(delta_ind, nu_ind);
}

