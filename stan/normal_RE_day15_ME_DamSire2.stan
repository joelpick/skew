
data {
	int<lower=0> N;                 
	int J;
  int<lower=0> N_nest;
  int<lower=0> N_ind;
  int<lower=0> N_dam_sire;

	vector[N] y;                    
  matrix [N, J] X;

  int<lower=0, upper=N_nest> nest_id[N];
  int<lower=0, upper=N_ind> ind_id[N];
  int<lower=0, upper=N_dam_sire> dam_id[N]; 
  int<lower=0, upper=N_dam_sire> sire_id[N];
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

  vector [N_nest] nest_scaled;
	real<lower=0> sigma_nest;

  vector [N_dam_sire] dam_sire_scaled;
  real<lower=0> sigma_dam_sire;

  vector [N_ind] ind_scaled;
  real<lower=0> sigma_ind;

  real<lower=0> sigma_E;
}
transformed parameters{

  vector [N_nest] nest_effects = nest_scaled * sigma_nest;
  vector [N_dam_sire] dam_sire_effects = dam_sire_scaled * sigma_dam_sire;
  vector [N_ind] ind_effects = ind_scaled * sigma_ind;
  
  vector [N] mu = Q * beta_tilde + nest_effects[nest_id] + dam_sire_effects[dam_id] + dam_sire_effects[sire_id] + ind_effects[ind_id];
}
model {
  beta_tilde ~ normal(0, 10);

	sigma_nest ~ cauchy(0, 10);
  sigma_dam_sire ~ cauchy(0, 10);
  sigma_ind ~ cauchy(0, 10);  
  sigma_E ~ cauchy(0, 10);
  
  dam_sire_scaled ~ normal(0,1);
  nest_scaled ~ normal(0,1);
  ind_scaled ~ normal(0,1);

  y ~ normal(mu,sigma_E);

}
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x
}

