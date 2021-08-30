
data {
  int<lower=0> N;                     // number of observations
  int<lower=1> J; // number of predictors
  int<lower=0> N_nest;
  int<lower=0> N_ind;
  int<lower=0> N_ped;                     // number of individuals in pedigree
  int<lower=0> N_NoParentsOffspring;
  int<lower=0> N_ParentsOffspring;
  int<lower=0> N_NoOffspring;

  vector[N] y;                        // response  
  matrix[N, J] X; // matrix of predictors
  vector[N_ped] MSV;                      // Mendelian sampling variance

  int<lower=0, upper=N_nest> nest_id[N];
  int<lower=0, upper=N_ind> ind_id[N];
  int<lower=0, upper=N_ped> animal_id[N];     // animal id relating y to animal in pedigree (refers to row number of dam and sire)
  int<lower=0, upper=N_ped> dam[N_ped];       // dam id in order of animal id
  int<lower=0, upper=N_ped> sire[N_ped];      // sire id  
  int<lower=0, upper=N_ped> NoParentsOffspring[N_NoParentsOffspring];
  int<lower=0, upper=N_ped> ParentsOffspring[N_ParentsOffspring];
  int<lower=0, upper=N_ped> NoOffspring[N_NoOffspring];
}

transformed data{
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);

  vector[N_ped] MSsd = sqrt(MSV);
  vector[N_ped] MSV_red; 
  MSV_red[NoParentsOffspring] = rep_vector(0.0,N_NoParentsOffspring);
  MSV_red[ParentsOffspring] = rep_vector(0.0,N_ParentsOffspring);
  MSV_red[NoOffspring] = MSV[NoOffspring];
}


parameters {
  vector [N_nest] nest_scaled;
  real<lower=0> sigma_nest;

  vector [N_ind] ind_scaled;
  real<lower=0> sigma_ind;

  vector[J] beta_tilde; // fixed effects
  vector [N_NoParentsOffspring] A_NoParentsOffspring;          // scaled breeding value
  vector [N_ParentsOffspring] A_ParentsOffspring;
  // vector [N_NoOffspring] A_NoOffspring;

  real <lower=0> sigma_A;          // sd of A
  real <lower=0> sigma_E;          // residual SD
}
 
// transformed parameters {
// }

model {
  vector [N] mu;
  vector [N_ped] A_scaled;          // scaled breeding value
  vector [N_ped] A;                    // breeding value
  vector [N_ped] sigma_ind_mix;
  int PO_i;
  int NO_j;
  real A_mu;
  real dam_A;
  real sire_A;
 
  A_NoParentsOffspring ~ normal( 0 , 1 );
  A_scaled[NoParentsOffspring] = A_NoParentsOffspring;

  for (i in 1:N_ParentsOffspring)
  {
    PO_i = ParentsOffspring[i];
 
    if (dam[PO_i] == 0 ) dam_A = 0;
    else  dam_A = A_scaled[dam[PO_i]];
    
    if (sire[PO_i] == 0 ) sire_A = 0;
    else  sire_A = A_scaled[sire[PO_i]];
    
    A_mu = (sire_A  + dam_A)*0.5;
    A_ParentsOffspring[i] ~ normal( A_mu, MSsd[PO_i]);
    A_scaled[PO_i] = A_ParentsOffspring[i];
  }

  for (j in 1:N_NoOffspring)
  {
    NO_j = NoOffspring[j];

    if (dam[NO_j] == 0 ) dam_A = 0;
    else  dam_A = A_scaled[dam[NO_j]];
    
    if (sire[NO_j] == 0 ) sire_A = 0;
    else  sire_A = A_scaled[sire[NO_j]];

    A_scaled[NO_j] = (sire_A  + dam_A)*0.5;
  }
  
  nest_scaled ~ normal( 0 , 1 );
  ind_scaled ~ normal( 0 , 1 );


  sigma_A ~ cauchy(0,10);
  sigma_E ~ cauchy(0,10);

  A = A_scaled * sigma_A;
  sigma_ind_mix = sqrt(MSV_red * sigma_A^2 + sigma_ind^2);

  mu = Q * beta_tilde + A[animal_id] + nest_scaled[nest_id] * sigma_nest  + ind_scaled[ind_id] .* sigma_ind_mix[animal_id];


  y ~ normal(mu, sigma_E) ;
  
  // y ~ normal(Q * beta_tilde + A[animal], sigma_E_mix[animal]) ;
  
}

generated quantities {
      vector[J] beta = R_inv * beta_tilde; // coefficients on x
}
