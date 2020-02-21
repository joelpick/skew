  data {
  int N;
  vector[N] y;
  int K;
  matrix[N, K] X;
}

parameters {
  vector[K] beta;
  real<lower=4> nu_E;
  real<lower=0>  sigma_E;
}

transformed parameters {
  vector [N] mu = X*beta;
  real omega_E = sigma_E * sqrt((nu_E - 2) / nu_E);
}

model{
  beta ~ normal(0,100);
  sigma_E ~ cauchy(0, 10);
  nu_E ~ uniform(4,40);
  y ~ student_t(nu_E, mu, omega_E);
}
