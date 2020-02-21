functions{
  real gamma_to_delta(real skew){
    real b = sqrt(2/pi());
    real a2 = cbrt(2*b^3-b)^2;
    return(cbrt(skew) / sqrt(a2 + b^2 * cbrt(skew)^2));
  }
} 

data {
  int<lower=0> N;  // number of observations
  vector[N] y;     // response
  int K;
  matrix [N, K] X ; 
}

transformed data{
  real b = sqrt(2/pi());
}

parameters {
  vector [K]beta;
  real <lower=0> sigma_E;  
  real <lower=-0.99, upper=0.99> gamma_E;
}

transformed parameters {
  real <lower=-1, upper=1> delta_E = gamma_to_delta(gamma_E);
  real alpha_E = delta_E/(sqrt(1-delta_E^2));
  real <lower=0>  omega_E = sigma_E / sqrt(1-(b*delta_E)^2);
  real xi_E = 0 - omega_E * delta_E * b;
  vector [N] mu = xi_E + X*beta;
  }

model {
  beta ~ normal(0, 100);
  sigma_E ~ cauchy(0,10);
  gamma_E ~ uniform(-0.99, 0.99);
  y ~ skew_normal(mu, omega_E, alpha_E);
}

