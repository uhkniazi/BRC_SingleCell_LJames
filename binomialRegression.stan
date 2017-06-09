data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    int y[Ntotal]; // response variable binomial distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
    real<lower=0.1> betaSigma; // standard deviation parameter for the joint prior for betas/coefficients
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X * betas; 
  mu = inv_logit(mu);
}
model {
  betaSigma ~ gamma(0.5, 0.0001);
  betas[1] ~ normal(0, 10); //prior for the betas
  betas[2:Ncol] ~ normal(0, betaSigma);
  // likelihood function
  y ~ bernoulli(mu);
}