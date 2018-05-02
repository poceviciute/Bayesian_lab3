#Question 2

library(rstan)

ar_process = '
data {
  int<lower=0> iter;
  real<lower=0> sigma2;
  real<lower=-1, upper=1> phi;
  real mu;
}
transformed data {
  real x[iter];
  x[1] = mu;
}
parameters {
  real epsilon;
}
model {
  for(i in 2:iter){
    x[i] ~ normal(mu + phi*(x[i-1]-mu), sqrt(sigma2));
  }
}'

burnin = 100
niter = 200
fit1<-stan(model_code=ar_process,
           data=list(iter=200, sigma2=2, phi=0.5, mu=10),
           warmup=burnin,
           iter=niter,
           chains=1)