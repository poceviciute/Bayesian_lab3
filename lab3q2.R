#Question 2

library(rstan)

ar_process = '
data {
  int<lower=0> iter;
  real mu;
  real<lower=0> sigma2;
  real<lower=-1, upper=1> phi;
}
parameters {
  real x[iter];
}
model {
  for(i in 2:iter){
    x[i] ~ normal(mu + phi*(x[i-1]-mu), sqrt(sigma2));
  }
}'

burnin <- 100
niter <- 200
mu <- 10
init_para <- list(list(x=c(mu,rep(0,niter-1))))
#needs to be a list of a list since each chain needs its own list
fit1<-stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=0.5, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1)
print(fit1)
fit1_postmean <- get_posterior_mean(fit1)[burnin:niter]
plot(fit1_postmean, type="l", xaxt="n", xlab="Draws", ylab="Posterior mean", main="Posterior mean of AR(1) process")
axis(1, at=seq(0,(n-burn_in), by=10), labels=seq(burn_in,n, by=10))
