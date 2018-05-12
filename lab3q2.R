#Question 2

library(rstan)

#a)

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

burnin <- 101
niter <- 300
mu <- 10
init_para <- list(list(x=c(mu,rep(0,niter-1))))
#needs to be a list of a list since each chain needs its own list
fit1<-stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=0.5, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1, control = list(max_treedepth = 15))
print(fit1)
fit1_postmean <- get_posterior_mean(fit1)[burnin:niter]
plot(fit1_postmean, type="l", xaxt="n", xlab="Draws", ylab="Posterior mean", main="Posterior mean, phi=0.5")
axis(1, at=seq(0,(niter-burnin), by=10), labels=seq(burnin,niter, by=10))

fit2<-stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=0.9, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1, control = list(max_treedepth = 15))
print(fit2)
fit2_postmean <- get_posterior_mean(fit2)[burnin:niter]
plot(fit2_postmean, type="l", xaxt="n", xlab="Draws", ylab="Posterior mean", main="Posterior mean, phi=0.9")
axis(1, at=seq(0,(niter-burnin), by=10), labels=seq(burnin,niter, by=10))

fit3<-stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=-0.9, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1)

print(fit3)
fit3_postmean <- get_posterior_mean(fit3)[burnin:niter]
plot(fit3_postmean, type="l", xaxt="n", xlab="Draws", ylab="Posterior mean", main="Posterior mean, phi=-0.9")
axis(1, at=seq(0,(niter-burnin), by=10), labels=seq(burnin,niter, by=10))

fit4<-stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=0, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1)
print(fit4)
fit4_postmean <- get_posterior_mean(fit4)[burnin:niter]
plot(fit4_postmean, type="l", xaxt="n", xlab="Draws", ylab="Posterior mean", main="Posterior mean, phi=0")
axis(1, at=seq(0,(niter-burnin), by=10), labels=seq(burnin,niter, by=10))


#b)

xT <- stan(model_code=ar_process,
           data=list(iter=niter,sigma2=2, phi=0.3, mu=10),
           init = init_para,
           warmup=burnin,
           iter=niter,
           chains = 1)
print(xT)
xT_postmean <- get_posterior_mean(xT)[burnin:niter]

yT <-  stan(model_code=ar_process,
            data=list(iter=niter,sigma2=2, phi=0.95, mu=10),
            init = init_para,
            warmup=burnin,
            iter=niter,
            chains = 1)
print(yT)
yT_postmean <- get_posterior_mean(yT)[burnin:niter]

ar_mcmc = '
data {
  int<lower=0> iter;
  vector[iter] x;
}
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=-1, upper=1> phi;
}
model {
  mu ~ uniform(0,20);
  sigma ~ uniform(0,5);
  phi ~ normal(0, 1);
  for(i in 2:iter){
    x[i] ~ normal(mu + phi*(x[i-1]-mu), sigma);
  }
}'


burnin <- 100
niter <- 300
ndraws <- 200
#init_prior <- list(list(mu=0, phi=0.5, sigma=2))

xT_fit <- stan(model_code=ar_mcmc,
           data=list(iter=ndraws,x=xT_postmean),
           warmup=burnin,
           iter=niter,
           chains = 1)
print(xT_fit)
para_postmean <- get_posterior_mean(xT_fit)[burnin:niter]
mu_post2 <- extract(xT_fit)$mu[burnin:niter]
phi_post2 <- extract(xT_fit)$phi[burnin:niter]
plot(mu_post2, phi_post2, pch=19, xlab="mu", ylab="phi", main="Joint posterior of mu and phi for x")
perc1x <- sapply(as.data.frame(xT_fit), FUN=quantile, probs=0.025)[1:3]
perc2x <- sapply(as.data.frame(xT_fit), FUN=quantile, probs=0.975)[1:3]
n_eff_x <- summary(xT_fit)$summary[,"n_eff"][1:3]

yT_fit <- stan(model_code=ar_mcmc,
               data=list(iter=ndraws,x=yT_postmean),
               warmup=burnin,
               iter=niter,
               chains = 1)
print(yT_fit)
Ypara_postmean <- get_posterior_mean(yT_fit)[burnin:niter]

mu_post <- extract(yT_fit)$mu[burnin:niter]
phi_post <- extract(yT_fit)$phi[burnin:niter]
plot(mu_post, phi_post, pch=19, xlab="mu", ylab="phi", main="Joint posterior of mu and phi for y")
perc1y <- sapply(as.data.frame(yT_fit), FUN=quantile, probs=0.025)[1:3]
perc2y <- sapply(as.data.frame(yT_fit), FUN=quantile, probs=0.975)[1:3]
n_eff_y<- summary(yT_fit)$summary[,"n_eff"][1:3]
#c)

campy<-read.table("D:/LiU/732A91/Bayesian_lab3/Campy.dat",header=TRUE)


campy_model = '
data {
  int<lower=0> iter;
  int c[iter];
  real mu;
  real<lower=0> sigma;
  real<lower=1> lambda;
  real<lower=-1, upper=1> phi;
}
parameters {
  real x[iter];
}

model {
  for(i in 2:iter){
    x[i] ~ normal(mu + phi*(x[i-1]-mu), sigma/lambda);
    c[i] ~ poisson(exp(x[i]));
  }
    
}'

c_fit <- stan(model_code=campy_model,
                     data=list(iter=140,c=campy$c, phi=0.5, sigma=sqrt(2), lambda=1, mu=10),
                     warmup=30,
                     iter=140,
                     chains = 1)

print(c_fit)
theta <- exp(get_posterior_mean(c_fit))
perc1 <- exp(sapply(as.data.frame(c_fit), FUN=quantile, probs=0.025))
perc2 <- exp(sapply(as.data.frame(c_fit), FUN=quantile, probs=0.975))
plot(theta, type="l", main="Campy")
lines(campy$c, col="red")
lines(perc1, col="blue")
lines(perc2, col="blue")
extract(c_fit)


#d)

d_fit <- stan(model_code=campy_model,
              data=list(iter=140,c=campy$c, phi=0.5, sigma=sqrt(2), lambda=3, mu=10),
              warmup=30,
              iter=140,
              chains = 1)

print(d_fit)
theta2 <- exp(get_posterior_mean(d_fit))
perc1d <- exp(sapply(as.data.frame(d_fit), FUN=quantile, probs=0.025))
perc2d <- exp(sapply(as.data.frame(d_fit), FUN=quantile, probs=0.975))
plot(theta2, type="l", main="Campy")
lines(campy$c, col="red")
lines(perc1d, col="blue")
lines(perc2d, col="blue")

plot(theta, type="l")
lines(theta2, col="red")
