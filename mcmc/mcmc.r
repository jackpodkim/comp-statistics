
# MH
n = 10000

mix_f = function(x){return(0.7*dnorm(x,7,0.52)+(1-0.7)*dnorm(x,10,0.52))}

mh = function(x){
  x_t1 = rnorm(1, x, 0.01) 
  r = min(1, mix_f(x_t1)/mix_f(x))
  if(r<1){
    ifelse(runif(1)<r, return(x_t1), return(x))
  }
  return(x_t1)
}

#######
x=rep(0, n)
x[1]=0 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1])
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')
par(mfrow=c(2,1))
plot(x)
hist(x, main='start from 0')

#######
x=rep(0, n)
x[1]=7 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1])
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')
par(mfrow=c(2,1))
plot(x)
hist(x, main='start from 7')

#######
x=rep(0, n)
x[1]=15 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1])
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')
par(mfrow=c(2,1))
plot(x)
hist(x, main='start from 15')


mh = function(x){
  x_t1 = runif(1, 0, 20) 
  r = min(1, mix_f(x_t1)/mix_f(x))
  if(r<1){
    ifelse(runif(1)<r, return(x_t1), return(x))
  }
  return(x_t1)
}
x=rep(0, n)
x[1]=7 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1])
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')
par(mfrow=c(2,1))
plot(x)
hist(x)


############ Mod6
# MH


# The Normal distribution mean = x(t), sd = .1, so N(x(t), .1^2) .
# The Normal distribution mean = x(t), sd = .01, so N(x(t), .01^2) .
# The Normal distribution mean = x(t), sd = .001, so N(x(t), .001^2) .
# The Uniform on (0,5)

# define func
mix_f = function(x){return(dexp(x, 1))}
mh = function(x, t_sd){
  x_t1 = rnorm(1, x, t_sd) 
  r = min(1, mix_f(x_t1)/mix_f(x))
  if(r<1){
    ifelse(runif(1)<r, return(x_t1), return(x))
  }
  return(x_t1)
}
cusum_ = function(x){
  n = length(x)
  theta_hat = mean(x)
  pre_cusum = rep(0, n)
  for(i in 1:n){
    pre_cusum[i] = (x[i]-theta_hat)
  }
  return(cumsum(pre_cusum))
}

# The Normal distribution mean = x(t), sd = .1, so N(x(t), .1^2) .
n = 10000
t_sd = 0.1

# loop
x=rep(0, n)
x[1]=0 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1], t_sd)
  if(x[i]!=x[i-1]){ar=ar+1}
}

par(mfrow=c(3,1))
plot(x, type='l',)
plot(cusum_(x))
acf(x, lag.max=n)
mtext("N(x(t), .1^2)", side = 3, line = - 2, outer = TRUE)

# The Normal distribution mean = x(t), sd = .01, so N(x(t), .01^2) .
n = 10000
t_sd = 0.01

# loop
x=rep(0, n)
x[1]=0 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1], t_sd)
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')

par(mfrow=c(3,1))
plot(x, type='l',)
plot(cusum_(x))
acf(x, lag.max=n)
mtext("N(x(t), .01^2)", side = 3, line = - 2, outer = TRUE)

# The Normal distribution mean = x(t), sd = .001, so N(x(t), .001^2) .
n = 10000
t_sd = 0.001

# loop
x=rep(0, n)
x[1]=0 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1], t_sd)
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')

par(mfrow=c(3,1))
plot(x, type='l',)
plot(cusum_(x))
acf(x, lag.max=n)
mtext("N(x(t), .001^2)", side = 3, line = - 2, outer = TRUE)

# The Uniform on (0,5)
mh = function(x, t_sd){
  x_t1 = runif(1, 0, 5) 
  r = min(1, mix_f(x_t1)/mix_f(x))
  if(r<1){
    ifelse(runif(1)<r, return(x_t1), return(x))
  }
  return(x_t1)
}

n = 10000
t_sd = 0.001

# loop
x=rep(0, n)
x[1]=0 # start value
ar = 0
for(i in 2:n){
  x[i]=mh(x[i-1], t_sd)
  if(x[i]!=x[i-1]){ar=ar+1}
}
cat('Acceptance ratio:', ar/n, '\n')

par(mfrow=c(3,1))
plot(x, type='l',)
plot(cusum_(x))
acf(x, lag.max=n)
mtext("Uniform on (0,5)", side = 3, line = - 2, outer = TRUE)


########### Gibbs
# set vars
n=10000
x1 = x2 = numeric(n)
mu1 = mu2 = 0
std1 = std2 = 1 
burnin = 500

corr = 0 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p0 = x1[-(1:burnin)]
x2.p0 = x2[-(1:burnin)]

corr = 0.1 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p1 = x1[-(1:burnin)]
x2.p1 = x2[-(1:burnin)]

corr = 0.2 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p2 = x1[-(1:burnin)]
x2.p2 = x2[-(1:burnin)]

corr = 0.3 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p3 = x1[-(1:burnin)]
x2.p3 = x2[-(1:burnin)]

corr = 0.4 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p4 = x1[-(1:burnin)]
x2.p4 = x2[-(1:burnin)]

corr = 0.5 # ρ = 0, .1, .2, .3, .4, .5
for(i in 2:n){
  x1[i] = rnorm(1, mu1+(x2[i-1]-mu2)*corr, std1^2*(1-corr^2))
  x2[i] = rnorm(1, mu2+(x1[i]-mu1)*corr, std2^2*(1-corr^2))
}
x1.p5 = x1[-(1:burnin)]
x2.p5 = x2[-(1:burnin)]

par(mfrow=c(6,2))
plot(x1.p0, type='l'); title('p=0'); hist(x1.p0)
plot(x2.p0, type='l'); title('p=0'); hist(x2.p0)
plot(x1.p1, type='l'); title('p=0.1'); hist(x1.p1)
plot(x2.p1, type='l'); title('p=0.1'); hist(x2.p1)
plot(x1.p2, type='l'); title('p=0.2'); hist(x1.p2)
plot(x2.p2, type='l'); title('p=0.2'); hist(x2.p2)
par(mfrow=c(6,2))
plot(x1.p3, type='l'); title('p=0.3'); hist(x1.p3)
plot(x2.p3, type='l'); title('p=0.3'); hist(x2.p3)
plot(x1.p4, type='l'); title('p=0.4'); hist(x1.p4)
plot(x2.p4, type='l'); title('p=0.4'); hist(x2.p4)
plot(x1.p5, type='l'); title('p=0.5'); hist(x1.p5)
plot(x2.p5, type='l'); title('p=0.5'); hist(x2.p5)


### 6.3
