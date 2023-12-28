
### jackknife

### 1.a
plot(x, y)

f = function(par, x){
    m = par[1]
    s = par[2]
    n = length(x)
    l =  -(n/2)*(log(2*pi*s^2))+(-1/(2*s^2)) * sum((x-m)^2)
    return(-l)
}

mle = optim(par=c(0, 2^0.5), fn=f, x=x)
mle$par
mean(x); sd(x)
# > mle$par
# [1] 0.1735848 1.4267681

g = dnorm(x=x, mean=mle$par[1], sd=mle$par[2])
R = abs(y-g)
E_pxy.hat = sum(R)/length(y)
E_pxy.hat

# > E_pxy.hat
# [1] 0.02524844

### 1.b
s1 = x[1:50]
s2 = x[51:100]
y1 = y[1:50]
y2 = y[51:100]

# 1/s1 sum_from_y2(R(y2-g1)), vice versa
g = dnorm(x=s1, mean=mle$par[1], sd=mle$par[2])
R = abs(y2-g)
E_pxy.hat1 = sum(R)/length(y1)
g = dnorm(x=s2, mean=mle$par[1], sd=mle$par[2])
R = abs(y1-g)
E_pxy.hat2 = sum(R)/length(y2)

E_pxy.hat.cv = (E_pxy.hat1+E_pxy.hat2)/2
E_pxy.hat.cv
# > E_pxy.hat.cv
# [1] 0.08432905

### 1.c
# CV's error estimate is actually bigger (0.025 vs 0.084) than without-split, 
# but might be a better estimator of true error if the splits are more granular.


### 3.b
# (b) Posted on the course Canvas site is the data set Jackknife.txt containing 100 obser-
# vations from a N (0, 1) distribution. Use this sample to calculate b2 and the jackknife
# estimate of the standard deviation for the cases k = 1 and k = 5.
y = scan('7.Computational_stat/mod8/Jackknife.txt')
y

b2 = function(data){
    n = length(data)
    m = mean(data)
    return(sum(data-m)^4 / (sum(data-m)^2)^2)
}

jk.std = function(data, k){
    n = length(data)
    j.res = numeric(n-k)
    for(i in 1:(n-k)){
        t_j_sp = data[-(i:(i+k-1))]
        t_j_mean = mean(t_j_sp)
        j.res[i] = ((1/(n-k))*(sum((t_j_sp - t_j_mean)^2)))^0.5
    }
    return(mean(j.res))
}
jk.std(y, 2)

for(k in 1:5){
    print(c(b2(y), jk.std(y, k)))
}


### bootstrap

stomach = c(25, 42, 45, 46, 51, 103, 124, 146, 340, 396, 412, 876, 1112)
breast = c(24, 40, 719, 727, 791, 1166, 1235, 1581, 1804, 3460, 3808)

stomach = log(stomach)
breast = log(breast)

calculate.bs.ci = function(data, m, alpha){
  bs.t = numeric(m)
  bs.se = numeric(m)
  for(i in 1:m){
    bs = sample(data, length(data), replace=TRUE)
    bs.se[i] = sum((bs-mean(bs))^2)^0.5/length(data)
    bs.t[i] = (mean(bs.mu)-mean(data))/bs.se[i]
  }
  bs.qt = quantile(bs.t, c(alpha/2, 1-alpha/2))
  bs.sd = sd(bs.se)
  return(c(mean(stomach)-bs.qt[2]*bs.sd, mean(stomach)-bs.qt[1]*bs.sd))
}

stomach.bs.ci = calculate.bs.ci(stomach, 10000, 0.05)
stomach.bs.ci

breast.bs.ci = calculate.bs.ci(breast, 10000, 0.05)
breast.bs.ci

### 1.2
# single elements exchanged
base = abs(mean(stomach)-mean(breast))

cnt = 0
for(i in 1:length(stomach)){
  for(j in 1:length(breast)){
    sim.st = stomach
    sim.br = breast
    # exchange
    sim.st[i] = sim.br[j] 
    sim.br[j] = sim.st[i]
    if(abs(mean(sim.st)-mean(sim.br))>=base){cnt=cnt+1}
  }
}

p.val = cnt/(length(stomach)*length(breast))
p.val


### 1.3
# exp
stomach.bs.ci.exp = exp(stomach.bs.ci)
breast.bs.ci.exp = exp(breast.bs.ci)

stomach = c(25, 42, 45, 46, 51, 103, 124, 146, 340, 396, 412, 876, 1112)
breast = c(24, 40, 719, 727, 791, 1166, 1235, 1581, 1804, 3460, 3808)
stomach.bs.ci.orig = calculate.bs.ci(stomach, 10000, 0.05)
breast.bs.ci.orig = calculate.bs.ci(breast, 10000, 0.05)


stomach.bs.ci
breast.bs.ci
stomach.bs.ci.exp
breast.bs.ci.exp
stomach.bs.ci.orig
breast.bs.ci.orig

### 2.

bs <- function(data, n){
  data = as.matrix(data)
  data.row = nrow(data)
  data.bs = replicate(
    n, 
    data[sample.int(data.row,replace=TRUE)]
    )
  bs.sd = apply(data.bs, 2, sd)
  return(list(sp=data.bs, se=bs.sd))
}

u = runif(1000,0,1)
bs.u = bs(u, 1000)

bs.u$se
hist(bs.u$se)

sd(bs.u$se)

a = c(1,2,3,4,5,6,7,8,9)
b = bs(a, 7)
b

b$sp
b$se

b[,3]

colMeans(b)

apply(b, 2, sd)


### 4.
## a.

sp = rnorm(100)
mu.hat = mean(sp)
mu.hat

# > mu.hat
# [1] -0.0947098

## b.

bs <- function(data, n){
  data = as.matrix(data)
  data.row = nrow(data)
  data.bs = replicate(
    n, 
    data[sample.int(data.row,replace=TRUE)]
    )
  bs.mu = apply(data.bs, 2, mean)
  return(list(sp=data.bs, mu=bs.mu))
}

b=10
bs.sp = bs(sp, b)
> bs.sp$mu
 [1] -0.039310598 -0.055443623 -0.126996449 -0.172284374 -0.010338145
 [6]  0.025548969 -0.131902034 -0.019532689 -0.114702431 -0.005179032

mu.bias = mean(bs.sp$mu)-mu.hat
# > mu.bias
# [1] 0.02969576

mu.var = sum((bs.sp$mu-mean(bs.sp$mu))^2)/(b-1)
# > mu.var
# [1] 0.004431537

## c.
unbiased.bs <- function(data, n, bias){
  data = as.matrix(data)
  data.row = nrow(data)
  data.bs = replicate(
    n, 
    data[sample.int(data.row,replace=TRUE)]
    )
  bs.mu = apply(data.bs, 2, mean)
  bs.mu = bs.mu - bias
  return(list(sp=data.bs, mu=bs.mu))
}
b=10
bs.sp = unbiased.bs(sp, b, mu.bias)
# > bs.sp$mu
#  [1] -0.17939485 -0.17936911 -0.19590414 -0.27883607  0.05863648 -0.04170347
#  [7] -0.08511879 -0.03678423 -0.12182360 -0.05233852

mu.bias = mean(bs.sp$mu)-mu.hat
# > mu.bias
# [1] -0.01655383

mu.var = sum((bs.sp$mu-mean(bs.sp$mu))^2)/(b-1)
# > mu.var
# [1] 0.009782648
