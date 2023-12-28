
############# orthogonal approximation, chevychef
### 4.a

chevf = function(x){1/(1-x^2)^0.5}

integrate(chevf, lower=-1, upper=1)

integrate(function(x){x^3/(1-x^2)^-0.5}, -1, 1)


### 4.b
getwd()
data = scan('Orthogonal.txt')
data

chev.sum = function(x){
  sum(c(1,x,x^2-0.5,x^2-3/4*x))
}

qk = numeric(length(data))

for(i in 1:length(data)){
  qk[i] = chev.sum(data[i])
}

ck = 1/length(data)*sum(qk)

x=seq.int(-1,1,0.001)
f=numeric(length(x))
for(i in 1:length(x)){
  f[i]=ck*chev.sum(x)
}

plot(density(f))

A = matrix(c(4,1,0,1,4,1,0,1,4), nrow=3)
colnames(A) <- paste0('z', 1:3)
b=matrix(c(2,5/6,0))
solve(A,b)


###################### kernel estimation

### mod 11

# 10.1. Sanders et al. provide a comprehensive dataset of infrared emissions 
# and other characteristics of objects beyond our galaxy [567]. These data are 
# available from the website for our book. Let X denote the log of the variable 
# labeled F12, which is the total 12-µm-band flux measurement on each object.

# a. Fit a normal kernel density estimate for X, using bandwidths derived from 
# XXX the UCV(h) criterion XXX, Silverman’s rule of thumb, the Sheather–Jones approach, 
# Terrell’s maximal smoothing principle, and any other approaches you wish.
# Comment on the apparent suitability of each bandwidth for these data.

# b. Fit kernel density estimates for X using uniform, normal, Epanechnikov, 
# and triweight kernels, each with bandwidth equivalent to the Sheather–Jones 
# bandwidth for a normal kernel. Comment.

data = scan('F12.txt')
data

xi = log(data)

### a.

# silverman
silv.h = function(x){
  n = length(x)
  sd.hat = min(var(x)^0.5, IQR(x)/1.35)
  return(4/(3*n)^(1/5)*sd.hat)
}


# try implementing SJ
sj.h = function(x, xi, h0){
  n = length(xi)
  sj = numeric(n)

  for(i in 1:n){
    param = (x-xi[i])/h0
    pre = 1/(n*h0^3)
    L.dx2 = D(D(expression((1/(sqrt(2*pi)))*exp(-0.5*param^2)), 'param'), 'param')
    f.dx2 = pre * eval(L.dx2)
    sj = cbind(sj, f.dx2)
  }
  sj.Rdx2 = rowSums(sj)
  R.k = 1/(2*pi^0.5)
  return((R.k*(1/(n*var(sj.Rdx2)^2*sj.Rdx2)))^(1/5))
}
sj.h(x, xi, silv.h(xi))


# tarrell
tr.h = function(x){
  n = length(x)
  RK = 1/(2*pi^0.5)
  sd = var(x)^0.5
  return(3*(RK/(35*n))^(1/5)*sd)
}

tr.h(xi) # 0.2524386

### result
silv.h(xi) # 0.6919854

bw.SJ(x = x, method = "dpi") # 0.3685721

tr.h(xi) # 0.2524386


### graph
# silverman
h = silv.h(xi)
x = seq(min(xi)-(3*h), max(xi)+(3*h), length.out=length(xi))

hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (exp(-z^2/2)/((2*pi)^0.5))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue', alpha=0.5)
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')


# sj
h = bw.SJ(x = x, method = "dpi")
x = seq(min(xi)-(3*h), max(xi)+(3*h), length.out=length(xi))

hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (exp(-z^2/2)/((2*pi)^0.5))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')

# tarrell
h = tr.h(xi)
x = seq(min(xi)-(3*h), max(xi)+(3*h), length.out=length(xi))

hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (exp(-z^2/2)/((2*pi)^0.5))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')


### b.
# uniform, normal, Epanechnikov, and triweight kernels
K.unif = function(z){
  z[abs(z)>=1] = 0
  z[z!=0] = 1
  return(0.5*z)
}
K.norm = function(z){
  return(exp(-z^2/2)/(2*pi)^0.5)
}
K.Epkv = function(z){
  f=numeric(length(z))
  for (i in 1:length(z)){
    if (abs(z[i])>=1){
      f[i]=0
    }
    else{
      f[i]= 3/4*(1-z[i]^2)
    }
  }
  return(f)
}
K.tri = function(z){
  f=numeric(length(z))
  for (i in 1:length(z)){
    if (abs(z[i])>=1){
      f[i]=0
    }
    else{
      f[i]=1-abs(z[i])
    }
  }
  return(f)
}
h = silv.h(xi)
x = seq(min(xi)-(3*h), max(xi)+(3*h), length.out=length(xi))

# uniform
hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (K.unif(z))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')

# normal
hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (K.norm(z))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')

# Epanechnikov
hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (K.Epkv(z))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')

# triweight
hist(xi, nclass=30, prob=TRUE)
rug(xi)

ke = numeric(length(xi))
n = length(xi)

for(i in 1:length(xi)){
  z = (x-xi[i])/h
  fhat = (K.tri(z))
  ke = cbind(ke, fhat) # 
  lines(x, fhat, lwd=1, col='blue')
}
kde = 1/(n*h)*rowSums(ke)
lines(x, kde, lwd = 3, col = 'red')

# histogram
hist(xi, nclass=30, prob=TRUE)

n = length(xi)
pk = numeric((range(xi)[2]-range(xi)[1])/h+1)
last.bin = min(xi)
x.bin = c()

for(i in 1:length(bins)){
  bin.i = xi[(last.bin<=xi)&(xi<last.bin+h)]
  pk[i] = length(bin.i)/(n)

  last.bin = last.bin + h
  x.bin[i] = last.bin
}
lines(x.bin, pk, lwd=3, col='red')



