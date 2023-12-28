
X = scan('smooth.txt', skip=2, nlines=20)
Y = scan('smooth.txt', skip=26, nlines=20)

data = data.frame(x=X,y=Y)
# order
data = data[order(data$x),]
plot(y~x, data)

# test
for(i in 1:5){print(data[i,])} # loops in sorted order, good
data[1:5,]
k = 5
n = length(y)
for(i in 1:5){print(max(i-(k-1)/2, 1))}
for(i in 95:100){print(min(i+(k-1)/2, 100))}

# main
csrm = function(y, k){
    n = length(y)
    s.hat = rep(0, length(y))
    s = c()
    band.min = max(i-(k-1)/2, 1)
    band.max = min(i+(k-1)/2, n)
    for(i in 1:n){
        band.min = max(i-(k-1)/2, 1)
        band.max = min(i+(k-1)/2, n)
        band = y[band.min:band.max]
        s.hat[i] = mean(band)
        # s matrix
        s.row = rep(0, length(y))
        s.row[band.min:band.max] = 1/length(band)
        s = rbind(s, s.row)
    }
    return(list(s=s, s.hat=s.hat))
}

cvrss = function(kmax, y){
    # k vector
    sp.k = seq(3, kmax, 2)
    n = length(y)
    cv.k = numeric(length(sp.k))
    for(i in 1:length(sp.k)){
        tmp = csrm(y, sp.k[i]) # tmp$s = s matrix, tmp$s.hat = s hat
        sii = diag(tmp$s)
        s.hat = tmp$s.hat
        cv.k[i] = mean(((y-s.hat)/(1-sii))^2)
    }
    return(list(sp.k=sp.k, cvrss=cv.k)) # sp.k = samples of k, cvrss = cvrss result
}

par(mfrow=c(2,1))
# CVRSS
cvres = cvrss(23, data$y) # 2k+1, k<=11 -> k=23
cvres$sp.k
cvres$cvrss
plot(cvres$cvrss~cvres$sp.k, type='b', pch=15)
abline(v=5, col='red') # k=5

# constant span running mean
csrm.res = csrm(data$y, 5)
plot(y~x, data)
lines(data$x, csrm.res$s.hat, col='red')


### 2 median
csrmed = function(y, k){
    n = length(y)
    s.hat = rep(0, length(y))
    s = c()
    band.min = max(i-(k-1)/2, 1)
    band.max = min(i+(k-1)/2, n)
    for(i in 1:n){
        band.min = max(i-(k-1)/2, 1)
        band.max = min(i+(k-1)/2, n)
        band = y[band.min:band.max]
        s.hat[i] = median(band)
        # s matrix
        s.row = rep(0, length(y))
        s.row[band.min:band.max] = 1/length(band)
        s = rbind(s, s.row)
    }
    return(list(s=s, s.hat=s.hat))
}

cvrss = function(kmax, y){
    # k vector
    sp.k = seq(3, kmax, 2)
    n = length(y)
    cv.k = numeric(length(sp.k))
    for(i in 1:length(sp.k)){
        tmp = csrmed(y, sp.k[i]) # tmp$s = s matrix, tmp$s.hat = s hat
        sii = diag(tmp$s)
        s.hat = tmp$s.hat
        cv.k[i] = mean(((y-s.hat)/(1-sii))^2)
    }
    return(list(sp.k=sp.k, cvrss=cv.k)) # sp.k = samples of k, cvrss = cvrss result
}

# CVRSS
par(mfrow=c(1,1))
cvres = cvrss(23, data$y) # 2k+1, k<=11 -> k=23
cvres$sp.k
cvres$cvrss
plot(cvres$cvrss~cvres$sp.k, type='b', pch=15)
abline(v=13, col='red') # k=13

# constant span running median
csrm.res = csrmed(data$y, 13)
plot(y~x, data)
lines(data$x, csrm.res$s.hat, col='red')

install.packages("ramify")
library(ramify)
print_mat = matrix(csrm.res$s, nrow=length(data$y))
pprint(print_mat)

### 3
# kernel normal
ker = function(z){exp(-z^2/2)/(2*pi)^0.5}
h.seq = seq(0.1,1.1,0.1)

# kernel smoother

# ks = function(x, y, h){
#     n = length(y)
#     s.mat = c()
#     s.hat = numeric(n)
#     for(i in 1:n){
#         w = ker((x-x[i])/h)/sum(ker((x-x[i])/h))
#         w = dnorm((x-x[i])/h)/sum(dnorm((x-x[i])/h))
#         s.mat = rbind(s.mat, w)
#         s.hat[i] = sum(y[i]*w)
#     }
#     return(list(mat=s.mat, hat=s.hat))
# }

ks = function(x, y, h){
    n = length(y)
    s.mat = c()
    s.hat = numeric(n)
    denominator = 0
    for(i in 1:n){
      denominator= denominator+ dnorm((x-x[i])/h)
    }
    for(i in 1:n){
        w = dnorm((x-x[i])/h)/denominator
        s.mat = rbind(s.mat, w)
        s.hat[i] = sum(y[i]*w)
    }
    return(list(mat=s.mat, hat=s.hat))
}

par(mfrow=c(6,2))
for(bw in h.seq){
    kernel.res = ks(data$x, data$y, bw)
    plot(y~x, data=data, main=paste('bw=',bw))
    lines(data$x, kernel.res$hat, col='red')
}

h.seq

h.seq[3] # 0.3
pprint(matrix(kernel.res$mat, nrow=length(data$y)))

matrix(kernel.res$mat, nrow=length(data$y))[,3]
