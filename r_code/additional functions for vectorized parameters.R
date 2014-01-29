com.log.vecsum<-function(x,y){
#This function is modified from the "com.log.sum" function to make it able
#to deal with vectorized parameters
x<-as.vector(x)
y<-as.vector(y)
if (length(x)>1 || length(y)>1){
if (length(x) != length(y))
  stop("the length of parameters should be the same")
n<-length(x)
keepy<-rep(0,n)
for (i in 1:n){
keepy[i]<-com.log.vecsum(x[i], y[i])
  }
return(as.matrix(keepy))
 }
else {
x<-as.numeric(x)
y<-as.numeric(y)
if (x == -Inf){
return(y)
  }
else if (y == -Inf){
return(x)
  }
else if (x > y){
return(x + log(1 + exp(y - x)))
  }
else {
return(y + log(1 + exp(x - y)))
  }
 }
}


com.compute.log.vecz<-function(lambda, nu, log.error=0.001){
#This function is modified from the function"com.compute.log.z" to deal
#with vectorized parameters.
lambda<-as.vector(lambda)
nu<-as.vector(nu)
if (length(lambda)>1 || length(nu)>1){
if (length(lambda) != length(nu))
  stop("the length of parameters should be the same")
n<-length(lambda)
keepy<-rep(0,n)
for(i in 1:n){
keepy[i]<-com.compute.log.vecz(lambda[i], nu[i])
  }
return(as.matrix(keepy))
 }
else {
lambda<-as.numeric(lambda)
nu<-as.numeric(nu)
if (lambda <= 0 || nu < 0) 
  stop("Invalid arguments, only defined for lambda > 0, nu >= 0")
z = -Inf
z.last = 0
j = 0
while (abs(z - z.last) > log.error) {
z.last = z
z = com.log.sum(z, j * log(lambda) - nu * com.log.factorial(j))
j = j + 1
  }
return(z)
 }
}


com.compute.vecz<-function(lambda, nu){
#This code is modiefied from the function "com.compute.z" to deal 
#with vectorized parameters.
return(exp(com.compute.log.vecz(lambda, nu)))
}


