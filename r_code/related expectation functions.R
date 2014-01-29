com.expect.compute.log.z<-function (lambda, nu, log.error=0.001){
#This code computes the logarithm of the numerator of the expression for the expectation of y(i), which is the value of the data in the ith group. 
#note that this function is meant to deal with vectorized parameters.
lambda<-as.vector(lambda)
nu<-as.vector(nu)
if (length(lambda) != length(nu))
  stop("parameters need to be the same length")
n<-length(lambda)
if (n > 1){
keepy<-rep(0,n)
for(i in 1:n){
keepy[i]<-com.expect.compute.log.z(lambda[i], nu[i])
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
j = 1
while (abs(z-z.last)>log.error) {
z.last = z
z = com.log.sum(z, log(j)+j * log(lambda) - nu * com.log.factorial(j))
j = j + 1
  }
return(z)
 }
}


com.expect.compute.z<-function (lambda, nu) 
#This code computes the value of the numerator of expression for the expectation of y(i), which is the value of the data in the ith group. 
{
    return(exp(com.expect.compute.log.z(lambda, nu)))
}


com.logexpect.compute.log.z<-function (lambda, nu, log.error=0.001){
#See the explanation of the code in "com.logexpect.compute.z", this code takes the logarithm of that expression.
#note that this function is meant to deal with vectorized parameters
lambda<-as.vector(lambda)
nu<-as.vector(nu)
if (length(lambda) != length(nu))
  stop("parameters need to be the same length")
n<-length(lambda)
if (n > 1){
keepy<-rep(0,n)
for(i in 1:n){
keepy[i]<-com.logexpect.compute.log.z(lambda[i], nu[i])
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
j = 2
while (abs(z-z.last)>log.error) {
z.last=z
z = com.log.sum(z, log(com.log.factorial(j))+j * log(lambda) - nu * com.log.factorial(j))
j = j + 1
  }
    return(z)
 }
}


com.logexpect.compute.z<-function (lambda, nu) 
#This code computes the expectation of the logarithm of the data in the ith group. See the expression in the derivation, 
#or more precisely, the term got by taking the derivatives of the likelihood function with respect to nu for more details.
{
    return(exp(com.logexpect.compute.log.z(lambda, nu)))
}
