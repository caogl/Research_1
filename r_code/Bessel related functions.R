com.compute.log.bessel <- function(z, alpha, nu, log.error=0.001){
#This code computes the log of the generalized version of the modified Bessel function of the 
#1st kind in light of the dispersion parameter, nu.  For the special case where nu=1, this is 
#precisely the log of the modified Bessel function of the 1st kind.
#note that this function is meant to deal with vectorized parameters.
alpha<-abs(alpha)
z<-as.vector(z)
alpha<-as.vector(alpha)
nu<-as.vector(nu)
if (length(z) != length(alpha) || length(z) != length(nu))
  stop ("parameters need to be the same length")
n<-length(z)
if (n > 1){
keepy<-rep(0,n)
for (i in 1:n){
keepy[i]<-com.compute.log.bessel(z[i], alpha[i], nu[i])
  }
return(as.matrix(keepy))
 }
else {
z<-as.numeric(z)
alpha<-as.numeric(alpha)
nu<-as.numeric(nu) 
if (z <= 0 || nu < 0) 
  stop("Invalid arguments, only defined z > 0 and nu >= 0 ")
I = -Inf
I.last = 0
j=0
while (abs(I - I.last) > log.error) {
I.last = I
I = com.log.sum(I, ((2*j)+alpha) * (log(z) - log(2)) - nu * (com.log.factorial(j+alpha) +  com.log.factorial(j)))
j = j + 1
   }
return(I)
 }
}


com.compute.bessel <- function(z, alpha, nu){
# This code computes the generalized version of the modified Bessel function of the 1st kind
# in light of the dispersion parameter, nu.  For the special case where nu=1, this is precisely
# the modified Bessel function of the 1st kind.
return(exp(com.compute.log.bessel(z, alpha, nu)))
}


com.compute.log.derivbessel<-function(z,alpha,nu,log.error=0.001){
#See the explanation of the code in "cms.compute.derivbessel", this code takes the logarithm of that expression.
#note that this function is meant to deal with vectorized parameters.
alpha<-abs(alpha)
z<-as.vector(z)
alpha<-as.vector(alpha)
nu<-as.vector(nu)
if (length(z) != length(alpha) || length(z) != length(nu))
         stop ("parameters need to be the same length")
n<-length(z)
if (n > 1){
keepy<-rep(0,n)
for (i in 1:n){
keepy[i]<-com.compute.log.derivbessel(z[i], alpha[i], nu[i])
  }
return(as.matrix(keepy))
 }
else {
z<-as.numeric(z)
alpha<-as.numeric(alpha)
nu<-as.numeric(nu)
if (z <= 0 || nu < 0) 
        stop("Invalid arguments, only defined for z > 0 and nu >= 0")
I=-Inf
I.last=0
if (alpha == 0){
j=1
    }
else {
j=0
    }
while(abs(I-I.last)>log.error){
I.last=I
I=com.log.sum(I,log(2*j+alpha)-nu*com.log.factorial(j+alpha)-nu*com.log.factorial(j)+(2*j+alpha-1)*log(z))
j=j+1
  }
return(I)
 }
}


com.compute.derivbessel<-function(z, alpha, nu){
#This code computes the numerator of the infinite series terms in the expression got by taking the partial derivatives of 
#the generalized form of the modified Bessel function of the first kind with respect to lambda1 or lambda2. See the expression 
#in the derivation, or more precisely, the term got by taking the derivatives of the likelihood function with respect to lambda1
#and lambda2 for more details. 
return(exp(com.compute.log.derivbessel(z, alpha, nu)))
}


com.compute.log.derivvbessel<-function(z, alpha, nu, log.error=0.001){
#See the explanation of the code in "cms.compute.derivvbessel", this code takes the logarithm of that expression.
#Note that this function is meant to deal with vectorized parameters rather than single values.
alpha<-abs(alpha)
z<-as.vector(z)
alpha<-as.vector(alpha)
nu<-as.vector(nu)
if (length(z) != length(alpha) || length(z) != length(nu))
         stop ("parameters need to be the same length")
n<-length(z)
if (n > 1){
keepy<-rep(0,n)
for (i in 1:n){
keepy[i]<-com.compute.log.derivvbessel(z[i], alpha[i], nu[i])
  }
return(as.matrix(keepy))
 }
else {
z<-as.numeric(z)
alpha<-as.numeric(alpha)
nu<-as.numeric(nu)
if (z <= 0 || nu < 0) 
        stop("Invalid arguments, only defined for z > 0 and nu >= 0")
I=-Inf
I.last=0
#note here the following code is for the correct calculation of the first term/terms of infinite series
if (alpha == 0){
j=2
   }
else if (alpha == 1){
j=1
   }
else {
j=0
   }
while(abs(I-I.last)>log.error){
I.last=I
I=com.log.sum(I,log(com.log.factorial(j+alpha)+com.log.factorial(j))-nu*com.log.factorial(j+alpha)
-nu*com.log.factorial(j)+(2*j+alpha)*log(z))
j=j+1
  }
return(I)
 }
}


com.compute.derivvbessel<-function(z, alpha, nu){
#This code computes the numerator term of the expression got by taking the partial derivatives of the generalized form 
#of the modified Bessel function of the first kind with respect to nu. See the expression in the derivation, or more precisely, 
#the term got by taking the derivatives of the likelihood function with respect to nu for more details.
return(exp(com.compute.log.derivvbessel(z, alpha, nu)))
}


