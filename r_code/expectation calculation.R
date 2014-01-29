#This function is used to calculate the expectation value of the COM-Poisson distribution 
#with the parameter values lambda and nu, respectively. Notice that when lambda < 1 and
#lambda > 10^nu, the normalizing constant z is possible to be infinity in R computation, thus 
#making the expectation computation incorrect. So when the pramameters fall to these bounds, 
#the expected value may have to be attained by approximation.  
expectation<-function(lambda,nu){
lambda<-as.vector(lambda)
nu<-as.vector(nu)
if (length(lambda) > 1 || length(nu) > 1){
if (length(lambda) != length(nu))
stop("parameter vectors should have the same length")
n<-length(lambda)
slot<-rep(0,n)
for (i in 1:n){
slot[i]<-expectation(lambda[i],nu[i])
  }
slot<-as.matrix(slot)
return(slot)
 }
else {
lambda<-as.numeric(lambda)
nu<-as.numeric(nu)
expect<-com.expect.compute.z(lambda,nu)/com.compute.vecz(lambda,nu)
return(expect)
 }
}
