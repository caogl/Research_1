#This code is used to generate n random numbers from the COM-Poisson distribution when the parameter 
#values are lambda and nu, respectively
#notice that this is different from most inherit random number generating functions.
#when the parameters lambda and nu are numbers, this function is meant to generate n random numbers follow
#the Conway-Maxwell-Poisson distribution with the parameter values alpha and nu, respectively.
#However, note that this random number generating function is specially designed to fit the "COM-Skellam 
#regression", when the parameters lambda and nu are vectors of the same length, which is assumed to be n,
#this functions is then meant to generate n random numbers, the first of them is generated from the 
#COM-Poisson distribution with parameter values, which are first number of the parameter vectors lambda 
#and nu, so is the second, third, etc.
#Notice that in this function, when the condtions that parameter value lambda > 10^nu and nu < 1 are both 
#satisfied, only the approximated expectation, rather than the random number is generated, because when the 
#paramaters fall into these value bounds, the approximation would be quiet good and more importantly, z tends 
#to be infinity or very large, making the random number generating algorithm lose its effect. 
makeCMPdata.gen <- function(n,lambda,nu){
lambda<-as.vector(lambda)
nu<-as.vector(nu)
if (length(lambda)>1 || length(nu)>1){
if (length(lambda) != length(nu))
stop ("the length of parameters must be the same")
if (n != length(lambda))
stop ("the numbers of generated random numbers should be equal to the length of the parameters to produce
one on one correspondence random number, for the EIM method")
rand <- rep(0,n)
for (j in 1:n){
rand[j]<-makeCMPdata.gen(1,lambda[j],nu[j])
  }
return(rand)
 }
else {
lambda<-as.numeric(lambda)
nu<-as.numeric(nu)
if (lambda<=0 || nu<0) 
stop ("Invalid argument, only defined for lambda > 0 and nu >= 0")
#generate n uniform distribution values
unifvals <- runif(n)
#Create space for simulated y's each time with the corresponding probability value
keepy <- rep(0,n)
for (i in 1:n){
# start counter for y
y <- 0
#Compute Z-inverse.  This equals P(Y=0).
zinv <- 1/com.compute.vecz(lambda, nu)
py <- zinv
while (py < unifvals[i]){
y <- y+1
py <- py + ((lambda^y)/(((factorial(y))^nu))*zinv)
    }
#Keep last value of y, where py >= unifvals[i]
keepy[i] <- y
  }
return(keepy)
 }
}


#This code is used to generate n random numbers from the COM-Poisson-Skellam distribution, with the parameter
# values mu1, mu2 and mu3
rcomskellam<-function(n,mu1,mu2,mu3){
makeCMPdata.gen(n,mu1,mu3)-makeCMPdata.gen(n,mu2,mu3)
}