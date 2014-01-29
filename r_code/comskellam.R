comskellam<-function (lmu1 = "loge", lmu2 = "loge", lmu3="loge", imu1 = NULL, imu2 = NULL, imu3= NULL,
    nsimEIM = 100, parallel = FALSE, zero = NULL) {
#return the unevaluated expression, the argument name itself.  
    lmu1 <- as.list(substitute(lmu1))
#result in an empty list, with the attribute of family as the characterized form of the variable
#in the bracket
    emu1 <- link2list(lmu1)
    lmu1 <- attr(emu1, "function.name")
    lmu2 <- as.list(substitute(lmu2))
    emu2 <- link2list(lmu2)
    lmu2 <- attr(emu2, "function.name")
#note here I use mu3 (and lmu3, emu3, etc) ranther than nu or mu(and lnu, enu, etc) as the common name
#of the dispersion parameter in the two data groups.
    lmu3 <- as.list(substitute(lmu3))
    emu3 <- link2list(lmu3)
    lmu3 <- attr(emu3, "function.name")
#To check that the input initial values for the parameters are positive
    if (length(imu1) && !is.Numeric(imu1, positive = TRUE)) 
        stop("bad input for argument 'imu1'")
    if (length(imu2) && !is.Numeric(imu2, positive = TRUE)) 
        stop("bad input for argument 'imu2'")
    if (length(imu3) && !is.Numeric(imu3, positive = TRUE)) 
        stop("bad input for argument 'imu3'")
#To check that the simulation time used to produce the approximate expected information 
#matrix(EIM), namely Wi in our regression algorithm, is integer valued and no less than 500.
    if (!is.Numeric(nsimEIM, allowable.length = 1, integer.valued = TRUE) || 
        nsimEIM <= 50) 
        stop("argument 'nsimEIM' should be an integer greater than 500")
#Creat new family functions in VGAM package and give the description of some properties 
#blurb gives a small description of the model.  
    new("vglmff", blurb = c("Conway-Maxwell-Skellam distribution\n\n", "Links:    ", 
#Only valid when fit this family function into the call "vglm"
#here namesof results s string of character to express the link function, don't be confused by 
#link="identity" as default value, we can imput the link function lmu1 here. Here lmu1="loge"
#is the link function (but not inputed or implemented in the above codes!, refer to "vglm") and
#"mu1" is the paramater after/to fit in the link. 
       namesof("mu1", lmu1, earg = emu1, tag = FALSE), ", ", 
       namesof("mu2", lmu2, earg = emu2, tag = FALSE), ", ", 
       namesof("mu3", lmu3, earg = emu3, tag = FALSE), ", ", "\n", 
       "Mean:     mu1^(1/mu3)-mu2^(1/mu3)"),
#constraints belongs to the object of class "expression" which sets up any constraint matrices defined 
#by arguments in the family function, or more intuitively, a name list of constraint matrices used in 
#the fitting.
       constraints = eval(substitute(expression({
#I simply writes the code here by analogizing to the Skellam Regression in "vglm" algorithm. cm.vgam 
#and cm.zero.vgam are undocumented in "VGAM" package and can not find any explaination in google. However, 
#I think there may be some explanations in the website of a tutorial, I will check if necessary, here I 
#think it is not a problem.
       constraints = cm.vgam(matrix(1, M, 1), x, .parallel, constraints, intercept.apply = TRUE)
       constraints = cm.zero.vgam(constraints, x, .zero, M)
    }), list(.parallel = parallel, .zero = zero))),
#Object of class "expression" used to perform error checking (especially for the variable y) and obtain starting 
#values for the model. In general, etastart or mustart are assigned values based on the variables y, x and w.
        initialize = eval(substitute(expression({
#written to check the integrity of prior weight w and response
       temp5 <- w.y.check(w = w, y = y, ncol.w.max = 1, ncol.y.max = 1, 
            Is.integer.y = TRUE, out.wy = TRUE, maximize = TRUE)
#weight here is an optional vector or matrix of (prior) weights to be used in the fitting process, with respect
#to the number of responses in vector generalized linear regression. if no weight is specified, then it corrspond
#to a vector with 1s as weights.
       w <- temp5$w
        y <- temp5$y
#the expression/form of the predictor link function
        predictors.names <- c(namesof("mu1", .lmu1, earg = .emu1, 
            tag = FALSE), namesof("mu2", .lmu2, earg = .emu2, 
            tag = FALSE), namesof("mu3", .lmu3, earg= .emu3,tag=FALSE))
#if no initial value is set at the begining, use the following method to find the initial value vector for the response 
#variable to begin the iteration.
        if (!length(etastart)) {
            junk = lm.wfit(x = x, y = c(y), w = c(w))
            var.y.est = sum(c(w) * junk$resid^2)/junk$df.residual
            mean.init = weighted.mean(y, w)
            mu1.init = max((var.y.est + mean.init)/2, 0.01)
            mu2.init = max((var.y.est - mean.init)/2, 0.01)
            mu3.init = 1
            mu1.init = rep(if (length(.imu1)) .imu1 else mu1.init, 
                length = n)
            mu2.init = rep(if (length(.imu2)) .imu2 else mu2.init, 
                length = n)
            mu3.init = rep(if (length(.imu3)) .imu3 else mu3.init,
                length = n)
#what follows etastart is the starting value for the response vector, by fitting the initial value of the response variable
#into the link function. 
#function "theta2eta" is used to check the consistency of the link function in the second argument and  
#the attribute of "function.name" in the third argument. then fit the mui.init vector into the link function to get 
#the initial value of the response vector.  
           etastart = cbind(theta2eta(mu1.init, .lmu1, earg = .emu1), 
                theta2eta(mu2.init, .lmu2, earg = .emu2), theta2eta(mu3.init, .lmu3, earg= .emu3))
        }
  }), list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3 = lmu3, .imu1 = imu1, .imu2 = imu2, .imu3 = imu3, 
        .emu1 = emu1, .emu2 = emu2, .emu3 = emu3))), 
#eta is a n*M matrix of linear preditors
#Linkinv: Object of class "function" which returns the fitted values, given the linear/additive predictors. The function
#must have arguments function(eta, extra = NULL). Attention: Don't be confused by the name "linkINV", it denotes that to 
#estimate the dependent variable, we need to calculate the inverse of the link function.
     linkinv = eval(substitute(function(eta, extra = NULL) {
#eta2theta is used to get the inverse of the link function and then apply the "eta" value into the inverse function, to get the 
#value of the response vector in the iteration process.
        mu1 = eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
        mu2 = eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
        mu3 = eta2theta(eta[, 3], link = .lmu3, earg = .emu3)
#the mean difference between the two goups. Don't know its use here. analogy to that in "Skellam" family function.       
        mu1^(1/mu3)-mu2^(1/mu3)
    }, list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3 = lmu3, .emu1 = emu1, .emu2 = emu2, .emu3 = emu3))), 
#Object of class "expression" to insert code at a special position (at the very end) of vglm.fit or vgam.fit. This code is 
#evaluated after the fitting. The list misc is often assigned components in this slot, which becomes the misc slot on the 
#fitted object.
            last = eval(substitute(expression({
            misc$link <- c(mu1 = .lmu1, mu2 = .lmu2, mu3 = .lmu3)
            misc$earg <- list(mu1 = .emu1, mu2 = .emu2, mu3= .emu3)
#Using the approximated expected information matrix (EIM) and indicate the simulation times.
            misc$expected = TRUE
            misc$nsimEIM = .nsimEIM
        }), list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3 = lmu3, .emu1 = emu1, .emu2 = emu2, .emu3 =emu3,
            .nsimEIM = nsimEIM))),
#return the log-likelihood of the model, mu is a n row matrix of fitted values. BUG APPEARS HERE!
            loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL) {
#eta2theta is used to get the inverse of the link function and then apply the "eta" value into the inverse function, to get the 
#value of the response vector in the iteration process.           
            mu1 = eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
            mu2 = eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
            mu3 = eta2theta(eta[, 3], link = .lmu3, earg = .emu3)
            if (residuals) stop("loglikelihood residuals not implemented yet") else {
#note here, different from the "skellam" family function, the "parallel" condition is not taken into consideration (no need to do that). 
            sum(c(w)* (-com.compute.log.vecz(mu1,mu3)-com.compute.log.vecz(mu2,mu3)+ 0.5 * y * log(mu1) - 0.5 * y * 
            log(mu2)+com.compute.log.bessel(2*sqrt(mu1*mu2),y,mu3)))
            }
        }, list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3= lmu3, .parallel = parallel, 
            .emu1 = emu1, .emu2 = emu2, .emu3 = emu3))),
#Object of class "character" giving class information about the family function.
            vfamily = c("Conway-Maxwell-Skellam"), 
#Object of class "expression" which returns a M-column matrix of first derivatives of the log-likelihood function with respect 
#to the linear/additive predictors, equivalently, the score vector.
            deriv = eval(substitute(expression({
#See previous
            mu1 = eta2theta(eta[, 1], link = .lmu1, earg = .emu1)
            mu2 = eta2theta(eta[, 2], link = .lmu2, earg = .emu2)
            mu3 = eta2theta(eta[, 3], link = .lmu3, earg = .emu3)
#First check the consistency of the second and third argument, the same as that is theta2eta, then take the first 
#order derivative of the link and gets its inverse.          
            dmu1.deta = dtheta.deta(mu1, link = .lmu1, earg = .emu1)
            dmu2.deta = dtheta.deta(mu2, link = .lmu2, earg = .emu2)
            dmu3.deta = dtheta.deta(mu3, link = .lmu3, earg = .emu3)
#Get the score vector by taking the derivatives with respect to the response variable.
                temp7 = (y/2-expectation(mu1,mu3))*(1/mu1)
                temp8 = (-y/2-expectation(mu2,mu3))*(1/mu2)
                temp9 = com.compute.derivbessel(sqrt(mu1*mu2),y,mu3)/com.compute.bessel(2*sqrt(mu1*mu2),y,mu3)
                temp6 = logexpectation(mu1,mu3)+logexpectation(mu2,mu3)
                temp10 = com.compute.derivvbessel(sqrt(mu1*mu2),y,mu3)/com.compute.bessel(2*sqrt(mu1*mu2),y,mu3)
                dl.dmu1 = temp7+0.5*sqrt(mu2/mu1)*temp9
                dl.dmu2 = temp8+0.5*sqrt(mu1/mu2)*temp9
                dl.dmu3 = temp6 + temp10
#Get the weighted vector score by taking the derivative with respect to the link function.           
          c(w) * cbind(dl.dmu1 * dmu1.deta, dl.dmu2 * dmu2.deta, dl.dmu3*dmu3.deta)
        }), list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3 = lmu3, .emu1 = emu1, .emu2 = emu2, .emu3 = emu3,
            .nsimEIM = nsimEIM))), 
#weight returns the second derivatives of the log-likelihood function. Here is the approximate expected information matrix (EIM) in 
#simulation method.
            weight = eval(substitute(expression({
            run.var = run.cov12 = run.cov23 =run.cov13 = 0
            for (ii in 1:(.nsimEIM)) {
#generate n (same with the sample size) random numbers belong to the COM-Skellam distribution with corresponding parameter values.
                ysim = rcomskellam(n, mu1 = mu1, mu2 = mu2, mu3=mu3)
                temp7 = (ysim/2-expectation(mu1,mu3))*(1/mu1)
                temp8 = (-ysim/2-expectation(mu2,mu3))*(1/mu2)
                temp9 = com.compute.derivbessel(sqrt(mu1*mu2),ysim,mu3)/com.compute.bessel(2*sqrt(mu1*mu2),ysim,mu3)
                temp6 = logexpectation(mu1,mu3)+logexpectation(mu2,mu3)
                temp10 = com.compute.derivvbessel(sqrt(mu1*mu2),ysim,mu3)/com.compute.bessel(2*sqrt(mu1*mu2),ysim,mu3)
                dl.dmu1 = temp7+0.5*sqrt(mu2/mu1)*temp9
                dl.dmu2 = temp8+0.5*sqrt(mu1/mu2)*temp9
                dl.dmu3 = temp6 + temp10
#remove object                
                rm(ysim)
#simulate the variance-covariance matrix, note that there is a unbiased degree modification in the process.
                temp3 = cbind(dl.dmu1, dl.dmu2, dl.dmu3)
                run.var = ((ii - 1) * run.var + temp3^2)/ii
                run.cov12 = ((ii - 1) * run.cov12 + temp3[, 1] * 
                  temp3[, 2])/ii
                run.cov23 = ((ii - 1) * run.cov23 + temp3[, 2] * 
                  temp3[, 3])/ii
                run.cov13 = ((ii - 1) * run.cov13 + temp3[, 1] * 
                  temp3[, 3])/ii                                
            }
#Here wz is the weigh Wi, a M*M matrix in the algorithm derivation.
            wz = if (intercept.only) matrix(colMeans(cbind(run.var, 
                run.cov12, run.cov23, run.cov13)), n, dimm(M), byrow = TRUE) else cbind(run.var, 
                run.cov12, run.cov23, run.cov13)
            dtheta.detas = cbind(dmu1.deta, dmu2.deta, dmu3.deta)
#Suppose we have n symmetric positive-definite square matrices, each M by M, and these are stored in an array of dimension c(n,M,M). 
#Then these can be more compactly represented by a matrix of dimension c(n,K) where K is an integer between M and M*(M+1)/2 inclusive. 
#The mapping between these two representations is given by this function. It firstly enumerates by the diagonal elements, followed by
#the band immediately above the diagonal, then the band above that one, etc. The last element is (1,M). This function performs the 
#mapping from elements (j,k) of symmetric positive-definite square matrices to the columns of another matrix representing such. This
#is called the matrix-band format and is used by the VGAM package.
            index0 = iam(NA, NA, M = M, both = TRUE, diag = TRUE)
            wz = wz * dtheta.detas[, index0$row] * dtheta.detas[, 
                index0$col]
#the weigted (refer to the prior) weight matrix
            c(w) * wz
        }), list(.lmu1 = lmu1, .lmu2 = lmu2, .lmu3 = lmu3, .emu1 = emu1, .emu2 = emu2, .emu3 = emu3,
            .nsimEIM = nsimEIM))))
}
