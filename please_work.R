###--------------------------------------------------------------------------###
### Author: Lucy Zheng
### Helper functions for dealing with missing data
###--------------------------------------------------------------------------###
library(LaplacesDemon)
######## Helper function to update parameters ###### 



#Helper function eps_posterior computes f(eps| m,y,x)the posterior up to a proportionality constant 
#d: data that consists of x,y,m1,m2
eps_posterior=function(d,eps0,eps1){
  #log_lik is the log scale likelihood from Bernoulli distribution 
  prob1=invlogit(eps0 + eps1*d$x)
  prob2=invlogit(eps0 + eps1*d$y)
  #the likelihood is invlogit(eps0 + eps1*d$x)^m1 [1-invlogit(eps0 + eps1*d$x)]^1-m1
  #and invlogit(eps0 + eps1*d$x)^m2 [1-invlogit(eps0 + eps1*d$y)]^1-m2
  log_lik1 = sum(dbern( d$m1, prob = prob1, log = T ))
  log_lik2 = sum(dbern( d$m2, prob = prob2, log = T ))
  #log priors for eps. we want tight priors because we might not have enough data to tune etas. 
  log_eps0=dnorm( eps0 , mean = 0, sd =1, log = T) 
  log_eps1=dnorm( eps1 , mean = 0, sd =1, log = T) 
  log_post_eval  = log_lik1+log_lik2 +log_eps0 + log_eps1
  
  return(log_post_eval)
}

#Helper function eta_posterior computes f(y|xi; eta)f(eta) 
#eta consists of beta0 and beta1
eta_posterior=function(d,beta0,beta1){
  prob=invlogit(beta0 + beta1*d$x)
  log_lik = sum(dbern( d$y, prob = prob, log = T ))
  #log priors for betas
  log_beta0=dnorm(beta0,mean=0,sd=sqrt(3),log=T)
  log_beta1=dnorm(beta1,mean=0,sd=sqrt(3),log=T)
  log_post_eval  = log_lik + log_beta0 +  log_beta1
  
  return(log_post_eval)
}

#Helper function to find the success probability for imputing x when 
#subjects are missing both x and y 
#d=d[d$m1==0 & d$m2==0,]
findxprob2=function(d,eps0,eps1,beta0,beta1,theta){
  prob1=invlogit(eps0 + eps1) #P(M1=0|y,X=1) 
  #because m1 = 0 so we don't need 1-invlogit 
  p_m1x1=dbern(d$m1, prob = prob1, log = T ) #vector with a length of xymiss 
  #This y is from previous draws? 
  #because m2=0 so we don't need 1-invlogit 
  p_m2x1=dbern(d$m2,invlogit(eps0 + eps1*d$y),log=T)#P(M2=0|y,X=1)  #vector with a length of xymiss 
  
  p_yx1=dbern(d$y,invlogit(beta0+beta1),log=T)#f(y|x=1)
  
  #a length d[d$m1==0 & d$m2==0,] vector that has 1 to represent x==1
  xvec1<-rep(1,nrow(d[d$m1==0 & d$m2==0,]))
  p_x1=dbern(xvec1, prob = theta, log = T ) #P(x=1)
  
  x1=p_m1x1+ p_m2x1+ p_yx1+ p_x1 # vector that stores our numerator 
  
  
  p_m1x0=dbern(d$m1, prob = invlogit(eps0), log = T )  #P(M1=0|y,x=0)
  
  p_m2x0=dbern(d$m2,invlogit(eps0 + eps1*d$y),log=T) #P(M2=0|y,X=0) 
  
  p_yx0=dbern(d$y,invlogit(beta0),log=T)#f(y|x=0)
  
  #a length d[d$m1==0 & d$m2==0,] vector that has 0 to represent x==0
  xvec0<-rep(0,nrow(d[d$m1==0 & d$m2==0,]))
  #because xvec0 is zero so we don't need 1-theta
  p_x0=dbern(xvec0, prob = theta, log = T ) #P(x=0)
  
  
  x0=p_m1x0+p_m2x0+p_yx0+p_x0
  
  denominator<-exp(x1)+exp(x0)
  
  successP<-exp(x1-log(denominator))
  
  return(successP)
}


# #### Testing findxprob2 ####
# x<-c(1,0,0,1)
# y<-c(1,1,0,0)
# m1<-c(0,0,0,0)
# m2<-c(0,0,0,0)
# d<-cbind(x,y,m1,m2)
# d<-data.frame(d)
# eps0=0.001
# eps1=0.002
# beta0=0.1
# beta1=0.3
# theta=0.2

#Helper function to find the success probability for imputing x when 
#subjects are missing only x 
#d=d[d$m1==1 & d$m2==0,]
findxprob1=function(d,eps0,eps1,beta0,beta1,theta){
  #m1=1 given x = 1, success prob = invlogit(eps0+eps1) 
  
  p_m1x1=dbern(d$m1, prob = invlogit(eps0 + eps1), log = T )#P(M1=1|y,X=1)
  
  #because all m2==0, we don't need 1-invlogit 
  p_m20x1=dbern(d$m2,prob=invlogit(eps0 + eps1*d$y),log=T)#P(M2=0|y,X=1)
  
  #f(y|x=1)
  py_x1=dbern(d$y,prob=invlogit(beta0+beta1),log=T)
  
  #p(x=1)
  #a length d[d$m1==1 & d$m2==0,] vector that has 1 to represent x==1
  xvec1<-rep(1,nrow(d[d$m1==1 & d$m2==0,]))
  p_x1=dbern(xvec1, prob = theta, log = T ) #P(x=1)
  
  x1<-p_m1x1+p_m20x1+py_x1+ p_x1
  
  #m1=1 given X=0
  p_m1x0=dbern(d$m1, prob = invlogit(eps0), log = T )#P(M1=1|y,X=0)
  
  p_m20x0=dbern(d$m2,prob=invlogit(eps0 + eps1*d$y),log=T)#P(M2=0|y,X=0)
  
  py_x0=dbern(d$y,prob=invlogit(beta0),log=T)#f(y|x=0)
  
  #p(x=0)
  #a length d[d$m1==1 & d$m2==0,] vector that has 0 to represent x==0
  xvec0<-rep(0,nrow(d[d$m1==1 & d$m2==0,]))
  p_x0=dbern(xvec0, prob = theta, log = T ) #P(x=0)
  
  x0<- p_m1x0+p_m20x0+py_x0+ p_x0
  
  denominator<-exp(x1)+exp(x0)
  
  successP<-exp(x1-log(denominator))
  
  return(successP)
  
}



# #### Testing findxprob1 ####
# x<-c(1,0,0,1)
# y<-c(1,1,0,0)
# m1<-c(1,1,1,1)
# m2<-c(0,0,0,0)
# d<-cbind(x,y,m1,m2)
# d<-data.frame(d)
# eps0=0.001
# eps1=0.002
# beta0=0.1
# beta1=0.3
# theta=0.2



#Helper function to find the success probability for imputing y when 
#subjects are missing  x and y 
#d=d[d$m1==0 & d$m2==0,]
findyprob2=function(d,eps0,eps1,beta0,beta1,theta){
  
  #P(M1=0|x) All m1=0 so we don't need 1-invlogit 
  p_m10x<-dbern(d$m1,prob =invlogit(eps0+eps1*d$x),log=T)
  
  #P(M2=0|Y=1) All m2=0 so we don't need 1-invlogit 
  p_m20y1<-dbern(d$m2,prob=invlogit(eps0+eps1),log=T)
  
  #f(Y=1|x)
  #create a vector of length d[d$m1==0 & d$m2==0,] to with 1 in it. 
  yvec1<-rep(1,nrow(d[d$m1==0 & d$m2==0,]))
  
  f_y1<-dbern(yvec1,prob=invlogit(beta0+beta1*d$x),log=T)
  
  #P(X)
  p_x=dbern(d$x, prob = theta, log = T ) #P(x=x)
  
  y1<-p_m10x+p_m20y1+f_y1+p_x
  
  #P(M2=0|Y=0) All m2=0 so we don't need 1-invlogit 
  p_m20y0<-dbern(d$m2,prob=invlogit(eps0),log=T)
  
  #f(Y=0|x)
  #create a vector of length d[d$m1==0 & d$m2==0,] to with 0 in it. 
  yvec0<-rep(0,nrow(d[d$m1==0 & d$m2==0,]))
  
  f_y0<-dbern(yvec0,prob=invlogit(beta0+beta1*d$x),log=T)
  
  y0<- p_m10x+p_m20y0+f_y0+p_x
  
  denominator<-exp(y1)+exp(y0)
  
  successP<-exp(y1-log(denominator))
  
  return(successP)
}



#### Testing findyprob2 ####
# y<-c(1,0,0,1)
# x<-c(1,1,0,0)
# m1<-c(0,0,0,0)
# m2<-c(0,0,0,0)
# d<-cbind(x,y,m1,m2)
# d<-data.frame(d)
# eps0=0.001
# eps1=0.002
# beta0=0.1
# beta1=0.3
# theta=0.2


#Helper function to find the success probability for imputing y when 
#subjects are missing only y 
#d=d[d$m1==0 & d$m2==1,]
findyprob1=function(d,eps0,eps1,beta0,beta1,theta){
  
  #P(M1=0|x) All m1=0 so we don't need 1-invlogit 
  p_m10x<-dbern(d$m1,prob =invlogit(eps0+eps1*d$x),log=T)
  
  #P(M2=1|Y=1) All m2=1 
  p_m21y1<-dbern(d$m2,prob=invlogit(eps0+eps1),log=T)
  
  #f(Y=1|x)
  #create a vector of length d[d$m1==0 & d$m2==1,] with 1 in it. 
  yvec1<-rep(1,nrow(d[d$m1==0 & d$m2==1,]))
  
  f_y1<-dbern(yvec1,prob=invlogit(beta0+beta1*d$x),log=T)
  
  #P(X)
  p_x=dbern(d$x, prob = theta, log = T ) #P(x=x)
  
  y1<-p_m10x+p_m21y1+f_y1+p_x
  
  #P(M2=1|Y=0) 
  p_m21y0<-dbern(d$m2,prob=invlogit(eps0),log=T)
  
  #f(Y=0|x)
  #create a vector of length d[d$m1==0 & d$m2==1,] with 0 in it. 
  yvec0<-rep(0,nrow(d[d$m1==0 & d$m2==1,]))
  
  f_y0<-dbern(yvec0,prob=invlogit(beta0+beta1*d$x),log=T)
  
  y0<- p_m10x+p_m21y0+f_y0+p_x
  
  denominator<-exp(y1)+exp(y0)
  
  successP<-exp(y1-log(denominator))
  
  return(successP)
}


###testing for findyprob1 ###


######## Metropolis-Hastings algorithm + gibbs sampling ########
run_sampler = function(data, x_miss,y_miss,M){
  # Local variables for 5 hyperparameters across M iterations: #eps0,eps1,beta0,beta1,theta
  eps= matrix(nrow=2, ncol = M) 
  #eta stores beta0 and beta1 
  eta=matrix(nrow=2,ncol=M)
  #theta
  thetaDraw=vector(length=M)
  
  #local variable d that stores the dataset for each iteration
  d<-data
  # setting initial values by parsing in missing covariates and missing outcomes. 
  d$x[d$m2==0]<-x_miss 
  d$y[d$m1==0]<-y_miss 
  #number of missing in x but no missing in y 
  n_xmiss=nrow(d[d$m2==0 & d$m1==1,])
  #number of missing in y but no missing in x 
  n_ymiss=nrow(d[d$m1==0 & d$m2==1,])
  #number of missing in y and x 
  n_xymiss=nrow(d[d$m1==0 & d$m2==0,])
  
  phi_eps = 0.003 ## proposal variance for eps 0.001 0.003, - 0.005 
  
  phi_eta= 0.004 # proposal variance for eta 0.004. 
  
  #set initial xi 
  ## propose candidate from independent normal
  eps[, 1]= rnorm(n = 2, mean = 0,sqrt(1) ) #from prior 
  
  #set initial eta
  eta[ ,1]=rnorm(n=2,mean = 0, sqrt(3)) #from prior 
  
  
  #update theta: theta is draw from the posterior from beta conjugacy prior 
  #We set prior to be theta~beta(2,2), therefore our posterior distribution has a closed form beta(sum(x)+2,n-sum(x)+2) 
  thetaDraw[1]=rbeta(1,sum(d$x)+2, nrow(d)-sum(d$x)+2)
  
  for(m in 2:M){
    
    #update xi 
    ## propose candidate from independent normal for eps0 and eps1
    eps_star = rnorm(n = 2, mean = eps[, m-1], sqrt(phi_eps) ) #centered around the previous draw
    
    ## evaluate un-normalized posterior density of candidate
    log_ftilde_eps_star = eps_posterior(d, eps_star[1], eps_star[2])
    ## evaluate un-normalized posterior density of previous value
    log_ftilde_eps_curr = eps_posterior(d, eps[1, m-1], eps[2, m-1] )
    
    ## compute acceptance probability: not proposal not needed (symmetry)
    r = exp( log_ftilde_eps_star - log_ftilde_eps_curr)
    alpha = min(1, r)
    
    if( runif(1) < alpha ){ ## accept with probability alpha
      eps[,m] = eps_star
    }else{
      eps[,m] = eps[,m-1] ## reject with probability 1-alpha
    }
    
    
    #update eta
    ## propose candidate from independent normal for beta0 and beta1
    eta_star = rnorm(n = 2, mean = eta[, m-1], sqrt(phi_eta) ) #centered around the previous draw
    
    ## evaluate un-normalized posterior density of candidate
    log_ftilde_eta_star = eta_posterior(d, eta_star[1], eta_star[2])
    ## evaluate un-normalized posterior density of previous value
    log_ftilde_eta_curr = eta_posterior(d, eta[1, m-1], eta[2, m-1] )
    
    ## compute acceptance probability: not proposal not needed (symmetry)
    r = exp( log_ftilde_eta_star - log_ftilde_eta_curr)
    alpha = min(1, r)
    
    if( runif(1) < alpha ){ ## accept with probability alpha
      eta[,m] = eta_star
    }else{
      eta[,m] = eta[,m-1] ## reject with probability 1-alpha
    }
    
    
    
    #update theta: theta is draw from the posterior from beta conjugacy prior 
    #We set prior to be theta~beta(2,2), therefore our posterior distribution has a closed form beta(sum(x)+2,n-sum(x)+2) 
    thetaDraw[m]=rbeta(1,sum(d$x)+2,nrow(d)-sum(d$x)+2)
    
    
    
    #y miss 
    d$y<-ifelse(d$m1==0,
          yes= ifelse(d$m2==0, #if missing both x and y values 
                       yes=rbern(n=n_xymiss,prob=findyprob2(d[d$m1==0 & d$m2==0,],eps[1,m],eps[2,m],eta[1,m],eta[2,m],thetaDraw[m])), 
                       no=rbern(n=n_ymiss,prob=findyprob1(d[d$m1==0 & d$m2==1,],eps[1,m],eps[2,m],eta[1,m],eta[2,m],thetaDraw[m]))), #if missing only y values 
          no=d$y )# if no missing, keep original value 
    
    
    #x miss
    d$x<-ifelse(d$m2==0, #if missing both x and y values
                yes=ifelse(d$m1==0,yes=rbern(n=n_xymiss,prob=findxprob2(d[d$m1==0 & d$m2==0,],eps0=eps[1,m],eps1=eps[2,m],eta[1,m],eta[2,m],thetaDraw[m])),
                                   no=rbern(n=n_xmiss,prob=findxprob1(d[d$m1==1 & d$m2==0,],eps[1,m],eps[2,m],eta[1,m],eta[2,m],thetaDraw[m]))), # if missing only x values 
                no=d$x) # if no missing, keep original value
  }
  
  #the X and Y missing imputation result for the last iteration for prediction analysis
  #We compute the posterior predictive checks here
  #prevalence of opioid use among individuals who are depressed 
  #percentage
  prevOpioidDepressDO=length(which(d$y== 1 & d$m1== 0 & d$m2== 0 & d$x==1))/length(which(d$x== 1 & d$m1== 0 & d$m2== 0)) 
  #number 
  noOpioid<-length(which(d$y== 1 & d$m1== 0 & d$m2== 0 & d$x==1))
  #prevalence of opioid use among individuals who are not depressed 
  prevOpioidNoDepressDO=length(which(d$y== 1 & d$m1== 0 & d$m2== 0 & d$x==0))/length(which(d$x== 0 & d$m1== 0 & d$m2== 0)) 
  noDepress<-which(d$y== 1 & d$m1== 0 & d$m2== 0 & d$x==0)
  #store all parameters into one matrix 
  omega= matrix(nrow=9, ncol = M) 
  omega<-rbind(eps,eta, thetaDraw, prevOpioidDepressDO,prevOpioidNoDepressDO,noOpioid,noDepress)
  #return all parameters
  return(omega)
  
}

