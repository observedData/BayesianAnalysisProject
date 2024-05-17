library(latex2exp)
library(haven)
library(tidyverse)
library(tableone) 
library(dplyr)
library(LaplacesDemon)

#Data cleaning process. We restrict our study population to be the PRAMS participants in the 7 
#states that participated in the OCBS. PRAMS survey weight was not used so generalizability to the
#total population is limited. 

data <- read_sas("/Users/linqing/Documents/Bayes/Final Project/data/bayesdataset.sas7bdat")

#We only include individuals who have participated in OCBS in the study population
#We filter by the seven states: "MA","PA","KY","LA","MO","UT", and "WV". Previously, we labeled them as 
#rawSample=1. Please refer to the SAS code in the Supplemental files. 
#After the filter, the total sample size is 8340. 
finalData<-data%>%
  select(pregDepression,postPartumTreat2Week)%>%
  filter(data$rawSample==1)

#We have created a binary variable pregDepression that stores the response of whether the individual have 
#depression during pregnancy. Please refer to the SAS code for details 

#We also created a binary variable called postPartumTreat2Week that stores the response of whether the individual 
#used opioids for two weeks postpartum. 

#we create a binary variable called m1 that stores whether there is missingness for outcome Y. When m1=1 means no missing, and vice versa 
#m2 denotes the missingness in our independent variable. m2=1 means no missing and m2=0 means missing in X. 


#We then create four cases to denote individuals with missingness in both variables, only missing x, only missing y, and no missing.
case4<-finalData%>%
  filter(is.na(finalData$postPartumTreat2Week) & is.na(finalData$pregDepression))%>%
  mutate(m1=0,m2=0)

case3<-finalData%>%
  filter(!is.na(finalData$postPartumTreat2Week) & is.na(finalData$pregDepression))%>%
  mutate(m1=1,m2=0)

case2<-finalData%>%
  filter(is.na(finalData$postPartumTreat2Week) & !is.na(finalData$pregDepression))%>%
  mutate(m1=0,m2=1)

case1<-finalData%>%
  filter(!is.na(finalData$postPartumTreat2Week) & !is.na(finalData$pregDepression))%>%
  mutate(m1=1,m2=1)

#We can combine the four cases back into one dataset 
finalCase<-rbind(case1, case2,case3,case4)

d<-finalCase

colnames(d)
#We change column names to make our life easier 
colnames(d) <- c('x', 'y', 'm1','m2')


#store the total number of data points in each cases
n_xyobs<-nrow(case1)
n_ymiss<-nrow(case2)
n_xmiss<-nrow(case3)
n_xymiss<-nrow(case4)

#Summary statistics 
#Prevalence of individuals who had prenatal depression over observed X data 
depressionNum<-length(which(d$x== 1 & d$m2== 1))#1352
prevDepress=length(which(d$x== 1 & d$m2== 1))/length(which(d$m2== 1)) #16.9%

#Prevalence of individuals who had persistent opioid use over observed Y data 
opioidNum<-length(which(d$y== 1 & d$m1== 1))#  231
prevOpioid=length(which(d$y== 1 & d$m1== 1))/length(which(d$m1== 1)) #11.9% 

#Prevalence of missing in only y outcome over the entire dataset 
prevYMiss<-n_ymiss/nrow(d) #75.3% 6281
#Prevalence of missing in only x outcome over the entire dataset 
prevXMiss<-n_xmiss/nrow(d)#2.5% 210
#Prevalence of missing in both x and y outcomes over the entire dataset 
prevXYMiss<-n_xymiss/nrow(d)#1.3% 111


#Set working directory. 
setwd("/Users/linqing/Documents/Bayes/Final Project/code")
#call our helper methods. 
source("please_work.R")

#initialize x and y missing values 

set.seed(1)

intXmiss<-rbern(n_xmiss+n_xymiss,0.49)

intYmiss<-rbern(n_ymiss+n_xymiss,0.49)


intXmiss1<-rbern(n_xmiss+n_xymiss,0.55)

intYmiss1<-rbern(n_ymiss+n_xymiss,0.55)


intXmiss2<-rbern(n_xmiss+n_xymiss,0.53)

intYmiss2<-rbern(n_ymiss+n_xymiss,0.53)

#burn-in phase: 0-200
#tunning phase: 201-2001
#sampling phase: 2001-5000
omega<-run_sampler(d,intXmiss,intYmiss,60000)
omega1<-run_sampler(d,intXmiss1,intYmiss1,60000)
omega2<-run_sampler(d,intXmiss2,intYmiss2,60000)
#after tunning, we then sample the parameters 
chain1=coda::mcmc(omega[,40001:60000]) 
chain2=coda::mcmc(omega1[,40001:60000]) 
chain3=coda::mcmc(omega2[,40001:60000]) 
class(chain1)
class(chain2)
class(chain3)
all_chains = coda::mcmc.list(chain1,chain2,chain3)
#eps0
chain1eps0=coda::mcmc(chain1[1,1:20000]) 
chain2eps0=coda::mcmc(chain2[1,1:20000]) 
chain3eps0=coda::mcmc(chain3[1,1:20000]) 
class(chain1eps0) 
class(chain2eps0) 
class(chain3eps0) 
all_chainseps0 = coda::mcmc.list(chain1eps0,chain2eps0,chain3eps0)
plot(all_chainseps0,main=TeX('$\\xi_0$'))

summary(all_chainseps0)
#posterior mean eps0 0.4836 (0.4493,0.5173)
coda::gelman.diag(all_chainseps0) #1

lapply(all_chainseps0, effectiveSize)

acfplot(all_chainseps0, ylim=c(-0.2,1.1))


#eps1
chain1eps1=coda::mcmc(chain1[2,1:20000]) 
chain2eps1=coda::mcmc(chain2[2,1:20000]) 
chain3eps1=coda::mcmc(chain3[2,1:20000]) 
class(chain1eps1) 
class(chain2eps1)
class(chain3eps1)
all_chainseps1 = coda::mcmc.list(chain1eps1, chain2eps1, chain3eps1)
plot(all_chainseps1,main=TeX(' $\\xi_1$'))

summary(all_chainseps1)
# -0.7365671 (-0.8655,  -0.6049  )
coda::gelman.diag(all_chainseps1) 1.01

lapply(all_chainseps1, effectiveSize) 

acfplot(all_chainseps1, ylim=c(-0.2,1.1))


#beta0
chain1beta0=coda::mcmc(chain1[3,1:20000]) 
chain2beta0=coda::mcmc(chain2[3,1:20000]) 
chain3beta0=coda::mcmc(chain3[3,1:20000]) 
class(chain1beta0) 
class(chain2beta0) 
class(chain3beta0) 
all_chainsbeta0 = coda::mcmc.list(chain1beta0, chain2beta0, chain3beta0)
plot(all_chainsbeta0,main=TeX(' $\\beta_0$'))

summary(all_chainsbeta0)
#mean -2.7345050, (-2.901, -2.564) 
coda::gelman.diag(all_chainsbeta0) 1

lapply(all_chainsbeta0, effectiveSize)

acfplot(all_chainsbeta0, ylim=c(-0.2,1.1))


chain1beta1=coda::mcmc(chain1[4,1:20000]) 
chain2beta1=coda::mcmc(chain2[4,1:20000]) 
chain3beta1=coda::mcmc(chain3[4,1:20000]) 
class(chain1beta1) 
class(chain2beta1) 
class(chain3beta1) 
all_chainsbeta1 = coda::mcmc.list(chain1beta1, chain2beta1, chain3beta1)
plot(all_chainsbeta1,main=TeX('$\\beta_1$'))

summary(all_chainsbeta1)
# 0.2236909 (-0.01056, 0.44829 )
coda::gelman.diag(all_chainsbeta1)#1.01

lapply(all_chainsbeta1, effectiveSize)

acfplot(all_chainsbeta1, ylim=c(-0.2,1.1))


chain1theta=coda::mcmc(chain1[5,1:20000]) 
chain2theta=coda::mcmc(chain2[5,1:20000]) 
chain3theta=coda::mcmc(chain3[5,1:20000]) 
class(chain1theta) 
class(chain2theta) 
class(chain3theta) 
all_chainstheta = coda::mcmc.list(chain1theta, chain2theta, chain3theta)
plot(all_chainstheta,main=TeX('Traceplots for $\\theta$'))
summary(all_chainstheta)
# 0.1.688 (0.1603, 0.1767)
coda::gelman.diag(all_chainstheta) #1 

lapply(all_chainstheta, effectiveSize)

acfplot(all_chainstheta, ylim=c(-0.2,1.1))

#compute posterior fit check 
#observed opioid use among prenatal depression
prevOpioidDepressDO=length(which(d$y== 1 & d$m1== 1 & d$m2== 1 & d$x==1))/length(which(d$x== 1 & d$m1== 1 & d$m2== 1)) #54 16.9%
#observed opioid use among no prenatal depression
prevOpioidNoDepressDO=length(which(d$y== 1 & d$m1== 1 & d$m2== 1 & d$x==0))/length(which(d$x== 0 & d$m1== 1 & d$m2== 1)) #157 11.1%

#imputed data 
#prevalence of persistent opioid use among individuals who had prenatal depression 
chain1prevOpioidDepressDO=coda::mcmc(chain1[6,1:20000]) 
chain2prevOpioidDepressDO=coda::mcmc(chain2[6,1:20000]) 
chain3prevOpioidDepressDO=coda::mcmc(chain3[6,1:20000]) 
class(chain1prevOpioidDepressDO) 
class(chain2prevOpioidDepressDO) 
class(chain3prevOpioidDepressDO) 
all_chainsprevOpioidDepressDO = coda::mcmc.list(chain1prevOpioidDepressDO, chain2prevOpioidDepressDO, chain3prevOpioidDepressDO)
summary(all_chainsprevOpioidDepressDO) #0.1693

#prevalence of persistent opioid use among individuals who did not have prenatal depression 
chain1prevOpioidNoDepressDO=coda::mcmc(chain1[7,1:20000]) 
chain2prevOpioidNoDepressDO=coda::mcmc(chain2[7,1:20000]) 
chain3prevOpioidNoDepressDO=coda::mcmc(chain3[7,1:20000]) 
class(chain1prevOpioidNoDepressDO) 
class(chain2prevOpioidNoDepressDO) 
class(chain3prevOpioidNoDepressDO) 
all_chainsprevOpioidNoDepressDO = coda::mcmc.list(chain1prevOpioidNoDepressDO, chain2prevOpioidNoDepressDO, chain3prevOpioidNoDepressDO)
summary(all_chainsprevOpioidNoDepressDO) #0.1106

#number of individuals who had persistent opioid use and had prenatal depression 
chain1noOpioid=coda::mcmc(chain1[8,1:20000]) 
chain2noOpioid=coda::mcmc(chain2[8,1:20000]) 
chain3noOpioid=coda::mcmc(chain3[8,1:20000]) 
class(chain1noOpioid) 
class(chain2noOpioid) 
class(chain3noOpioid) 
all_chainsnoOpioid = coda::mcmc.list(chain1noOpioid, chain2noOpioid, chain3noOpioid)
summary(all_chainsnoOpioid)

#number of individuals who had persistent opioid use and had no prenatal depression 
chain1noDepress=coda::mcmc(chain1[9,1:20000]) 
chain2noDepress=coda::mcmc(chain2[9,1:20000]) 
chain3noDepress=coda::mcmc(chain3[9,1:20000]) 
class(chain1noDepress) 
class(chain2noDepress) 
class(chain3noDepress) 
all_chainsnoDepress = coda::mcmc.list(chain1noDepress, chain2noDepress, chain3noDepress)
summary(all_chainsnoDepress)

####Second run 

nomega<-run_sampler(d,intXmiss,intYmiss,80000)
nomega1<-run_sampler(d,intXmiss1,intYmiss1,80000)
nomega2<-run_sampler(d,intXmiss2,intYmiss2,80000)
#pilot runs to set the step size
nchain1=coda::mcmc(nomega[,40001:80000]) 
nchain2=coda::mcmc(nomega1[,40001:80000]) 
nchain3=coda::mcmc(nomega2[,40001:80000]) 

nchain1eps0=coda::mcmc(nchain1[1,1:20000]) 
nchain2eps0=coda::mcmc(nchain2[1,1:20000]) 
nchain3eps0=coda::mcmc(nchain3[1,1:20000]) 
class(nchain1eps0) 
class(nchain2eps0) 
class(nchain3eps0) 
all_nchainseps0 = coda::mcmc.list(nchain1eps0,nchain2eps0,nchain3eps0)
plot(all_nchainseps0,main=TeX('$\\xi_0$'))

summary(all_nchainseps0)
#posterior mean eps0 0.330 (-0.027, 0.514)
coda::gelman.diag(all_nchainseps0) #1

lapply(all_nchainseps0, effectiveSize)

acfplot(all_nchainseps0, ylim=c(-0.2,1.1))


#eps1
nchain1eps1=coda::mcmc(nchain1[2,1:20000]) 
nchain2eps1=coda::mcmc(nchain2[2,1:20000]) 
nchain3eps1=coda::mcmc(nchain3[2,1:20000]) 
class(nchain1eps1) 
class(nchain2eps1)
class(nchain3eps1)
all_nchainseps1 = coda::mcmc.list(nchain1eps1, nchain2eps1, nchain3eps1)
plot(all_nchainseps1,main=TeX(' $\\xi_1$'))

summary(all_nchainseps1)
#  -0.0491 (-0.856 ,  1.441 )
coda::gelman.diag(all_nchainseps1)# 

lapply(all_nchainseps1, effectiveSize) 

acfplot(all_nchainseps1, ylim=c(-0.2,1.1))

#beta0
chain1nbeta0=coda::mcmc(nchain1[3,1:20000]) 
chain2nbeta0=coda::mcmc(nchain2[3,1:20000]) 
chain3nbeta0=coda::mcmc(nchain3[3,1:20000]) 
class(chain1nbeta0) 
class(chain2nbeta0) 
class(chain3nbeta0) 
all_chainsnbeta0 = coda::mcmc.list(chain1nbeta0, chain2nbeta0, chain3nbeta0)
plot(all_chainsnbeta0,main=TeX(' $\\beta_0$'))

summary(all_chainsnbeta0)
#mean -1.870, (-2.896, -0.002) 

lapply(all_chainsnbeta0, effectiveSize)

acfplot(all_chainsnbeta0, ylim=c(-0.2,1.1))


nchain1beta1=coda::mcmc(nchain1[4,1:20000]) 
nchain2beta1=coda::mcmc(nchain2[4,1:20000]) 
nchain3beta1=coda::mcmc(nchain3[4,1:20000]) 
class(nchain1beta1) 
class(nchain2beta1) 
class(nchain3beta1) 
nall_chainsbeta1 = coda::mcmc.list(nchain1beta1, nchain2beta1, nchain3beta1)
plot(nall_chainsbeta1,main=TeX('$\\beta_1$'))

summary(nall_chainsbeta1)
# 0.152 (-0.10, 0.432 )

lapply(nall_chainsbeta1, effectiveSize)

acfplot(nall_chainsbeta1, ylim=c(-0.2,1.1))


nchain1theta=coda::mcmc(nchain1[5,1:20000]) 
nchain2theta=coda::mcmc(nchain2[5,1:20000]) 
nchain3theta=coda::mcmc(nchain3[5,1:20000]) 
class(nchain1theta) 
class(nchain2theta) 
class(nchain3theta) 
all_nchainstheta = coda::mcmc.list(nchain1theta, nchain2theta, nchain3theta)
plot(all_nchainstheta,main=TeX('$\\theta$'))
summary(all_nchainstheta)
# 0.1.688 (0.160, 0.177)

lapply(all_nchainstheta, effectiveSize)

acfplot(all_nchainstheta, ylim=c(-0.2,1.1))

#compute posterior fit check 
#observed opioid use among prenatal depression
prevOpioidDepressDO=length(which(d$y== 1 & d$m1== 1 & d$m2== 1 & d$x==1))/length(which(d$x== 1 & d$m1== 1 & d$m2== 1)) #54 16.9%
#observed opioid use among no prenatal depression
prevOpioidNoDepressDO=length(which(d$y== 1 & d$m1== 1 & d$m2== 1 & d$x==0))/length(which(d$x== 0 & d$m1== 1 & d$m2== 1)) #157 11.1%

#imputed data 
#prevalence of persistent opioid use among individuals who had prenatal depression 
chain1prevOpioidDepressDO=coda::mcmc(nchain1[6,1:20000]) 
chain2prevOpioidDepressDO=coda::mcmc(nchain2[6,1:20000]) 
chain3prevOpioidDepressDO=coda::mcmc(nchain3[6,1:20000]) 
class(chain1prevOpioidDepressDO) 
class(chain2prevOpioidDepressDO) 
class(chain3prevOpioidDepressDO) 
all_chainsprevOpioidDepressDO = coda::mcmc.list(chain1prevOpioidDepressDO, chain2prevOpioidDepressDO, chain3prevOpioidDepressDO)
summary(all_chainsprevOpioidDepressDO) #0.1693

#prevalence of persistent opioid use among individuals who did not have prenatal depression 
chain1prevOpioidNoDepressDO=coda::mcmc(nchain1[7,1:20000]) 
chain2prevOpioidNoDepressDO=coda::mcmc(nchain2[7,1:20000]) 
chain3prevOpioidNoDepressDO=coda::mcmc(nchain3[7,1:20000]) 
class(chain1prevOpioidNoDepressDO) 
class(chain2prevOpioidNoDepressDO) 
class(chain3prevOpioidNoDepressDO) 
all_chainsprevOpioidNoDepressDO = coda::mcmc.list(chain1prevOpioidNoDepressDO, chain2prevOpioidNoDepressDO, chain3prevOpioidNoDepressDO)
summary(all_chainsprevOpioidNoDepressDO) #0.1106

#number of individuals who had persistent opioid use and had prenatal depression 
chain1noOpioid=coda::mcmc(nchain1[8,1:20000]) 
chain2noOpioid=coda::mcmc(nchain2[8,1:20000]) 
chain3noOpioid=coda::mcmc(nchain3[8,1:20000]) 
class(chain1noOpioid) 
class(chain2noOpioid) 
class(chain3noOpioid) 
all_chainsnoOpioid = coda::mcmc.list(chain1noOpioid, chain2noOpioid, chain3noOpioid)
summary(all_chainsnoOpioid) #54

#number of individuals who had persistent opioid use and had no prenatal depression 
chain1noDepress=coda::mcmc(nchain1[9,1:20000]) 
chain2noDepress=coda::mcmc(nchain2[9,1:20000]) 
chain3noDepress=coda::mcmc(nchain3[9,1:20000]) 
class(chain1noDepress) 
class(chain2noDepress) 
class(chain3noDepress) 
all_chainsnoDepress = coda::mcmc.list(chain1noDepress, chain2noDepress, chain3noDepress)
summary(all_chainsnoDepress) #753
