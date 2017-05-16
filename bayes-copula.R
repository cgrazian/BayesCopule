rm(list=ls())

#"Quadratic" function in Owen (2001) (Section 3.14)

Owen <- function(x, thre, gradiente)
{
  grad <- t(gradiente)
  par <- t(as.matrix(x))
  eps <- 1/thre
  
  z <- 1+par%*%grad
  ans <- z
  lo <- ((z)<eps)
  
  ans[ lo  ] <- log(eps) - 1.5 + 2*(z[lo])/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  -sum(ans)
}

Schennach <- function(x, gradiente)
{
  grad <- t(gradiente)
  par <- t(as.matrix(x))
  
  z <- exp(par%*%grad)/length(gradiente)
  sum(z)
}

#Below the function for the computation of empirical likelihood

##input: 
####### 1) SC --> individual contributions of the estimating function. 
####### It can be either an n x 1 vector (scalar parameter) or an n x p matrix 
####### (vector valued parameter)

##output: a list with elements
######## 1) w0 --> weights associated to each unit
######## 2) conv --> convergence of the algorithm (available only for vector valued parameter)
######## 3) elr --> empirical loglikelihood ratio
######## 4) el --> empirical loglikelihood

EL <- function(SC)
{
  n <- NROW(SC)
  p <- NCOL(SC)
  
  ##find Lagrange multiplier
  if(p!=1)
  {
    OBJ <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, control=list(maxit=1e4))
    molt_lagr_0 <- OBJ$pa
    conv <- OBJ$conv
  }
  else
  {
    molt_lagr_0 <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, method="Brent", lower=-1e1, upper=1e1)$pa
  }
  
  ##weights
  w_emp_0 <- as.numeric( 1/(n*(1+molt_lagr_0%*%t(SC))) )
  
  if(p!=1)
  {
    list(w0=w_emp_0, conv=conv, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0))
  }
  else
  {
    list(w0=w_emp_0, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0) )
  }
}

EL.Schennach <- function(SC)
{
  n <- NROW(SC)
  p <- NCOL(SC)
  
  ##find Lagrange multiplier

	molt_lagr_0 <- optim(par=rep(0,p), fn=Schennach, gradiente=SC, 
					method="Brent", lower=-1e1, upper=1e1)$pa
 
 	##weights
	const=as.numeric(sum(exp(molt_lagr_0%*%t(SC))))
	w_emp_0 <- as.numeric( exp(molt_lagr_0%*%t(SC)) / const )
	el=log(w_emp_0)
	elr=-2*sum(log(w_emp_0))-2*n*log(n)
  
	return( list(w=w_emp_0,el=el,elr=elr ) )
	
}

### Population mean example

n <- 30
mu <- 1
set.seed(1)
data <- rnorm(n, mu)
est.func <- data-mu
EL(est.func)
EL.Schennach(est.func)

###########################
### ABC for copula models
###########################

install.packages("gmm")
require(gmm)

install.packages("copula")
require(copula)

### Data simulation

theta=1.076 # Spearman =.5
ssize=100
S=1000  #MC sample size

cc<-claytonCopula(theta)
uu<-rCopula(ssize, cc)

### ELCOP function

ELCOP<- function(sample){
	n<-dim(sample)[1] 
	p<-dim(sample)[2]
	rr<-rank(sample[,1])
	ss<-rank(sample[,2])
	vv<-rr*ss
	#print(vv)
	cor(uu[,1], uu[,2], method="spearman")
	S=1000
	rh<-runif(S, -1,1)
	omega<-rep(0,S)
	for (s in 1:S) {
		  #print(rh[s])
		  estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
		  omega[s]<-exp(-EL(estim)$elr)
 		 #print(omega[s])
	}
	psam<-sample(rh, size=S, rep=T, prob=omega)
	sintesi<-c(quantile(psam,.05), quantile(psam,.5),  quantile(psam,.95))
	return(as.numeric(sintesi))  
}

ELCOP.Schennach<- function(sample){
	n<-dim(sample)[1] 
	p<-dim(sample)[2]
	rr<-rank(sample[,1])
	ss<-rank(sample[,2])
	vv<-rr*ss
	#print(vv)
	cor(uu[,1], uu[,2], method="spearman")
	S=1000
	rh<-runif(S, -1,1)
	omega<-rep(0,S)
	for (s in 1:S) {
		  #print(rh[s])
		  estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
		  omega[s]<-exp(-EL.Schennach(estim)$elr)
 		 #print(omega[s])
	}
	psam<-sample(rh, size=S, rep=T, prob=omega)
	sintesi<-c(quantile(psam,.05), quantile(psam,.5),  quantile(psam,.95))
	return(as.numeric(sintesi))  
}


as.numeric(ELCOP(uu))
as.numeric(ELCOP.Schennach(uu))

### Repeat simulations for 100 samples

require(gmm)
require(copula)

sim.cop.clay<-function(theta=1.076, Z=10,ssize=100, 
S=1000){
# ssize= Sample Size
# S=  MC sample size
cc<-claytonCopula(theta)
sintesi.finale<-matrix(rep(0, 3*Z), nrow=Z)
#
for (i in 1:Z){
  uu<-rCopula(ssize, cc)
  sintesi.finale[i,]<-ELCOP(uu)
  }
boxplot(sintesi.finale, xlab="5th, 50th and 95th percentiles of the posterior of rho")
title("Boxplots: Clayton Copula: true Rho=0.5")
return(sintesi.finale)
}

system.time(prova.clay1<-sim.cop.clay(ssize=50, Z=500, S=1000))

save.image("copula-Clayton.RData")

### Comparison between ABC and classical methods

ELCOP.vs.freq<- function(sample,S=1000){

	n<-dim(sample)[1] 
	p<-dim(sample)[2]

	# Ranghi
	rr<-rank(sample[,1])
	ss<-rank(sample[,2])
	vv<-rr*ss

	# Frequentist estimate
	rhon=12 * sum(vv) /(n*(n+1)*(n-1)) - 3*(n+1)/(n-1)
	
	# Frequentist estimate of the variance
	An=sum(rr/(n+1)*ss/(n+1))/n
	Bn=sum((rr/(n+1))^2*(ss/(n+1))^2)/n

	Cn=0
	for(i in 1:length(rr))
	{
		const=rr[i]/(n+1)*ss[i]/(n+1)
		for(j in 1:length(rr))
		{
			for(k in 1:length(33))
			{
				Cn=Cn+const*(rr[k]<=rr[i])*(ss[k]<=ss[j])				
			}
		}
	}	
	Cn=Cn/n^3+0.25-An

	Dn=0
	for(i in 1:length(rr))
	{
		for(j in 1:length(ss))
		{
			Dn=Dn+ss[i]*ss[j]/(n+1)^2*max(rr[i]/(n+1),rr[j]/(n+1))
		}
	}
	Dn=Dn/n^2

	En=0
	for(i in 1:length(rr))
	{
		for(j in 1:length(ss))
		{
			En=En+rr[i]*rr[j]/(n+1)^2*max(ss[i]/(n+1),ss[j]/(n+1))
		}
	}
	En=En/n^2

	sigman=144*(-9*An^2+Bn+2*Cn+2*Dn+2*En)

	ICn=c()
	if(sigman>=0)
	{
		ICn=c(rhon-1.96*sqrt(sigman/n),rhon+1.96*sqrt(sigman/n))
	} else {
		ICn=rep(NA,2)
	}

	### ABC 

	rh<-runif(S, -1,1)
	omega<-rep(0,S)
	for (s in 1:S) {
		  estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
		  omega[s]<-exp(-EL(estim)$elr)
	}
	psam<-sample(rh, size=S, rep=T, prob=omega)
	sintesi<-c(quantile(psam,.025), quantile(psam,.5),  quantile(psam,.975))
	return(list(rhon=rhon,ICn=ICn,ABC=sintesi))  
}

ELCOP.vs.freq.Sch<- function(sample){

	n<-dim(sample)[1] 
	p<-dim(sample)[2]

	# Ranghi
	rr<-rank(sample[,1])
	ss<-rank(sample[,2])
	vv<-rr*ss

	# Frequentist estimate
	rhon=12 * sum(vv) /(n*(n+1)*(n-1)) - 3*(n+1)/(n-1)
	
	# Frequentist estimate of the variance
	An=sum(rr/(n+1)*ss/(n+1))/n
	Bn=sum((rr/(n+1))^2*(ss/(n+1))^2)/n

	Cn=0
	for(i in 1:length(rr))
	{
		const=rr[i]/(n+1)*ss[i]/(n+1)
		for(j in 1:length(rr))
		{
			for(k in 1:length(33))
			{
				Cn=Cn+const*(rr[k]<=rr[i])*(ss[k]<=ss[j])				
			}
		}
	}	
	Cn=Cn/n^3+0.25-An

	Dn=0
	for(i in 1:length(rr))
	{
		for(j in 1:length(ss))
		{
			Dn=Dn+ss[i]*ss[j]/(n+1)^2*max(rr[i]/(n+1),rr[j]/(n+1))
		}
	}
	Dn=Dn/n^2

	En=0
	for(i in 1:length(rr))
	{
		for(j in 1:length(ss))
		{
			En=En+rr[i]*rr[j]/(n+1)^2*max(ss[i]/(n+1),ss[j]/(n+1))
		}
	}
	En=En/n^2

	sigman=144*(-9*An^2+Bn+2*Cn+2*Dn+2*En)

	ICn=c()
	if(sigman>=0)
	{
		ICn=c(rhon-1.96*sqrt(sigman/n),rhon+1.96*sqrt(sigman/n))
	} else {
		ICn=rep(NA,2)
	}

	### ABC 

	S=1000
	rh<-runif(S, -1,1)
	omega<-rep(0,S)
	for (s in 1:S) {
		  estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
		  omega[s]<-exp(-EL.Schennach(estim)$elr)
	}
	psam<-sample(rh, size=S, rep=T, prob=omega)
	sintesi<-c(quantile(psam,.025), quantile(psam,.5),  quantile(psam,.975))
	return(list(rhon=rhon,ICn=ICn,ABC=sintesi))  
}



ELCOP.vs.freq(uu)
ELCOP.vs.freq.Sch(uu)

### Repeat simulations for Z=500 samples

require(gmm)
require(copula)

sim.cop.clay.comp<-function(theta=1.076, Z=10,ssize=100, 
	S=10^4){

cc<-claytonCopula(theta)

sintesi.finale.Sch<-matrix(rep(0, 6*Z), nrow=Z)
sintesi.finale<-matrix(rep(0, 6*Z), nrow=Z)
#
for (i in 1:Z){

  uu<-rCopula(ssize, cc)

  obj.Sch=ELCOP.vs.freq.Sch(uu)
  obj=ELCOP.vs.freq(uu)

  sintesi.finale[i,1]=obj[[2]][1]
  sintesi.finale[i,2]=obj[[1]]
  sintesi.finale[i,3]=obj[[2]][2]
  sintesi.finale[i,4:6]=obj[[3]]
  
  sintesi.finale.Sch[i,1]=obj.Sch[[2]][1]
  sintesi.finale.Sch[i,2]=obj.Sch[[1]]
  sintesi.finale.Sch[i,3]=obj.Sch[[2]][2]
  sintesi.finale.Sch[i,4:6]=obj.Sch[[3]]

  }
return(list(Schennach=sintesi.finale.Sch,Owen=sintesi.finale))
}


#  user  system elapsed 
#1395.43    1.56 1478.47 

prova.comp<-sim.cop.clay.comp(ssize=1000, Z=500, S=10^4)
prova.comp.Sch<-sim.cop.clay.comp(ssize=1000, Z=500, S=10^4)

#   user  system elapsed 
#1553.36    0.78 1677.63 

par(mfrow=c(1,2))
plot(prova.comp[,2],type="l",ylim=c(0,1),xlab="Experiment",ylab=expression(rho),
	main="Clayton",col="blue",lty=1,lwd=1)
lines(prova.comp[,1],col="blue",lty=3,lwd=0.5)
lines(prova.comp[,3],col="blue",lty=3,lwd=0.5)

lines(prova.comp[,5],col="green",lty=1,lwd=1)
lines(prova.comp[,4],col="green",lty=3,lwd=0.5)
lines(prova.comp[,6],col="green",lty=3,lwd=0.5)

abline(h=0.5,col="red",lwd=1.5)

sum(is.na(prova.comp[,1])) # 0
mean(prova.comp[,3]-prova.comp[,1]) #0.2663587
mean(prova.comp[,6]-prova.comp[,4]) #0.2596814

le.freq=prova.comp[,3]-prova.comp[,1]
le.ABC=prova.comp[,6]-prova.comp[,4]

length(le.ABC)
sum(le.freq>le.ABC)/length(le.ABC)
# 0.588

# Coverage
rho.comp=rho(claytonCopula(1.076))
sum(prova.comp[,3]>rho.comp&prova.comp[,1]<rho.comp)/500 #99.8%
sum(prova.comp[,6]>rho.comp&prova.comp[,4]<rho.comp)/500 #100%

save.image("copula-Clayton.RData")

#####################
### Multivariate rho
#####################

ELCOP.vs.freq.multi<- function(sam,var.IC,S=10^3){

	n<-dim(sam)[1] 
	d<-dim(sam)[2]

	# Ranks
	U.hat=matrix(NA,ncol=d,nrow=n)
	for(i in 1:d)
	{
		U.hat[,i]=rank(sam[,i])/n
	}
	VV1=apply(1-U.hat,1,prod)
	VV2=apply(U.hat,1,prod)

	# Frequentist estimate
	const=(d+1)/(2^d-(d+1))
	estim1=const*(2^d/n*sum(VV1)-1)
	estim2=const*(2^d/n*sum(VV2)-1)
	estim3=mean(estim1,estim2)
	
	# Frequentist estimate of the variance
	sigman=var.IC
	IC1n=c()
	IC2n=c()
	if(sigman>=0)
	{
		IC1n=c(estim1-1.96*sqrt(sigman/n),estim1+1.96*sqrt(sigman/n))
		IC2n=c(estim2-1.96*sqrt(sigman/n),estim2+1.96*sqrt(sigman/n))
		IC3n=c(estim3-1.96*sqrt(sigman/n),estim3+1.96*sqrt(sigman/n))

	} else {
		IC1n=rep(NA,2)
		IC2n=rep(NA,2)
	}

	### ABC 

	lbound=(2^d-factorial(d+1))/(factorial(d)*(2^d-(d+1)))
	rho<-runif(S, lbound,1)
	omega1<-rep(0,S)
	omega2<-rep(0,S)
	omega3<-rep(0,S)
	for (s in 1:S) {
		est3=estim3-rho[s]
		est1=const*(2^d/n*sum(VV1)-1) - rho[s]
		est2=const*(2^d/n*sum(VV2)-1) - rho[s]
		omega1[s]<-exp(-EL(est1)$elr)
		omega2[s]<-exp(-EL(est2)$elr)
		omega3[s]<-exp(-EL(est3)$elr)
	}
	psam1<-sample(rho, size=S, rep=T, prob=omega1)
	psam2<-sample(rho, size=S, rep=T, prob=omega2)
	psam3<-sample(rho, size=S, rep=T, prob=omega3)

	return(list(rhon1=estim1,rhon2=estim2,rhon3=estim3,
			IC1n=IC1n,IC2n=IC2n,IC3n=IC3n,
				ABC1=psam1,ABC2=psam2,ABC3=psam3))  
}

uu=rCopula(10^3,cc)
d=2
var.IC=((d+1)^2*(3*(4/3)^d-d-3))/(3*(1+d-2^d)^2)
try1=ELCOP.vs.freq.multi(uu,var.IC)

print("end - Multivariate rho")

# Let's repeat the same procedure 100 times

print("start - Multivariate rho - independence")

multirho.BAYvsFREQ=function(n,cc,var.IC,S)
{
	uu=rCopula(n,cc)
	try=ELCOP.vs.freq.multi(sam=uu,var.IC=var.IC,S=10^3)
	return(try)
}

# Let's try with d=2 and the variance in case of independence

final.multirho.freq <- matrix(NA,nrow=500,ncol=9)
final.multirho.bay1 <- matrix(NA,ncol=500,nrow=10^4)
final.multirho.bay2 <- matrix(NA,ncol=500,nrow=10^4)
final.multirho.bay3 <- matrix(NA,ncol=500,nrow=10^4)

for(i in 1:500)
{
	try=multirho.BAYvsFREQ(1000,cc,var.IC,S=10^4)

	final.multirho.freq[i,1]=try[[1]]
	final.multirho.freq[i,2]=try[[2]]
	final.multirho.freq[i,3]=try[[3]]

	final.multirho.freq[i,4:5]=try$IC1n
	final.multirho.freq[i,6:7]=try$IC2n
	final.multirho.freq[i,8:9]=try$IC3n

	final.multirho.bay1[,i]=try$ABC1
	final.multirho.bay2[,i]=try$ABC2
	final.multirho.bay3[,i]=try$ABC3
	
	print(i)

}

new.mat=cbind(final.multirho[,4],final.multirho[,1],final.multirho[,5],
		final.multirho[,6],final.multirho[,2],final.multirho[,7],
		final.multirho[,8],final.multirho[,3],final.multirho[,9],
		final.multirho[,10:18])
boxplot(new.mat,names=expression(IC1[0.025],rho[1],IC1[0.975],
		IC2[0.025],rho[2],IC2[0.975],
		IC3[0.025],rho[3],IC3[0.975],
		q1[0.025],q1[0.50],q1[0.975],
		q2[0.025],q2[0.50],q2[0.975],
		q3[0.025],q3[0.50],q3[0.975]),las=2)

ICcov.freq=matrix(NA,500,3)
ICcov.bay=matrix(NA,500,3)

for(i in 1:500)
{
		
	ICcov.freq[i,1]=new.mat[i,1]<0.5 & new.mat[i,3]>0.5
	ICcov.freq[i,2]=new.mat[i,4]<0.5 & new.mat[i,6]>0.5
	ICcov.freq[i,3]=new.mat[i,7]<0.5 & new.mat[i,9]>0.5
	
	ICcov.bay[i,1]=new.mat[i,10]<0.5 & new.mat[i,12]>0.5
	ICcov.bay[i,2]=new.mat[i,13]<0.5 & new.mat[i,15]>0.5
	ICcov.bay[i,3]=new.mat[i,16]<0.5 & new.mat[i,18]>0.5

}

apply(ICcov.freq,2,sum)/500
# 0.97, 0.97, 0.97

apply(ICcov.bay,2,sum)/500
# 1, 1, 1

save.image("copula-Clayton.RData")

# By supposing independence the coverage is the one expected

# Let's try with d=6

library(copula)
cc.multi=claytonCopula(1.076,dim=6)

final.multirho6.freq <- matrix(NA,nrow=500,ncol=9)
final.multirho6.bay1 <- matrix(NA,ncol=500,nrow=10^4)
final.multirho6.bay2 <- matrix(NA,ncol=500,nrow=10^4)
final.multirho6.bay3 <- matrix(NA,ncol=500,nrow=10^4)

for(i in 1:500)
{
	try6=multirho.BAYvsFREQ(1000,cc.multi,var.IC,S=10^4)

	final.multirho6.freq[i,1]=try6[[1]]
	final.multirho6.freq[i,2]=try6[[2]]
	final.multirho6.freq[i,3]=try6[[3]]

	final.multirho6.freq[i,4:5]=try6$IC1n
	final.multirho6.freq[i,6:7]=try6$IC2n
	final.multirho6.freq[i,8:9]=try6$IC3n

	final.multirho6.bay1[,i]=try6$ABC1
	final.multirho6.bay2[,i]=try6$ABC2
	final.multirho6.bay3[,i]=try6$ABC3
	
	print(i)

}

save.image("copula-Clayton.RData")

new.mat6=cbind(final.multirho6[,4],final.multirho6[,1],final.multirho6[,5],
		final.multirho6[,6],final.multirho6[,2],final.multirho6[,7],
		final.multirho6[,8],final.multirho6[,3],final.multirho6[,9],
		final.multirho6[,10:18])
boxplot(new.mat6,names=expression(IC1[0.025],rho[1],IC1[0.975],
		IC2[0.025],rho[2],IC2[0.975],
		IC3[0.025],rho[3],IC3[0.975],
		q1[0.025],q1[0.50],q1[0.975],
		q2[0.025],q2[0.50],q2[0.975],
		q3[0.025],q3[0.50],q3[0.975]),las=2)

ICcov6.freq=matrix(NA,500,3)
ICcov6.bay=matrix(NA,500,3)

for(i in 1:500)
{
		
	ICcov6.freq[i,1]=new.mat6[i,1]<0.5 & new.mat6[i,3]>0.5
	ICcov6.freq[i,2]=new.mat6[i,4]<0.5 & new.mat6[i,6]>0.5
	ICcov6.freq[i,3]=new.mat6[i,7]<0.5 & new.mat6[i,9]>0.5
	
	ICcov6.bay[i,1]=new.mat6[i,10]<0.5 & new.mat6[i,12]>0.5
	ICcov6.bay[i,2]=new.mat6[i,13]<0.5 & new.mat6[i,15]>0.5
	ICcov6.bay[i,3]=new.mat6[i,16]<0.5 & new.mat6[i,18]>0.5

}

apply(ICcov6.freq,2,sum)/500
# 0.988, 0.000, 0.988

apply(ICcov6.bay,2,sum)/500
# 1, 1, 1

## Not all the intervals contains the 0.5 value for the estimates, 
## the intervals for rho2 never contain the value 0.5,
## the covarage of the intervals for rho2 is the one expected
## the covarage of the intervals for rho3 is the one expected

save.image("copula-Clayton.RData")

#######################################
### Bootstrap estimate of the variance
#######################################

require(gmm)
require(copula)

bootvar.multirho<- function(sam,B=10000){

	n<-dim(sam)[1] 
	d<-dim(sam)[2]
	const=(d+1)/(2^d-(d+1))

	estim1=c()
	estim2=c()
	estim3=c()

	for(j in 1:B)
	{
		U.hat=matrix(NA,ncol=d,nrow=n)
		n.sam=sample(1:n,size=n,rep=T)
		sam.boot=sam[n.sam,]
		
		for(i in 1:d)
		{
			U.hat[,i]=rank(sam.boot[,i])/n
		}
		VV1=apply(1-U.hat,1,prod)
		VV2=apply(U.hat,1,prod)

		# Frequentist estimate
		estim1[j]=const*(2^d/n*sum(VV1)-1)
		estim2[j]=const*(2^d/n*sum(VV2)-1)
		estim3[j]=mean(estim1[j],estim2[j])
	
	}
	
	return(cbind(estim1,estim2,estim3))  

}

ssize=1000
cc.multi=claytonCopula(1.076,dim=6)
uu=rCopula(ssize,cc.multi)
try6.boot=bootvar.multirho(uu)

var6.boot=apply(try6.boot,2,var)
var6.boot

ssize=100
cc=claytonCopula(1.076,dim=2)
uu=rCopula(ssize,cc)
try2.boot=bootvar.multirho(uu)

var2.boot=apply(try2.boot,2,var)
var2.boot

# Monte Carlo estimate of the variance

MCvar.multirho<- function(cc,n,B=10000){

	d<-dim(cc)
	const=(d+1)/(2^d-(d+1))

	estim1=c()
	estim2=c()
	estim3=c()

	for(j in 1:B)
	{
		U.hat=matrix(NA,ncol=d,nrow=n)
		uu=rCopula(n,cc)
				
		for(i in 1:d)
		{
			U.hat[,i]=rank(uu[,i])/n
		}
		VV1=apply(1-U.hat,1,prod)
		VV2=apply(U.hat,1,prod)

		# Frequentist estimate
		estim1[j]=const*(2^d/n*sum(VV1)-1)
		estim2[j]=const*(2^d/n*sum(VV2)-1)
	
	}
	
	return(cbind(estim1,estim2))  

}

ssize=1000
cc.multi=claytonCopula(1.076,dim=6)

try6.MC=MCvar.multirho(cc.multi,n=ssize)

var6.MC=apply(try6.MC,2,var)
var6.MC


# Let's repeat the same procedure 500 times

finalboot.multirho=matrix(NA,nrow=500,ncol=18)
for(i in 1:nrow(finalboot.multirho))
{
	try.boot=multirho.BAYvsFREQ(1000,cc,max(var2.boot),S=10^4)

	finalboot.multirho[i,1]=try.boot[[1]]
	finalboot.multirho[i,2]=try.boot[[2]]
	finalboot.multirho[i,3]=try.boot[[3]]

	finalboot.multirho[i,4:5]=try.boot$IC1n
	finalboot.multirho[i,6:7]=try.boot$IC2n
	finalboot.multirho[i,8:9]=try.boot$IC3n

	finalboot.multirho[i,10:12]=try.boot$ABC1
	finalboot.multirho[i,13:15]=try.boot$ABC2
	finalboot.multirho[i,16:18]=try.boot$ABC3

}


new.matboot=cbind(finalboot.multirho[,4],finalboot.multirho[,1],finalboot.multirho[,5],
		finalboot.multirho[,6],finalboot.multirho[,2],finalboot.multirho[,7],
		finalboot.multirho[,8],finalboot.multirho[,3],finalboot.multirho[,9],
		finalboot.multirho[,10:18])
boxplot(new.matboot,names=expression(IC1[0.025],rho[1],IC1[0.975],
		IC2[0.025],rho[2],IC2[0.975],
		IC3[0.025],rho[3],IC3[0.975],
		q1[0.025],q1[0.50],q1[0.975],
		q2[0.025],q2[0.50],q2[0.975],
		q3[0.025],q3[0.50],q3[0.975]),las=2)

ICcovboot.freq=matrix(NA,500,3)
ICcovboot.bay=matrix(NA,500,3)

for(i in 1:500)
{
		
	ICcovboot.freq[i,1]=new.matboot[i,1]<0.5 & new.matboot[i,3]>0.5
	ICcovboot.freq[i,2]=new.matboot[i,4]<0.5 & new.matboot[i,6]>0.5
	ICcovboot.freq[i,3]=new.matboot[i,7]<0.5 & new.matboot[i,9]>0.5
	
	ICcovboot.bay[i,1]=new.matboot[i,10]<0.5 & new.matboot[i,12]>0.5
	ICcovboot.bay[i,2]=new.matboot[i,13]<0.5 & new.matboot[i,15]>0.5
	ICcovboot.bay[i,3]=new.matboot[i,16]<0.5 & new.matboot[i,18]>0.5

}

apply(ICcovboot.freq,2,sum)/500

# 0.202, 0.208, 0.202
# rho1 rho2 rho3
# Coverage under 0.95 for all the statistics

apply(ICcovboot.bay,2,sum)/500

# 1.00 1.00 1.00
# rho1 rho2 rho3

save.image("copula-Clayton.RData")

#### Let's try with dimension 6

library(copula)
cc.multi=claytonCopula(1.076,dim=6)

finalboot6.multirho.freq <- matrix(NA,nrow=500,ncol=9)
finalboot6.multirho.bay1 <- matrix(NA,ncol=500,nrow=10^4)
finalboot6.multirho.bay2 <- matrix(NA,ncol=500,nrow=10^4)
finalboot6.multirho.bay3 <- matrix(NA,ncol=500,nrow=10^4)

for(i in 1:500)
{
	try6.boot=multirho.BAYvsFREQ(1000,cc.multi,max(var6.boot),S=10^4)

	finalboot6.multirho.freq[i,1]=try6.boot[[1]]
	finalboot6.multirho.freq[i,2]=try6.boot[[2]]
	finalboot6.multirho.freq[i,3]=try6.boot[[3]]

	finalboot6.multirho.freq[i,4:5]=try6.boot$IC1n
	finalboot6.multirho.freq[i,6:7]=try6.boot$IC2n
	finalboot6.multirho.freq[i,8:9]=try6.boot$IC3n

	finalboot6.multirho.bay1[,i]=try6.boot$ABC1
	finalboot6.multirho.bay2[,i]=try6.boot$ABC2
	finalboot6.multirho.bay3[,i]=try6.boot$ABC3
	
	print(i)

}

save.image("copula-Clayton.RData")

new.matboot6=cbind(finalboot6.multirho[,4],finalboot6.multirho[,1],finalboot6.multirho[,5],
		finalboot6.multirho[,6],finalboot6.multirho[,2],finalboot6.multirho[,7],
		finalboot6.multirho[,8],finalboot6.multirho[,3],finalboot6.multirho[,9],
		finalboot6.multirho[,10:18])

### What is the true rho which we are simulating from?

install.packages("cubature")
library(cubature)

# Rho1
f.clayton=function(x,theta)
{
	d=length(x)
	temp=(sum(x^(0-theta))-d+1)^(-1/theta)
	return(temp)
}

int.clayton=adaptIntegrate(f.clayton,theta=1.076,lowerLimit=rep(0,6),upperLimit=rep(1,6))

d=6
rho.clayton1=((d+1)/(2^d-d-1)) * (2^d*int.clayton$integral-1)

rho.clayton1
# 0.5144599

# Rho2
f.rho2=function(x,cop)
{
	d=length(x)
	ind.cop=sum(log(x))
	cop.dens=dCopula(x,cop,log=TRUE)
	temp=ind.cop+cop.dens
	return(exp(temp))
}

int2.clayton=adaptIntegrate(f2.clayton,cop=cc.multi,
				lowerLimit=rep(0,6),upperLimit=rep(1,6))

d=6
rho.clayton2=((d+1)/(2^d-d-1)) * (2^d*int2.clayton$integral-1)

rho.clayton2
# 0.346185

ICcovboot6.freq=matrix(NA,500,2)
ICcovboot6.bay=matrix(NA,500,2)

for(i in 1:500)
{
		
	ICcovboot6.freq[i,1]=new.matboot6[i,1]<rho.clayton1 & new.matboot6[i,3]>rho.clayton1
	ICcovboot6.freq[i,2]=new.matboot6[i,4]<rho.clayton2 & new.matboot6[i,6]>rho.clayton2
	
	ICcovboot6.bay[i,1]=new.matboot6[i,10]<rho.clayton1 & new.matboot6[i,12]>rho.clayton1
	ICcovboot6.bay[i,2]=new.matboot6[i,13]<rho.clayton2 & new.matboot6[i,15]>rho.clayton2

}

apply(ICcovboot6.freq,2,sum)/500

# 0.054 0.050
# rho1  rho2
# Coverage under 0.95 for all the statistics

apply(ICcovboot6.bay,2,sum)/500

# 1.00 1.00
# rho1 rho2

save.image("copula-Clayton.RData")

## Rho1, frequentist and Bayesian estimates
par(mfrow=c(1,2))
plot(new.matboot6[,2],type="l",ylim=c(-1,1),xlab="Experiment",ylab=expression(rho[1]),
	main="Clayton 1",col="blue")
lines(new.matboot6[,1],col="blue",lty=3)
lines(new.matboot6[,3],col="blue",lty=3)

lines(new.matboot6[,11],ylim=c(-1,1),col="green")
lines(new.matboot6[,10],col="green",lty=3)
lines(new.matboot6[,12],col="green",lty=3)
abline(h=rho.clayton1,col="red",lwd=1.5)

for(i in 1:length(new.matboot6[,2]))
{
	segments(x0=i,y0=new.matboot6[i,1],x1=i,y=new.matboot6[i,3],col="blue",lwd=2)
}
points(new.matboot6[,1],col="blue")
points(new.matboot6[,3],col="blue")
abline(h=0.5,col="red",lwd=1.5)

plot(new.matboot6[,11],pch=20,ylim=c(-1,1),xlab="Experiment",ylab=expression(rho[1]),
	main="Bayesian")
for(i in 1:length(new.matboot6[,11]))
{
	segments(x0=i,y0=new.matboot6[i,10],x1=i,y=new.matboot6[i,12],col="blue",lwd=2)
}
points(new.matboot6[,10],col="blue")
points(new.matboot6[,12],col="blue")
abline(h=rho.clayton1,col="red",lwd=1.5)


# Let's try to see if the length of the confidence intervals varies
# with the dimensionality

d.incr=c(2:10,25,50)
rho.incr.dim=list()

for(i in 1:length(d.incr))
{
	mat.rho.incrd=matrix(NA,ncol=12,nrow=50)
	cc.incrd=claytonCopula(1.076,dim=d.incr[i])
	uu=rCopula(ssize,cc.incrd)
	try.boot=bootvar.multirho(uu)
	var.boot=apply(try.boot,2,var)
	for(j in 1:50)
	{
		try6.boot=multirho.BAYvsFREQ(1000,cc.incrd,max(var.boot),S=10^4)

		mat.rho.incrd[j,3]=try6.boot[[1]]
		mat.rho.incrd[j,1:2]=try6.boot$IC1n

		mat.rho.incrd[j,6]=try6.boot[[2]]
		mat.rho.incrd[j,4:5]=try6.boot$IC2n

		mat.rho.incrd[j,7:9]=try6.boot$ABC1
		mat.rho.incrd[j,10:12]=try6.boot$ABC2
		print(paste(d.incr[i],j,sep="_"))

	}
	rho.incr.dim[[i]]=mat.rho.incrd
}

le.incrd=matrix(NA,ncol=4,nrow=length(d.incr))
for(i in 1:length(d.incr))
{
	le.incrd[i,1]=mean(rho.incr.dim[[i]][,2]-rho.incr.dim[[i]][,1])
	le.incrd[i,2]=mean(rho.incr.dim[[i]][,5]-rho.incr.dim[[i]][,4])

	le.incrd[i,3]=mean(rho.incr.dim[[i]][,9]-rho.incr.dim[[i]][,7])
	le.incrd[i,4]=mean(rho.incr.dim[[i]][,12]-rho.incr.dim[[i]][,10])	
}

le.incrd

save.image("copula-Clayton.RData")

# Let's see if the value of rho(1) is fixed by varying the dimension d
# for any value of theta

theta.seq=seq(0.1,10,0.1)
rho1clayton.mat=matrix(NA,nrow=length(theta.seq),ncol=5)

for (i in 1:length(theta.seq))
{
	for(j in 1:5)
	{
		d=d.incr[j]
		int.clayton=adaptIntegrate(f.clayton,theta=theta.seq[i],
					lowerLimit=rep(0,d),upperLimit=rep(1,d))			
		rho1clayton.mat[i,j]=((d+1)/(2^d-d-1)) * 
								(2^d*int.clayton$integral-1)
		print(paste(d,i,sep="_"))
	}	
}


save.image("copula-Clayton.RData")

###############################
### Comparison with a Bayesian 
### parametric method 
###############################

bayes.copula <- function(x,nsim,burn=0.1*nsim,cop.type)
{
	# x:	sample
	# nsim:	number of simulations
	# burn: burnin
	# cop.type:	type of Copula
	
	# Pseudo-observations
	u.obs <- pobs(x)
	
	theta.post <- c()
	lik <- c()

	# Initialization
	post.den <- -Inf
	while(post.den==(-Inf))
	{
		if(cop.type=="Clayton"){
			theta.inits <- rgamma(1,0.1,0.1)		
			cop <- claytonCopula(dim=2,param=theta.inits)	
			lprior <- log(dgamma(theta.inits,0.1,0.1))	
		} else {
			if(cop.type=="Gumbel"){
				theta.inits <- rtruncnorm(1,a=1,b=Inf,mean=1,sd=10)		
				cop <- gumbelCopula(dim=2,param=theta.inits)		
				lprior <- log(dtruncnorm(theta.inits,a=1,b=Inf,
							mean=1,sd=10))			
			} else {
				theta.inits <- rnorm(1,0,sd=10)		
				cop <- frankCopula(dim=2,param=theta.inits)		
				lprior <- log(dnorm(theta.inits,0,10))		
			}
		}
		llik <- sum(log(dCopula(u.obs,cop)))
		post.den <- llik + lprior
	}
	theta.post[1] <- theta.inits
	lik[1] <- llik
		
	for(count in 2:nsim)
	{
		# Propose a new value
		theta.old <- theta.post[count-1]

		if(cop.type=="Clayton"){
			theta.prop <- rtruncnorm(1, a=0, b=Inf, 
						mean = theta.old, sd = 0.1)				
			cop <- claytonCopula(dim=2,param=theta.prop)	
			lprior <- log(dgamma(theta.prop,0.1,0.1))
			
			q.num <- log(dtruncnorm(theta.old,a=0,b=Inf,
							mean=theta.prop,sd=0.1))
			q.den <- log(dtruncnorm(theta.prop,a=0,b=Inf,
							mean=theta.old,sd=0.1))
		
		} else {
			if(cop.type=="Gumbel"){
				theta.prop <- rtruncnorm(1, a=1, b=Inf, 
							mean = theta.old, sd = 0.1)				
				cop <- gumbelCopula(dim=2,param=theta.prop)		
				lprior <- log(dtruncnorm(theta.prop,a=1,b=Inf,
							mean=1,sd=10))			

				q.num <- log(dtruncnorm(theta.old,a=1,b=Inf,
							mean=theta.prop,sd=0.1))
				q.den <- log(dtruncnorm(theta.prop,a=1,b=Inf,
							mean=theta.old,sd=0.1))

			} else {
				theta.prop <- rnorm(1,theta.old,sd=0.1)		
				cop <- frankCopula(dim=2,param=theta.prop)	
				lprior <- log(dnorm(theta.inits,0,10))	
				
				q.num <- log(dnorm(theta.old,mean=theta.prop,sd=0.1))
				q.den <- log(dnorm(theta.prop,mean=theta.old,sd=0.1))
				
			}
		}

		llik <- sum(log(dCopula(u.obs,cop)))
	
		post.num <- llik + lprior
		
		den <- post.den + q.den
		num <- post.num + q.num
					
		alpha <- runif(1)
		if(exp(num-den)>alpha)
		{
			theta.post[count] <- theta.prop
			post.den <- post.num
			lik[count] <- llik
		} else {
			theta.post[count] <- theta.post[count-1]
			lik[count] <- lik[count-1]
		}
		
	}
	
	return(theta.post)
	
}

###############
## Dimension 2
###############

# Simulations from a Clayton

bayes.np.mat <- matrix(NA,ncol=3,nrow=500)
bayes.param.mat.clay <- matrix(NA,ncol=3,nrow=500)
bayes.param.mat.frank <- matrix(NA,ncol=3,nrow=500)
bayes.param.mat.gumbel <- matrix(NA,ncol=3,nrow=500)

for(count in 1:500)
{
	uu.clayton2=rCopula(ssize,cc)

	bayes.np <- ELCOP.vs.freq(sam=uu.clayton2)
	print(paste(count,"np",sep="-"))
	bayes.param.clay <- bayes.copula(uu.clayton2,nsim=10^4,
					cop.type="Clayton")
	print(paste(count,"clayton",sep="-"))
	bayes.param.frank <- bayes.copula(uu.clayton2,nsim=10^4,
					cop.type="Frank")
	print(paste(count,"frank",sep="-"))
	bayes.param.gumbel <- bayes.copula(uu.clayton2,nsim=10^4,
					cop.type="Gumbel")
	print(paste(count,"gumbel",sep="-"))
	
	bayes.np.mat[count,] <- bayes.np$ABC
	bayes.param.mat.clay[count,] <- quantile(bayes.param.clay[1000:10^4],
										probs=c(0.025,0.5,0.975))
	bayes.param.mat.frank[count,] <- quantile(bayes.param.frank[1000:10^4],
										probs=c(0.025,0.5,0.975))
	bayes.param.mat.gumbel[count,] <- 
								quantile(bayes.param.gumbel[1000:10^4],
										probs=c(0.025,0.5,0.975))
										
	save.image("copula-Clayton.RData")
	
}

rho(cc)

plot(bayes.np.mat[,2],pch=19,cex=0.5,ylim=c(0,1),
	ylab=expression(rho),xlab="Experiment",
	main="Simulation from Clayton copula")
for(i in 1:500)
{
	lines(x=c(i,i),y=c(bayes.np.mat[i,1],bayes.np.mat[i,3]),lwd=0.3)
}
clay.bayes.clay <- matrix(NA,nrow=nrow(bayes.param.mat.clay),
					ncol=ncol(bayes.param.mat.clay))
frank.bayes.clay <- matrix(NA,nrow=nrow(bayes.param.mat.frank),
					ncol=ncol(bayes.param.mat.frank))
gumbel.bayes.clay <- matrix(NA,nrow=nrow(bayes.param.mat.gumbel),
					ncol=ncol(bayes.param.mat.gumbel))
for(i in 1:nrow(bayes.param.mat.clay))
{

	cc.temp <- claytonCopula(param=bayes.param.mat.clay[i,1])
	clay.bayes.clay[i,1] <- rho(cc.temp)
	cc.temp <- claytonCopula(param=bayes.param.mat.clay[i,2])
	clay.bayes.clay[i,2] <- rho(cc.temp)
	cc.temp <- claytonCopula(param=bayes.param.mat.clay[i,3])
	clay.bayes.clay[i,3] <- rho(cc.temp)

	cf.temp <- frankCopula(param=bayes.param.mat.frank[i,1])
	frank.bayes.clay[i,1] <- rho(cf.temp)
	cf.temp <- frankCopula(param=bayes.param.mat.frank[i,2])
	frank.bayes.clay[i,2] <- rho(cf.temp)
	cf.temp <- frankCopula(param=bayes.param.mat.frank[i,3])
	frank.bayes.clay[i,3] <- rho(cf.temp)

	cg.temp <- gumbelCopula(param=bayes.param.mat.gumbel[i,1])
	gumbel.bayes.clay[i,1] <- rho(cg.temp)
	cg.temp <- gumbelCopula(param=bayes.param.mat.gumbel[i,2])
	gumbel.bayes.clay[i,2] <- rho(cg.temp)
	cg.temp <- gumbelCopula(param=bayes.param.mat.gumbel[i,3])
	gumbel.bayes.clay[i,3] <- rho(cg.temp)

}

points(clay.bayes.clay[,2],pch=19,cex=0.5,col="orange")
for(i in 1:500)
{
	lines(x=c(i,i),y=c(clay.bayes.clay[i,1],clay.bayes.clay[i,3]),
			lwd=0.5,col="orange")
}

points(frank.bayes.clay[,2],pch=19,cex=0.5,col="blue")
for(i in 1:500)
{
	lines(x=c(i,i),y=c(frank.bayes.clay[i,1],frank.bayes.clay[i,3]),
			lwd=0.3,col="blue")
}

points(gumbel.bayes.clay[,2],pch=19,cex=0.5,col="green")
for(i in 1:500)
{
	lines(x=c(i,i),y=c(gumbel.bayes.clay[i,1],gumbel.bayes.clay[i,3]),
			lwd=0.5,col="green")
}

abline(h=rho(cc),col="red",lwd=2)

########################################
########################################
###### Tail dependence
########################################
########################################

### bivariate variable with unknown distribution,
###		use of nonparametric estimation of lambda_U.

require(gmm)
require(copula)

t=sqrt(n)
u.mat=matrix(rep((1-t/n),2),ncol=2)
C.n(u=u.mat,U=uu.clayton)

2-((1-C.n(u=u.mat,U=uu))/(t/n))
2-2*exp(1/n*sum(log(sqrt(log(1/uu[,1])*log(1/uu[,2]))/log(1/(max(uu[,1],uu[,2]))^2))))

ELCOP.vs.freq.NPTDC<- function(X,varU.freq,varL.freq,S)
{

	n=dim(X)[1] 
	t=sqrt(n)
	p=dim(X)[2]

	### Frequentist estimate
	lambdaU.mat <- matrix(NA,ncol=p,nrow=n)
	lambdaL.mat <- matrix(NA,ncol=p,nrow=n)
	for(j in 1:p)
	{
		lambdaU.mat[,j] <- rank(X[,j])>n*(1-t/n)
		lambdaL.mat[,j] <- rank(X[,j])<t
	}
	lambda_U <- (sum(apply(lambdaU.mat,1,sum)==p)/n) / (t/n)
	lambda_L <- (sum(apply(lambdaL.mat,1,sum)==p)/n) / (t/n)
	# uU.mat=matrix(rep((1-t/n),p),ncol=p)
	# uL.mat=matrix(rep(t/n,p),ncol=p)
	# lambda_U=2-((1-C.n(u=uU.mat,X=X))/(t/n))
	# lambda_L=C.n(u=uL.mat,X=X)/(t/n)

	# Frequentist intervals
	lambdaU_IC=lambda_U+c(-1.96,1.96)*sqrt(varU.freq/n)
	lambdaL_IC=lambda_L+c(-1.96,1.96)*sqrt(varL.freq/n)

	### ABC 

	lambdaUh=runif(S, 0,1)
	lambdaLh=runif(S, 0,1)
	omegaU=rep(0,S)
	omegaL=rep(0,S)
	for (s in 1:S) {
		  estimU= lambda_U - lambdaUh[s]
		  estimL= lambda_L - lambdaLh[s]
		  omegaU[s]=exp(-EL(estimU)$elr)
		  omegaL[s]=exp(-EL(estimL)$elr)
	}
	psamU=sample(lambdaUh, size=S, rep=T, prob=omegaU)
	psamL=sample(lambdaLh, size=S, rep=T, prob=omegaL)
	outU=c(quantile(psamU,.025), quantile(psamU,.5),  quantile(psamU,.975))
	outL=c(quantile(psamL,.025), quantile(psamL,.5),  quantile(psamL,.975))
	return(list(lambda.freq=c(lambda_U,lambda_L),postU=psamU,postL=psamL))  
}

bootvar.TDC<- function(sam,B=10000){

	n<-dim(sam)[1] 
	t=sqrt(n)
	d<-dim(sam)[2]

	U.mat=matrix(rep((1-t/n),d),ncol=d)
	L.mat=matrix(rep((t/n),d),ncol=d)
	lambda_U=c()
	lambda_L=c()

	for(j in 1:B)
	{
		# upperTDC estimator
		idx=sample(1:n,size=n,rep=T)
		sam.boot=sam[idx,]
		lambda_U[j]=2-((1-C.n(u=U.mat,U=sam.boot))/(t/n))
		lambda_L[j]=C.n(u=L.mat,U=sam.boot)/(t/n)

	}
	
	return(cbind(lambda_U,lambda_L))  

}

save.image("copula-Clayton.RData")

### Clayton copula

bootvar.TDC(uu.clayton)
var.TDC.clayton=apply(bootvar.TDC(uu.clayton),2,var)
ELCOP.vs.freq.NPTDC(uu.clayton,varU.freq=var.TDC.clayton[1],
				varL.freq=var.TDC.clayton[2])

Z=500
outNPTDC2.clayton.freq=matrix(NA,ncol=2,nrow=Z)
outNPTDC2.clayton.bayU <- matrix(NA,ncol=Z,nrow=10^4)
outNPTDC2.clayton.bayL <- matrix(NA,ncol=Z,nrow=10^4)

for (i in 1:Z)
{
	uu=rCopula(ssize, cc)
	obj.clayton=ELCOP.vs.freq.NPTDC(uu,
							var.TDC.clayton[1],
  							var.TDC.clayton[2],S=10^4)
  	outNPTDC2.clayton.freq[i,]=obj.clayton[[1]]
	outNPTDC2.clayton.bayU[,i] <- obj.clayton[[2]]
	outNPTDC2.clayton.bayL[,i] <- obj.clayton[[3]]
	
	print(i)
}

save.image("copula-Clayton.RData")

# Multivariate TDC
uu.cc.multi <- rCopula(ssize,cc.multi)

Z=500
outNPTDC2.clayton.multi.freq=matrix(NA,ncol=2,nrow=Z)
outNPTDC2.clayton.multi.bayU <- matrix(NA,ncol=Z,nrow=10^4)
outNPTDC2.clayton.multi.bayL <- matrix(NA,ncol=Z,nrow=10^4)

for (i in 1:Z)
{
	uu.multi=rCopula(ssize, cc.multi)
	obj.clayton=ELCOP.vs.freq.NPTDC(uu.multi,
							var.TDC.clayton[1],
  							var.TDC.clayton[2],S=10^4)
  	outNPTDC2.clayton.multi[i,]=obj.clayton[[1]]
	outNPTDC2.clayton.multi.bayU[,i] <- obj.clayton[[2]]
	outNPTDC2.clayton.multi.bayL[,i] <- obj.clayton[[3]]
	
	print(i)
}

save.image("copula-Clayton.RData")


# For the Clayton copula the upper tail dependence is 0, then to 
# look for the coverage

theta.clayton=1.076
true.lambdaU.clayton=0
true.lambdaL.clayton=2^(-1/theta.clayton)
sum(outNPTDC2.clayton[,1]<true.lambdaU.clayton&
					outNPTDC2.clayton[,3]>true.lambdaU.clayton)/Z 
#0.036 Frequentist
sum(outNPTDC2.clayton[,4]<true.lambdaL.clayton&
					outNPTDC2.clayton[,6]>true.lambdaL.clayton)/Z 
#0.00 Frequentist
sum(outNPTDC2.clayton[,7]<true.lambdaU.clayton&
					outNPTDC2.clayton[,9]>true.lambdaU.clayton)/Z 
#1 Bayesian
sum(outNPTDC2.clayton[,10]<true.lambdaL.clayton&
					outNPTDC2.clayton[,12]>true.lambdaL.clayton)/Z 
#1 Bayesian

# Coverage ok for the Bayesian estimator, too low for the frequentist one. 

# Let's see the length of the IC

# Frequentist
mean(outNPTDC2.clayton[,3]-outNPTDC2.clayton[,1])
# 0.03178564
mean(outNPTDC2.clayton[,6]-outNPTDC2.clayton[,4])
# 0.01657418

# Bayesian
mean(outNPTDC2.clayton[,9]-outNPTDC2.clayton[,7])
# 1.239105
mean(outNPTDC2.clayton[,12]-outNPTDC2.clayton[,10])
# 1.19679

# The length of Bayesian intervals is higher than the length of frequentist 
# intervals, but the coverage is not the one wanted. In this case, the Bayesian procedure seems to have a 
# better behavior. 

save.image("copula-Clayton.RData")

print("end - Tail dependence")

par(mfrow=c(2,2))
plot(outNPTDC2.clayton[,2],pch=20,
		ylim=c(-1,1),xlab="Experiment",ylab=expression(lambda[U]),
		main="Frequentist")
for(i in 1:length(outNPTDC2.clayton[,2]))
{
	segments(x0=i,y0=outNPTDC2.clayton[i,1],x1=i,
				y=outNPTDC2.clayton[i,3],col="blue",lwd=2)
}
points(outNPTDC2.clayton[,1],col="blue")
points(outNPTDC2.clayton[,3],col="blue")
abline(h=0.0,col="red",lwd=1.5)

plot(outNPTDC2.clayton[,8],pch=20,ylim=c(-1,1),
	xlab="Experiment",ylab=expression(lambda[U]),
	main="Bayesian")
for(i in 1:length(outNPTDC2.clayton[,8]))
{
	segments(x0=i,y0=outNPTDC2.clayton[i,7],x1=i,
				y=outNPTDC2.clayton[i,9],col="green",lwd=0.5)
}
points(outNPTDC2.clayton[,7],col="green")
points(outNPTDC2.clayton[,9],col="green")
abline(h=0.0,col="red",lwd=1.5)

plot(outNPTDC2.clayton[,5],pch=20,ylim=c(-1,1),xlab="Experiment",
		ylab=expression(lambda[L]),
		main="Frequentist")
for(i in 1:length(outNPTDC2.clayton[,5]))
{
	segments(x0=i,y0=outNPTDC2.clayton[i,4],x1=i,
				y=outNPTDC2.clayton[i,6],col="blue",lwd=0.5)
}
points(outNPTDC2.clayton[,4],col="blue")
points(outNPTDC2.clayton[,6],col="blue")
abline(h=2^(0-1/1.076),col="red",lwd=1.5)

plot(outNPTDC2.clayton[,11],pch=20,ylim=c(-1,1),xlab="Experiment",
		ylab=expression(lambda[L]),
		main="Bayesian")
for(i in 1:length(outNPTDC2.clayton[,11]))
{
	segments(x0=i,y0=outNPTDC2.clayton[i,10],x1=i,
				y=outNPTDC2.clayton[i,12],col="green",lwd=0.5)
}
points(outNPTDC2.clayton[,10],col="green")
points(outNPTDC2.clayton[,12],col="green")
abline(h=2^(0-1/1.076),col="red",lwd=1.5)

