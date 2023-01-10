#####################################################################################################################################
### optimal sample sizes per treatment condition in an incompelte within subject design with treatment combinations 01 and 02
### Author: Mirjam Moerbeek
### last update: June 29, 2022
#####################################################################################################################################

### cost constraint
C.s=500
C.0=000
C.1=423.6
C.2=635.4
C.s01=C.s+C.0+C.1
C.s02=C.s+C.0+C.2
B=100000

### covariance matrix of repeated measures within a subject
var.0=100
var.1=125
var.2=150
covar.01=0.3*sqrt(var.0*var.1)
covar.02=0.3*sqrt(var.0*var.2)
covar.12=0.3*sqrt(var.1*var.2)
D=matrix(c(var.0,covar.01,covar.02,covar.01,var.1,covar.12,covar.02,covar.12,var.2),nrow=3)

### design matrices and covariance matrix of responses for complete design
X.012=diag(3)
Z.012=diag(3)
V.012=Z.012%*%D%*%t(Z.012)

### design matrices and covariance matrix of responses for incomplete design: treatment combination placebo and treatment 1
X.01=X.012[1:2,]
Z.01=Z.012[1:2,]
V.01=Z.01%*%D%*%t(Z.01)

### design matrices and covariance matrix of resonses for incomplete design: treatment combination placebo and treatment 2
X.02=X.012[c(1,3),]
Z.02=Z.012[c(1,3),]
V.02=Z.02%*%D%*%t(Z.02)

### objective function for compound optimal design
fun1=function(x)
{
  # calculate covariance matrix of means 
  covmat=solve(x[1]* t(X.01)%*%solve(V.01)%*%X.01 + x[2]* t(X.02)%*%solve(V.02)%*%X.02)
 
  # calculate variances of both contrasts 
  a.01=c(1,-1,0)
  var.contrast.01=t(a.01)%*%covmat%*%a.01
  a.02=c(1,0,-1)
  var.contrast.02=t(a.02)%*%covmat%*%a.02
print(c(var.contrast.01,var.contrast.02))
  
  # objective function
  criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
  criterion
}

### objective function for D_A optimal design
fun2=function(x)
{
  covmat=solve(x[1]* t(X.01)%*%solve(V.01)%*%X.01 + x[2]* t(X.02)%*%solve(V.02)%*%X.02)
  
  # matrix with coefficients of both contrasts
  a.01=c(1,-1,0)
  a.02=c(1,0,-1)
  A=t(rbind(a.01,a.02))
  
  # objective function
  criterion=det(t(A)%*%covmat%*%A)
  criterion
}

#####################################################################################################################################
### evaluate all combinations of integer N.01 and N.02 >0 for which budget constraint holds
#####################################################################################################################################
N.01.max=ceiling(B/C.s01)
N.01.vector=seq(0,N.01.max)
results.mat=matrix(NA,nrow=(N.01.max),ncol=8)
#

for(ii in 2:(N.01.max-1))
{
  #print(ii)
  N.01=N.01.vector[ii]
  N.02=floor((B-C.s01*N.01)/C.s02)
  costs=C.s01*N.01+C.s02*N.02

  ### compound optimal design with lambda.01=0.25
  lambda.01=0.25
  lambda.02=1-lambda.01
  objective1a=fun1(c(N.01,N.02))

  ### compound optimal design with lambda.01=0.5
  lambda.01=0.5
  lambda.02=1-lambda.01
  objective1b=fun1(c(N.01,N.02))

  ### compound optimal design with lambda.01=0.75
  lambda.01=0.75
  lambda.02=1-lambda.01
  objective1c=fun1(c(N.01,N.02))

  ### D_A optimal design
  objective2=fun2(c(N.01,N.02))
  results.mat[ii,]=c(N.01,N.02,N.01+N.02,costs,objective1a,objective1b,objective1c,objective2)
}

# compound optimal design with lambda.01=0.25
optdes1a=which.min(results.mat[,5])
results.mat[optdes1a,1:3]
# compound optimal design with lambda.01=0.5
optdes1b=which.min(results.mat[,6])
results.mat[optdes1b,1:3]
# compound optimal design with lambda.01=0.75
optdes1c=which.min(results.mat[,7])
results.mat[optdes1c,1:3]
# D_A optimal design
optdes2=which.min(results.mat[,8])
results.mat[optdes2,1:3]

#####################################################################################################################################
### calculate relative efficiency of uniform design: compound optimal design and D_A optimal design
#####################################################################################################################################
# sample size in both treatment combinations
N.u=floor(B/(C.s01+C.s02))
costs=C.s01*N.u+C.s02*N.u

# RE for compound optimal design with lambda.01=0.25
lambda.01=0.25
lambda.02=1-lambda.01
RE.1u=results.mat[optdes1a,5]/fun1(c(N.u,N.u))
round(RE.1u,3)

# RE for compound optimal design with lambda.01=0.5
lambda.01=0.5
lambda.02=1-lambda.01
RE.1u=results.mat[optdes1b,6]/fun1(c(N.u,N.u))
round(RE.1u,3)

# RE for compound optimal design with lambda.01=0.75
lambda.01=0.75
lambda.02=1-lambda.01
RE.1u=results.mat[optdes1c,7]/fun1(c(N.u,N.u))
round(RE.1u,3)

# RE for D_A optimal design
RE.2u=(results.mat[optdes2,8]/fun2(c(N.u,N.u)))^0.5
round(RE.2u,3)

#####################################################################################################################################
### calculate relative efficiency of complete design: compound optimal design
#####################################################################################################################################
# calculate number of subjects in complete design
N.012=floor(B/(C.s+C.0+C.1+C.2))

# calculate covariance matrix of means 
covmat=solve(N.012*t(X.012)%*%solve(V.012)%*%X.012)

# calculate variances of both contrasts 
a.01=c(1,-1,0)
var.contrast.01=t(a.01)%*%covmat%*%a.01
a.02=c(1,0,-1)
var.contrast.02=t(a.02)%*%covmat%*%a.02

# calculate RE for compound optimal design with lambda.01=0.25
lambda.01=0.25
lambda.02=1-lambda.01
criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
RE.1a=criterion/results.mat[optdes1a,5]
round(RE.1a,3)

# calculate RE for compound optimal design with lambda.01=0.5
lambda.01=0.5
lambda.02=1-lambda.01
criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
RE.1b=criterion/results.mat[optdes1b,6]
round(RE.1b,3)

# calculate RE for compound optimal design with lambda.01=0.75
lambda.01=0.75
lambda.02=1-lambda.01
criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
RE.1c=criterion/results.mat[optdes1c,7]
round(RE.1c,3)

#####################################################################################################################################
### calculate relative efficiency of complete design: D_A optimal design
#####################################################################################################################################
# calculate number of subjects in complete design
N.012=floor(B/(C.s+C.0+C.1+C.2))

# calculate covariance matrix of means 
covmat=solve(N.012*t(X.012)%*%solve(V.012)%*%X.012)

# matrix with coefficients of both contrasts
a.01=c(1,-1,0)
a.02=c(1,0,-1)
A=t(rbind(a.01,a.02))
 
# objective function
criterion=det(t(A)%*%covmat%*%A)

# relative efficiency
RE.2=(criterion/results.mat[optdes2,8])^0.5
round(RE.2,3)








#####################################################################################################################################
### plots with sample sizes as a function of lambda_1
#####################################################################################################################################

lambda1.vec=seq(0,1,length=101)
N1=rep(NA,101)
N2=rep(NA,101)
N=rep(NA,101)

for(jj in 2:100)
	{
	lambda.01=lambda1.vec[jj]
	lambda.02=1-lambda.01

	N.01.max=ceiling(B/C.s01)
	N.01.vector=seq(0,N.01.max)
	results.mat=matrix(NA,nrow=(N.01.max+1),ncol=6)

	for(ii in 2:(N.01.max-1))
	{
		N.01=N.01.vector[ii]
  		N.02=floor((B-C.s01*N.01)/C.s02)
  		costs=C.s01*N.01+C.s02*N.02
  		objective1=fun1(c(N.01,N.02))
  		objective2=fun2(c(N.01,N.02))
  		results.mat[ii,]=c(N.01,N.02,N.01+N.02,costs,objective1,objective2)
	}

	# compound optimal design
	optdes1=which.min(results.mat[,5])
	results.mat[optdes1,]

	N1[jj]=results.mat[optdes1,1]
	N2[jj]=results.mat[optdes1,2]
	N[jj]=results.mat[optdes1,3]
	}

N1[1]=0
N2[1]=floor(B/(C.s+C.0+C.2))

N1[101]=floor(B/(C.s+C.0+C.1))
N2[101]=0

N[1]=N1[1]+N2[1]
N[101]=N1[101]+N2[101]

#dev.new(width=16/2.54, height=16/2.54)
plot(lambda1.vec,N1,ylim=c(0,110),pch=16,col="black",xlab=(~paste(lambda[1])),ylab="Sample size", cex.axis=1.25,cex.lab=1.25)
points(lambda1.vec,N2,pch=17,col="black")
points(lambda1.vec,N,pch=18,col="black")
text(0.11,7,expression(paste(N["01"]  ^"*" )),cex=1.25)
text(0.9,7,expression(paste(N["02"]  ^"*" )),cex=1.25)
text(0.5,92.5,expression(paste(N ^"*" )),cex=1.25)




