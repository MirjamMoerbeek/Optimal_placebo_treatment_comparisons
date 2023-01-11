#####################################################################################################################################
### optimal sample sizes per treatment condition in an incompelte within subject design with treatment combinations 01 and 02
### Author: Mirjam Moerbeek
### last update: January 11, 2023
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
criterion=function(x)
{
  # calculate covariance matrix of means 
  covmat=solve(x[1]* t(X.01)%*%solve(V.01)%*%X.01 + x[2]* t(X.02)%*%solve(V.02)%*%X.02)
 
  # calculate variances of both contrasts 
  a.01=c(1,-1,0)
  var.contrast.01=t(a.01)%*%covmat%*%a.01
  a.02=c(1,0,-1)
  var.contrast.02=t(a.02)%*%covmat%*%a.02
  
  # objective function
  criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
  criterion
}

#####################################################################################################################################
### plot with sample sizes as a function of lambda_1
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
	results.mat=matrix(NA,nrow=(N.01.max+1),ncol=5)

	for(ii in 2:(N.01.max-1))
	{
		N.01=N.01.vector[ii]
  		N.02=floor((B-C.s01*N.01)/C.s02)
  		costs=C.s01*N.01+C.s02*N.02
  		objective1=criterion(c(N.01,N.02))
  		results.mat[ii,]=c(N.01,N.02,N.01+N.02,costs,objective1)
	}

	# compound optimal design
	optdes=which.min(results.mat[,5])
	results.mat[optdes,]

	N1[jj]=results.mat[optdes,1]
	N2[jj]=results.mat[optdes,2]
	N[jj]=results.mat[optdes,3]
	}

N1[1]=0
N2[1]=floor(B/(C.s+C.0+C.2))

N1[101]=floor(B/(C.s+C.0+C.1))
N2[101]=0

N[1]=N1[1]+N2[1]
N[101]=N1[101]+N2[101]

dev.new(width=16/2.54, height=16/2.54)
plot(lambda1.vec,N1,ylim=c(0,110),pch=16,col="black",xlab=(~paste(lambda[1])),ylab="Sample size", cex.axis=1.25,cex.lab=1.25)
points(lambda1.vec,N2,pch=17,col="black")
points(lambda1.vec,N,pch=18,col="black")
text(0.11,7,expression(paste(N["01"]  ^"*" )),cex=1.25)
text(0.9,7,expression(paste(N["02"]  ^"*" )),cex=1.25)
text(0.5,92.5,expression(paste(N ^"*" )),cex=1.25)

#####################################################################################################################################
### evaluate all combinations of natural N.01 and N.02 for which budget constraint holds and calculate objective with lambda_1=0.75
#####################################################################################################################################

N.01.max=ceiling(B/C.s01)
N.01.vector=seq(0,N.01.max)
results.mat=matrix(NA,nrow=(N.01.max),ncol=5) # matrix to store results

for(ii in 2:(N.01.max-1))
{
  N.01=N.01.vector[ii]
  N.02=floor((B-C.s01*N.01)/C.s02)
  costs=C.s01*N.01+C.s02*N.02

  ### compound optimal design with lambda.01=0.75
  lambda.01=0.75
  lambda.02=1-lambda.01
  objective=criterion(c(N.01,N.02))

  ### store results in matrix	
  results.mat[ii,]=c(N.01,N.02,N.01+N.02,costs,objective)
}

#####################################################################################################################################
### derive compound optimal design with lambda_1=0.75 based on findings in previous step
#####################################################################################################################################

optdes=which.min(results.mat[,5])
results.mat[optdes,1:3]

#####################################################################################################################################
### calculate relative efficiency of uniform design versus compound optimal design with lambda_1=0.75
#####################################################################################################################################

# sample size in both treatment combinations
N.u=floor(B/(C.s01+C.s02))
costs=C.s01*N.u+C.s02*N.u

lambda.01=0.75
lambda.02=1-lambda.01
RE.1u=results.mat[optdes,5]/criterion(c(N.u,N.u))
round(RE.1u,3)

#####################################################################################################################################
### effect size as a function of power
#####################################################################################################################################

# variance of first contrast
lambda.01=1
lambda.02=1-lambda.01
variance01=criterion(c(64,36))

# variance of second contrast
lambda.01=0
lambda.02=1-lambda.01
variance02=criterion(c(64,36))

# effect size for first contrast for power=0.8 and type I error rate=0.05
ES01=(qnorm(0.95)+qnorm(0.8))*sqrt(variance01)

# effect size for second contrast for power=0.8 and type I error rate=0.05
ES02=(qnorm(0.95)+qnorm(0.8))*sqrt(variance02)

# standardized effect size for first contrast for power=0.8 and type I error rate=0.05
sd.pooled=sqrt(var.0/2+var.1/2)
ES01/sd.pooled

# standardized effect size for second contrast for power=0.8 and type I error rate=0.05
sd.pooled=sqrt(var.0/2+var.2/2)
ES02/sd.pooled

#####################################################################################################################################
### power as a function of effect size (Cohen's d)
#####################################################################################################################################

# Cohen's d
d=0.5

# variance of first contrast
lambda.01=1
lambda.02=1-lambda.01
variance01=criterion(c(64,36))

# variance of second contrast
lambda.01=0
lambda.02=1-lambda.01
variance02=criterion(c(64,36))

# unstandardized effect size for first contrast
sd.pooled=sqrt(var.0/2+var.1/2)
ES01=d*sd.pooled

# unstandardized effect size for second contrast
sd.pooled=sqrt(var.0/2+var.2/2)
ES02=d*sd.pooled


# power for first contrast
pnorm(ES01/sqrt(variance01)-qnorm(0.95))

# power for second contrast
pnorm(ES02/sqrt(variance02)-qnorm(0.95))
