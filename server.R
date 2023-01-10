server <- function(input, output) {

  
###########################################################################################################################################################################
### Results for the incomplete within-subject design with scenario 1
###########################################################################################################################################################################
  output$table1 <- renderTable({
  
  ### budgetary constraint    
  C.s=input$C.s
  C.0=input$C.0
  C.1=input$C.1
  C.2=input$C.2
  C.s01=C.s+C.0+C.1
  C.s02=C.s+C.0+C.2
  B=input$B
  
  ### variances
  var.0=input$var.0
  var.1=input$var.1
  var.2=input$var.2
  
  ### covariances
  covar.01=input$covar.01
  covar.02=input$covar.02
  covar.12=input$covar.12
  
  ### covariance matrix
  D=matrix(c(var.0,covar.01,covar.02,covar.01,var.1,covar.12,covar.02,covar.12,var.2),nrow=3)
  
  ### optimality criterion
  optcrit=input$optcrit
  lambda.01=input$lambda.01
  lambda.02=1-lambda.01

  ### design matrices and covariance matrix of responses for complete within-subject design
  X.012=diag(3)
  Z.012=diag(3)
  V.012=Z.012%*%D%*%t(Z.012)
  
  ### design matrices and covariance matrix of responses for incomplete within-subject design: treatment combination placebo and treatment 1
  X.01=X.012[1:2,]
  Z.01=Z.012[1:2,]
  V.01=Z.01%*%D%*%t(Z.01)
  
  ### design matrices and covariance matrix of responses for incomplete within-subject design: treatment combination placebo and treatment 2
  X.02=X.012[c(1,3),]
  Z.02=Z.012[c(1,3),]
  V.02=Z.02%*%D%*%t(Z.02)
  
  ### objective function for compound optimal design
  fun.compound=function(x)
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
  
  ### objective function for D_A optimal design
  fun.DA=function(x)
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
  
  ### evaluate all combinations of integer N.01 and N.02 > 0 for which budget constraint holds
  N.01.max=ceiling(B/C.s01) # maximum number of subjects who receive placebo with treatment 1
  N.01.vector=seq(0,N.01.max)
  
  results.mat=matrix(NA,nrow=(N.01.max),ncol=5)
  for(ii in 2:(N.01.max-1))
  {
    N.01=N.01.vector[ii]
    N.02=floor((B-C.s01*N.01)/C.s02)
    costs=C.s01*N.01+C.s02*N.02
    
    ### compound optimal 
    if(optcrit==1)
      objective=fun.compound(c(N.01,N.02))
  
    ### D_A optimal design
    if(optcrit==2)
      objective=fun.DA(c(N.01,N.02))
    
    results.mat[ii,]=c(N.01,N.02,N.01+N.02,costs,objective)
  }
  
  # optimal design
  optdesign=which.min(results.mat[,5])
  N01=results.mat[optdesign,1]
  N02=results.mat[optdesign,2]
  N.total=results.mat[optdesign,3]
  costs=results.mat[optdesign,4]

  ### calculate relative efficiency of uniform design
  # sample size in both treatment combinations for uniform design
  N.u=floor(B/(C.s01+C.s02))

  # RE for compound optimal design 
  if(optcrit==1)
    RE.uniform=results.mat[optdesign,5]/fun.compound(c(N.u,N.u))

  # RE for D_A optimal design
  if(optcrit==2)
    RE.uniform=(results.mat[optdesign,5]/fun.DA(c(N.u,N.u)))^0.5

  ### calculate relative efficiency of complete within-subject design
  # calculate number of subjects in complete design
  N.012=floor(B/(C.s+C.0+C.1+C.2))

  # calculate covariance matrix of means 
  covmat=solve(N.012*t(X.012)%*%solve(V.012)%*%X.012)
  
  if(optcrit==1)
  {
  # calculate variances of both contrasts 
  a.01=c(1,-1,0)
  var.contrast.01=t(a.01)%*%covmat%*%a.01
  a.02=c(1,0,-1)
  var.contrast.02=t(a.02)%*%covmat%*%a.02
  
  # objective function
  criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
  
  # relative efficiency
  RE.complete=criterion/results.mat[optdesign,5]
  }
  
  if(optcrit==2)
  {
  # matrix with coefficients of both contrasts
  a.01=c(1,-1,0)
  a.02=c(1,0,-1)
  A=t(rbind(a.01,a.02))
  
  # objective function
  criterion=det(t(A)%*%covmat%*%A)
  
  # relative efficiency
  RE.complete=(criterion/results.mat[optdesign,5])^0.5
  }
  
  # table output to be printed to screen
  values=c(N01,N02,N.total,costs,round(RE.uniform,3),round(RE.complete,3))
  description=c("Number of subjects receiving placebo and treatment 1 (N01)","Number of subjects receiving placebo and treatment 2 (N02)","Total number of subjects (N01+N02)","Costs","Efficiency of uniform  versus optimal allocation","Efficiency of incomplete versus complete within-subject design")
  output=as.data.frame(cbind(description,values))
  names(output)=NULL
  output
  })
  

  ###########################################################################################################################################################################
  ### Results for the incomplete within-subject design with scenario 2
  ###########################################################################################################################################################################
  output$table2 <- renderTable({
    
    ### budgetary constraint    
    C.s=input$C.s
    C.0=input$C.0
    C.1=input$C.1
    C.2=input$C.2
    C.s01=C.s+C.0+C.1
    C.s02=C.s+C.0+C.2
    C.s12=C.s+C.1+C.2
    B=input$B
    
    ### variances
    var.0=input$var.0
    var.1=input$var.1
    var.2=input$var.2
    
    ### covariances
    covar.01=input$covar.01
    covar.02=input$covar.02
    covar.12=input$covar.12
    
    ### covariance matrix
    D=matrix(c(var.0,covar.01,covar.02,covar.01,var.1,covar.12,covar.02,covar.12,var.2),nrow=3)
    
    ### optimality criterion
    optcrit=input$optcrit
    lambda.01=input$lambda.01
    lambda.02=1-lambda.01
    
    ### design matrices and covariance matrix of responses for complete design
    X.012=diag(3)
    Z.012=diag(3)
    V.012=Z.012%*%D%*%t(Z.012)
    
    ### design matrices and covariance matrix of responses for incomplete design: treatment combination placebo and treatment 1
    X.01=X.012[1:2,]
    Z.01=Z.012[1:2,]
    V.01=Z.01%*%D%*%t(Z.01)
    
    ### design matrices and covariance matrix of responses for incomplete design: treatment combination placebo and treatment 2
    X.02=X.012[c(1,3),]
    Z.02=Z.012[c(1,3),]
    V.02=Z.02%*%D%*%t(Z.02)
    
    ### design matrices and covariance matrix of resonses for incomplete design: treatment combination treatment 1 and treatment 2
    X.12=X.012[2:3,]
    Z.12=Z.012[2:3,]
    V.12=Z.12%*%D%*%t(Z.12)
    
    ### objective function for compound optimal design
    fun.compound=function(x)
    {
      # calculate covariance matrix of means 
      covmat=solve(x[1]* t(X.01)%*%solve(V.01)%*%X.01 + x[2]* t(X.02)%*%solve(V.02)%*%X.02 + x[3]* t(X.12)%*%solve(V.12)%*%X.12)
      
      # calculate variances of both contrasts 
      a.01=c(1,-1,0)
      var.contrast.01=t(a.01)%*%covmat%*%a.01
      a.02=c(1,0,-1)
      var.contrast.02=t(a.02)%*%covmat%*%a.02
      
      # objective function
      criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
      criterion
    }
    
    ### objective function for D_A optimal design
    fun.DA=function(x)
    {
      covmat=solve(x[1]* t(X.01)%*%solve(V.01)%*%X.01 + x[2]* t(X.02)%*%solve(V.02)%*%X.02 + x[3]* t(X.12)%*%solve(V.12)%*%X.12)
      
      # matrix with coefficients of both contrasts
      a.01=c(1,-1,0)
      a.02=c(1,0,-1)
      A=t(rbind(a.01,a.02))
      
      # objective function
      criterion=det(t(A)%*%covmat%*%A)
      criterion
    }

    ### evaluate all combinations of integer N.01 and N.02 and N.12 > 0 for which budget constraint holds
    results.mat=matrix(NA,nrow=1,ncol=6)
    N.01.max=floor(B/C.s01)
    N.01.vector=seq(1,N.01.max)
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {    
    for (ii in 1:(N.01.max-1))
    {
      N.01=N.01.vector[ii]
      B.12=B-N.01*C.s01
      N.02.max=floor(B.12/C.s02)
      N.02.vector=seq(1,N.02.max)
      for(jj in 1:(N.02.max-1))
      {
        N.02=N.02.vector[jj]
        N.12=floor((B-N.01*C.s01-N.02*C.s02)/C.s12)
        costs=N.01*C.s01+N.02*C.s02+N.12*C.s12
        
        ### compound optimal 
        if(optcrit==1)
          objective=fun.compound(c(N.01,N.02,N.12))
        
        ### D_A optimal design
        if(optcrit==2)
          objective=fun.DA(c(N.01,N.02,N.12))
        
        results.mat=rbind(results.mat,c(N.01,N.02,N.12,N.01+N.02+N.12,costs,objective))
      }
      incProgress(1/length(N.01.vector))
      Sys.sleep(0.05)
    }
                 })
    
    # optimal design
    optdesign=which.min(results.mat[,6])
    N01=results.mat[optdesign,1]
    N02=results.mat[optdesign,2]
    N12=results.mat[optdesign,3]
    N.total=results.mat[optdesign,4]
    costs=results.mat[optdesign,5]

    ### calculate relative efficiency of uniform design
    # sample size in both treatment combinations
    N.u=floor(B/(C.s01+C.s02+C.s12))

    # RE for compound optimal design 
    if(optcrit==1)
      RE.uniform=results.mat[optdesign,6]/fun.compound(c(N.u,N.u,N.u))
    
    # RE for D_A optimal design
    if(optcrit==2)
      RE.uniform=(results.mat[optdesign,6]/fun.DA(c(N.u,N.u,N.u)))^0.5
    
    ### calculate relative efficiency of complete within-subject design
    # calculate number of subjects in complete design
    N.012=floor(B/(C.s+C.0+C.1+C.2))
    
    # calculate covariance matrix of means 
    covmat=solve(N.012*t(X.012)%*%solve(V.012)%*%X.012)
    
    if(optcrit==1)
    {
    # calculate variances of both contrasts 
    a.01=c(1,-1,0)
    var.contrast.01=t(a.01)%*%covmat%*%a.01
    a.02=c(1,0,-1)
    var.contrast.02=t(a.02)%*%covmat%*%a.02
    
    # objective function
    criterion=lambda.01*var.contrast.01+lambda.02*var.contrast.02
    # relative efficiency
    RE.complete=criterion/results.mat[optdesign,6]
    }
    
    if(optcrit==2)
    {
    # matrix with coefficients of both contrasts
    a.01=c(1,-1,0)
    a.02=c(1,0,-1)
    A=t(rbind(a.01,a.02))
    
    # objective function
    criterion=det(t(A)%*%covmat%*%A)
    
    # relative efficiency
    RE.complete=(criterion/results.mat[optdesign,6])^0.5
    }
    
    # table output to be printed to screen
    values=c(N01,N02,N12,N.total,costs,round(RE.uniform,3),round(RE.complete,3))
    description=c("Number of subjects receiving placebo and treatment 1 (N01)","Number of subjects receiving placebo and treatment 2 (N02)","Number of subjects receiving treatments 1 and 2 (N12)","Total number of subjects (N01+N02+N12)","Costs","Efficiency of uniform  versus optimal allocation","Efficiency of incomplete versus complete within-subject design")
    output=as.data.frame(cbind(description,values))
    names(output)=NULL
    output
    
  })
  

}