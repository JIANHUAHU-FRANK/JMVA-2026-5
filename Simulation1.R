
#sink("SRPS-Example2.txt", append=TRUE,split=TRUE)

#install.packages("glmnet")
#install.packages("rrpack") 
#install.packages("MASS") 
#install.packages("Matrix") 

library(Matrix)
library(glmnet)
library(rrpack)
library(MASS)

nreplication<-2

sample.size<-c(200)         

q<-2000     # The number of response variables
q1<-10     # The number of active response variables 
p<-100      # The number of predictor variables
p1<-10      # Then number of active predictor variables

theta.h<-6.5   # normal, n=200, 8.5(2.7551),7.5()
normal<-TRUE    # t(3), n=200, 10(1.5580), 8.5(0.8427)


newdata<-c()
#########################################################################
# Notations:
#########################################################################
# Group Response best subset selector: 
# Y-nxq observation matrix, alpha: significant level,
# X-nxp predictor matrix, p<n is required,
# number.groups: the number of response variables,
# A: Gorup list,
# steplength.alpha: searching step length of alpha in (0,0.05], 0.00125 is recommended.
#########################################################################

######################################################################################
# GRBS with line search of alpha in (0,0.05], gamma=0.95+2*alpha is from 0.95 to 1.15. 
######################################################################################
GRBS<-function(Y, X, number.groups,A,steplength.alpha)
{
  n<-dim(Y)[1]
  q<-dim(Y)[2]
  p<-dim(X)[2]
  
  d<-matrix(0,1,number.groups)
  for (i in 1:number.groups)
  {d[i]<-length(A[[i]])}
  
  qq<-q-sum(d)+number.groups
  
  one<-matrix(1,n,1)
  projection1<-one%*%t(one)/n
  
  Y<-(diag(n)-projection1)%*%Y
  X<-(diag(n)-projection1)%*%X
  
  projection<-X%*%ginv(t(X)%*% X)%*% t(X)
  hat.sigma<-matrix(0,q,q)
  
  for (j in 1:q)
  {
    hat.sigma[j,j]<-(t(Y[,j]) %*% (diag(n)-projection) %*% Y[,j]/(n-p))^(1/2) 
  }
  Y<-Y %*% ginv(hat.sigma)
  
  alpha.set<-seq(0,0.1,steplength.alpha)
  bb<-floor(0.1/steplength.alpha)
  
  cv<-c(matrix(0,1,bb))
  gcv<-c(matrix(0,1,bb)) 
  delta<-matrix(0, bb, q)
  
  for (kk in 1:bb) 
  {
    alpha<-alpha.set[kk+1]
    gamma<-0.95+2*alpha
    c<-qnorm(1-alpha, mean=0,sd=1,lower.tail = TRUE, log.p = FALSE)
    
    d1<-matrix(0,1,qq)
    A1<-matrix(0,1,qq)
    B1<-matrix(0,1,qq)
    C1<-matrix(0,1,qq)
    
    for (l in 1:qq)
    {    
      for (r in A[[l]])
      {
        A1[l]<-A1[l]+t(Y[,r]) %*% (diag(n)-projection) %*% Y[,r]
        B1[l]<-B1[l]+t(Y[,r]) %*% projection %*% Y[,r]
        C1[l]<-C1[l]+(t(Y[,r]) %*% (diag(n)-projection) %*% Y[,r])^2
      }
      
      d1[l]<-p+c*sqrt(2*n*p/(n-p))*sqrt(C1[l])/A1[l]
      
      lambda<-(n-p)/d1[l]^(gamma)
      
      if (A1[l]/(B1[l]^(gamma))<=lambda)
      {
        for (r in A[[l]])
        {
          delta[kk,r]<-1
        }
      }
    }
    
    Ps<-diag(c(delta[kk,]))
    residual<-(Y-projection%*%Y%*%Ps)
    
    J1<-c()
    for (j in 1:q){
      if (delta[kk,j]==1){J1<-cbind(J1,j)}
    }
    dd<-n*q*(1-p*length(J1)/(n*q))^2
    
    CV.matrix<-matrix(0,n,q)
    for (ii in 1:n)
    {
      for (jj in 1:q)
      {
        CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-projection[ii,ii]*Ps[jj,jj])^2
      }
    }
    cv[kk]<-sum(CV.matrix)/(n*q)
    gcv[kk]<-norm((Y-projection%*%Y%*%Ps),"F")^2/dd
  }
  nnumber<-max(which.min(cv),which.min(gcv))
  
  decision.y<-delta[nnumber,]
  bestset.y<-c()
  
  for (j in 1:q){
    if (delta[nnumber,j]==1){bestset.y<-cbind(bestset.y,j)}
  }
  result<-list(bestset=bestset.y, bestdecision=decision.y)
  return(result)
}

#############################################################################
# Response best-subst selector (RBS)
###################################################################################
RBS<-function(Y,X,steplength.alpha)
{
  n<-dim(Y)[1]
  q<-dim(Y)[2]
  p<-dim(X)[2]
  
  one<-matrix(1,n,1)
  projection1<-one%*%t(one)/n
  
  Y<-(diag(n)-projection1)%*%Y
  X<-(diag(n)-projection1)%*%X
  
  projection<-X%*%ginv(t(X)%*% X)%*% t(X)
  hat.sigma<-matrix(0,q,q)
  
  for (j in 1:q)
  {
    hat.sigma[j,j]<-(t(Y[,j]) %*% (diag(n)-projection) %*% Y[,j]/(n-p))^(1/2) 
  }
  Y<-Y %*% ginv(hat.sigma)
  
  bb<-floor(0.1/steplength.alpha)
  
  cv<-matrix(0,1,bb)
  gcv<-matrix(0,1,bb)
  alpha.set<-seq(0,0.1,steplength.alpha)
  
  delta<-matrix(0, bb, q)
  
  for (kk in 1:bb)
  {    
    alpha<-alpha.set[kk+1]
    gamma<-0.95+2*alpha
    c<-qnorm(1-alpha, mean=0,sd=1,lower.tail = TRUE, log.p = FALSE)
    d<-p+c*sqrt(2*n*p/(n-p))     
    lambda<-(n-p)/d^(gamma)
    
    for (jj in 1:q) 
    {
      A<-t(Y[,jj]) %*% (diag(n)-projection) %*% Y[,jj]
      B<-(t(Y[,jj]) %*% projection %*% Y[,jj])^(gamma)
      if (A/B<=lambda) 
      {delta[kk,jj]<-1} 
    }
    
    Ps<-diag(c(delta[kk,]))
    residual<-(Y-projection%*%Y%*%Ps)
    
    J1<-c()
    for (j in 1:q)
    {
      if (delta[kk,j]==1){J1<-cbind(J1,j)}
    }
    
    dd<-n*q*(1-p*length(J1)/(n*q))^2
    
    CV.matrix<-matrix(0,n,q)
    for (ii in 1:n)
    {
      for (jj in 1:q)
      {
        CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-projection[ii,ii]*Ps[jj,jj])^2
      }
    }
    cv[kk]<-sum(CV.matrix)/(n*q)
    gcv[kk]<-norm((Y-projection%*%Y%*%Ps),"F")^2/dd
  }
  nnumber<-max(which.min(cv),which.min(gcv))
  
  decision.y<-delta[nnumber,]
  bestset.y<-c()
  
  for (j in 1:q){
    if (delta[nnumber,j]==1){bestset.y<-cbind(bestset.y,j)}
  } 
  
  result<-list(bestset=bestset.y, bestdecision=decision.y)
  return(result)
}

######################################################################################
# Predictor variable selection:
# X: nxp predictor matrix, Y: nxq observation matrix, 
# G=diag(d_1,...,d_q), d_i=0,1 indicator, d_i=1 means y_i selected,
# number.lambda: the number of lambda taking values, for example, 100,
######################################################################################
predictor.choice<-function(X,Y,G,number.lambda)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(Y)[2]
  q1<-sum(G)
  
  if (q1>1){
  projection<-X%*%ginv(t(X)%*%X) %*%t(X)
  
  hat.sigma<-matrix(0,q,q)
  Y.originalset<-c()
  
  for (j in 1:q)
  {
    if (G[j,j]==1) {Y.originalset<-cbind(Y.originalset, j)}
  }
  
  for (j in 1:q)
  {
    hat.sigma[j,j]<-(t(Y[,j]) %*% (diag(n)-projection) %*% Y[,j]/(n-p))^(1/2) 
  }
  Y<-Y %*% solve(hat.sigma) %*% G
  values<-matrix(0,1,p)
  
  for(j in 1:p)
  {
   values[j]<-norm(t(X[,j]) %*% Y,"2")/n/q1^(1/2)
  }   
  lambda.max<-2*max(values)
  
  mid.lambda<-floor(number.lambda/2)
  
  step.size<-lambda.max/number.lambda
  
  lambda.set<-seq(0,lambda.max,step.size)
  
  cv<-c(matrix(99,1,number.lambda))
  gcv<-c(matrix(99,1,number.lambda)) 
  
  x_hat<-matrix(0,number.lambda,p)
  
  Theta<-list()
  for (kk in 1:number.lambda)
  {
    Theta[[kk]]<-matrix(0,p,q)
  }

  for (kk in mid.lambda:1)
  {
 #   print(kk)
    
    lambda<-lambda.set[kk+1]
    initial.estimator<-solve(t(X) %*% X+0.001*diag(p)) %*% t(X) %*% Y
    for (j in 1:p)
    {
      R<-Y-X[,-j] %*% initial.estimator[-j,]
      omega0<-q1^(1/2)
      s<-t(X[,j]) %*% R/(n*lambda*omega0)
      
      if (norm(s,"2")>1) 
      {
        for (k in Y.originalset)
        {
          Theta[[kk]][j,k]<-(t(X[,j])%*%(Y[,k]-X%*%initial.estimator[,k])/n+
            initial.estimator[j,k])/(1+lambda*omega0/norm(initial.estimator[j,],"2"))
        }
      }
      else 
      {Theta[[kk]][j,]<-initial.estimator[j,]}
    }
    
    for (i in 1:p)
    {
      R<-Y-X[,-i]%*%Theta[[kk]][-i,]
      omega<-q1^(1/2)
      ss<-t(X[,i]) %*% R/(n*lambda*omega)
      if (norm(ss,"2")>1) 
      {x_hat[kk,i]<-1}
      else
      {Theta[[kk]][i,]<-matrix(0,1,q)}
    }
    
    I_hat<-c()
    for (i in 1:p){
      if (x_hat[kk,i]==1){I_hat<-cbind(I_hat,i)}
    }
    
 #   print(c(I_hat))
    
    if(length(I_hat)!=p)
    {
      if(length(I_hat)!=0)
      {
      Px<-X[,I_hat]%*%ginv(t(X[,I_hat])%*%X[,I_hat]) %*%t(X[,I_hat])
      residual<-(Y-Px%*%Y)
      dd<-n*q1*(1-length(I_hat)/n)^2
      
      CV.matrix<-matrix(0,n,q)
      for (ii in 1:n)
      {
        for (jj in Y.originalset)
        {
          CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-Px[ii,ii]*G[jj,jj])^2
        }
      }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(residual,"F")^2/dd
      }else{
        residual<-Y
        dd<-n*q1
        
        CV.matrix<-matrix(0,n,q)
        for (ii in 1:n)
        {
          for (jj in Y.originalset)
          {
            CV.matrix[ii,jj]<-residual[ii,jj]^2
          }
        }
        cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
        gcv[kk]<-norm(residual,"F")^2/dd
      }
    }else{
      residual<-(Y-projection%*%Y)
      dd<-n*q1*(1-p/n)^2
      
      CV.matrix<-matrix(0,n,q)
      for (ii in 1:n)
      {
        for (jj in Y.originalset)
        {
          CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-projection[ii,ii]*G[jj,jj])^2
        }
      }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(residual,"F")^2/dd
      break}
  }
  
  start.point<-mid.lambda+1 
  for (kk in start.point:number.lambda)
  {
 #   print(kk)
    
    lambda<-lambda.set[kk+1]
    initial.estimator<-solve(t(X) %*% X+0.001*diag(p)) %*% t(X) %*% Y
    for (j in 1:p)
    {
      R<-Y-X[,-j] %*% initial.estimator[-j,]
      omega0<-q1^(1/2)
      s<-t(X[,j]) %*% R/(n*lambda*omega0)
      
      if (norm(s,"2")>1) 
      {
        for (k in Y.originalset)
        {
          Theta[[kk]][j,k]<-(t(X[,j])%*%(Y[,k]-X%*%initial.estimator[,k])/n+
          initial.estimator[j,k])/(1+lambda*omega0/norm(initial.estimator[j,],"2"))
        }
      }
      else 
        {Theta[[kk]][j,]<-initial.estimator[j,]}
    }
    
    for (i in 1:p)
    {
      R<-Y-X[,-i]%*%Theta[[kk]][-i,]
      omega<-q1^(1/2)
      ss<-t(X[,i]) %*% R/(n*lambda*omega)
      if (norm(ss,"2")>1) 
      {x_hat[kk,i]<-1}
      else
        {Theta[[kk]][i,]<-matrix(0,1,q)}
    }

    I_hat<-c()
    for (i in 1:p){
    if (x_hat[kk,i]==1){I_hat<-cbind(I_hat,i)}
    }
    
 #   print(c(I_hat))
    
    if(length(I_hat)!=0)
    {
    Px<-X[,I_hat]%*%ginv(t(X[,I_hat])%*%X[,I_hat]) %*%t(X[,I_hat])
    residual<-(Y-Px%*%Y)
    dd<-n*q1*(1-length(I_hat)/n)^2
    
     CV.matrix<-matrix(0,n,q)
     for (ii in 1:n)
     {
       for (jj in Y.originalset)
       {
         CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-Px[ii,ii]*G[jj,jj])^2
       }
     }
    cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
    gcv[kk]<-norm(residual,"F")^2/dd
    }else{
      CV.matrix<-matrix(0,n,q)
      for (ii in 1:n)
      {
        for (jj in Y.originalset)
        {
          CV.matrix[ii,jj]<-Y[ii,jj]^2
        }
      }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(residual,"F")^2/(n*q1)
      break}
  }
  nnumber<-min(which.min(cv),which.min(gcv))

    I_hat<-c()
    for (i in 1:p){
      if (x_hat[nnumber,i]==1){I_hat<-cbind(I_hat,i)}
    }
    Beta<-Theta[[nnumber]] %*% hat.sigma
    
  result<-list(Beta=Beta, x.bestdecision=x_hat[nnumber,], x.bestset=I_hat)
  }else
  {
  result<-list(Beta=matrix(0,p,q), x.bestdecision=matrix(0,1,p), x.bestset=c())
  }
  return(result)
}

############################################################
# Predictor-variable selection with group patern:
# number.group: the number of groups,
# B: group list,
#number.lambda: the number of lambda raking values
#############################################################

predictor.group<-function(X,Y,G,number.group,B,number.lambda)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  q<-dim(Y)[2]
  q1<-sum(G)
  projection<-X%*%ginv(t(X)%*%X) %*%t(X)
 
  if (q1>0)
  {
  Y.originalset<-c()
  for (j in 1:q)
  {
    if (G[j,j]==1) {Y.originalset<-cbind(Y.originalset, j)}
  }
  
  hat.sigma<-matrix(0,q,q)
  for (j in 1:q)
  {
    hat.sigma[j,j]<-(t(Y[,j]) %*% (diag(n)-projection) %*% Y[,j]/(n-p))^(-1/2) 
  }
  
  Y<-Y %*% hat.sigma %*% G
  
  Theta<-list()
  for (kk in 1:number.lambda)
  {
    Theta[[kk]]<-matrix(0,p,q)
  }
  
  values<-matrix(0,1,p)
  for(j in 1:p)
  {
    values[j]<-norm(t(X[,j]) %*% Y,"2")/(n*q1^(1/2))
  }
  u1<-max(values)
  
  values.groups<-matrix(0,1,number.group)
  for (j in 1:number.group)
  {
    v<-(length(B[[j]])*q1)^(1/2)
    values.groups[j]<-norm(t(X[,B[[j]]]) %*% Y,"F")/(n*v)
  }
  u2<-max(values.groups)
  lambda.max<-2*max(u1,u2)
  
  step.size<-lambda.max/number.lambda
  lambda.set<-seq(0,lambda.max,step.size)
  
  cv<-c(matrix(99,1,number.lambda))
  gcv<-c(matrix(99,1,number.lambda)) 
  x_hat<-matrix(0,number.lambda,p)
  
  mid.lambda<-floor(number.lambda/2)

  for (kk in mid.lambda:1)
  {
    initial.estimator<-ginv(t(X) %*% X+0.001*diag(p)) %*% t(X) %*% Y
    
    lambda<-lambda.set[kk+1]
    for (j in 1:number.group)
    {
      R<-Y-X[,-B[[j]]] %*% initial.estimator[-B[[j]],]
      omega0<-(length(B[[j]])*q1)^(-1/2)
      s<-t(X[,B[[j]]]) %*% R/(n*lambda*omega0)
      
      if (norm(s,"F")>1)
      {
        for (k in Y.originalset)
        {
          for (l in B[[j]])
          {
            Theta[[kk]][l,k]<-(t(X[,l])%*%(Y[,k]-X%*%initial.estimator[,k])/n+
                 initial.estimator[l,k])/(1+lambda*omega0/norm(initial.estimator[B[[j]],],"F"))
          }
        }
      }
      else
      {
        for (l in B[[j]])
        {
          Theta[[kk]][l,]<-initial.estimator[l,]
        }
      }
    }
  
    for (i in 1:number.group)
    {
      R<-Y-X[,-B[[i]]]%*%Theta[[kk]][-B[[i]],]
      omega<-(length(B[[i]])*q1)^(1/2)
      ss<-t(X[,B[[i]]]) %*% R/(n*lambda*omega)
      
      if (norm(ss,"F")>1)
      {
        for (l in B[[i]])
        {x_hat[kk,l]<-1}
      }
      else
      {
        for (l in B[[i]])
        {
          Theta[[kk]][l,]<-matrix(0,1,q)
        }
      }
    }
    
    I_hat<-c()
    for (i in 1:p){
      if (x_hat[kk,i]==1){I_hat<-cbind(I_hat,i)}
    }
    
  #  print(c(I_hat))

    if(length(I_hat)!=p)
    {
      if (length(I_hat)!=0)
      {
      Px<-X[,I_hat]%*%ginv(t(X[,I_hat])%*%X[,I_hat]) %*%t(X[,I_hat])
      residual<-(Y-Px%*%Y)
      
      dd<-n*q1*(1-length(I_hat)/n)^2
      
      CV.matrix<-matrix(0,n,q)
      for (ii in 1:n)
      {
        for (jj in Y.originalset)
        {
          CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-Px[ii,ii]*G[jj,jj])^2
        }
      }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(residual,"F")^2/dd
      }else{
        residual<-Y
        dd<-n*q1
        CV.matrix<-matrix(0,n,q)
        for (ii in 1:n)
        {
          for (jj in Y.originalset)
          {
            CV.matrix[ii,jj]<-residual[ii,jj]^2
          }
        }
        cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
        gcv[kk]<-norm(residual,"F")^2/dd
        
      }
    }else{ 
    residual<-(Y-projection%*%Y)
    
    dd<-n*q1*(1-p/n)^2
    
    CV.matrix<-matrix(0,n,q)
    for (ii in 1:n)
    {
      for (jj in Y.originalset)
      {
        CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-projection[ii,ii]*G[jj,jj])^2
      }
    }
    cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
    gcv[kk]<-norm(residual,"F")^2/dd
    break
    }
  }
  
  start.point<-mid.lambda+1
  
  for (kk in start.point:number.lambda)
  {
    initial.estimator<-solve(t(X) %*% X+0.001*diag(p)) %*% t(X) %*% Y
    lambda<-lambda.set[kk+1]
    
    for (j in 1:number.group)
    {
      R<-Y-X[,-B[[j]]] %*% initial.estimator[-B[[j]],]
      omega0<-(length(B[[j]])*q1)^(1/2)
      s<-t(X[,B[[j]]]) %*% R/(n*lambda*omega0)
      
      if (norm(s,"F")>1) 
      {
        for (k in Y.originalset)
        {
          for (l in B[[j]]) 
          { 
            Theta[[kk]][l,k]<-(t(X[,l])%*%(Y[,k]-X%*%initial.estimator[,k])/n
                               +initial.estimator[l,k])/(1+lambda*omega0/norm(initial.estimator[B[[j]],],"F"))
          }
        }
      }
      else 
      { 
        for (l in B[[j]]) 
        { 
          Theta[[kk]][l,]<-initial.estimator[l,]
        }
      }
    }
    
    for (i in 1:number.group)
    {
      R<-Y-X[,-B[[i]]]%*%Theta[[kk]][-B[[i]],]
      omega<-(length(B[[i]])*q1)^(1/2)
      ss<-t(X[,B[[i]]]) %*% R/(n*lambda*omega)
      if (norm(ss,"F")>1)
      {
        for (l in B[[i]])
        {x_hat[kk,l]<-1}
      }
      else 
      {
        for (l in B[[i]])
        {
          Theta[[kk]][l,]<-matrix(0,1,q)
        }
      }
    }
    
    I_hat<-c()
    for (i in 1:p){
      if (x_hat[kk,i]==1){I_hat<-cbind(I_hat,i)}
    }
    
  #  print(c(I_hat))
  
    if(length(I_hat)!=0)
    {
      Px<-X[,I_hat]%*%ginv(t(X[,I_hat])%*%X[,I_hat]) %*%t(X[,I_hat])
      residual<-(Y-Px%*%Y)
    
       dd<-n*q1*(1-length(I_hat)/n)^2
    
     CV.matrix<-matrix(0,n,q)
     for (ii in 1:n)
     {
       for (jj in Y.originalset)
       {
         CV.matrix[ii,jj]<-residual[ii,jj]^2/(1-Px[ii,ii]*G[jj,jj])^2
       }
     }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(residual,"F")^2/dd
    }else{
      CV.matrix<-matrix(0,n,q)
      for (ii in 1:n)
      {
        for (jj in Y.originalset)
        {
          CV.matrix[ii,jj]<-Y[ii,jj]^2
        }
      }
      cv[kk]<-sum(CV.matrix, na.rm =TRUE)/(n*q1)
      gcv[kk]<-norm(Y,"F")^2/(n*q1)
      break
      }
  }
  nnumber<-min(which.min(cv),which.min(gcv))
  I_hat<-c()
  for (i in 1:p){
    if (x_hat[nnumber,i]==1){I_hat<-cbind(I_hat,i)}
  }
  
  Beta<-Theta[[nnumber]] %*% solve(hat.sigma)
  
  result<-list(Beta=Beta, x.bestdecision=x_hat[nnumber,], x.bestset=I_hat)
  }else
  {result<-list(Beta=matrix(0,p,q), x.bestdecision=matrix(0,1,p), x.bestset=c())}
  return(result)
}
##########################################################################################################################
# Alternative 1??Lasso-based multiple response regression: Y=X*beta--Lasso penalty is imposed on each element of beta
# Input: X--n*p predictors ; Y--n*q response; nlambda--number of tuning parameters; lambda--sequence of tuning parameter
##########################################################################################################################
multi.Lasso<-function(X,Y,nlambda=50,lambda.min.ratio=0.1,lambda=NULL){
  n=dim(X)[1]
  p=dim(X)[2]
  q=dim(Y)[2]
  
  if (is.null(lambda)){
    lambda_max=matrix(0,q,1)
    lambda_min=matrix(0,q,1)
    for (j in 1:q)	
    {
      fit1<- glmnet(X, Y[,j],nlambda=nlambda); 
      lambda_max[j]<-max(fit1$lambda)
      lambda_min[j]<-min(fit1$lambda)
    }
    lambda_max=max(lambda_max)
    lambda_min=min(lambda_min)
    e=(log(lambda_max)-log(lambda_max*lambda.min.ratio))/(nlambda-1);
    lambda=exp(seq(from=log(lambda_max),to=log(lambda_max*lambda.min.ratio),by=-e))
  }
  
  beta=list()
  for (j in 1:nlambda){
    beta[[j]]=matrix(0,p,q)
  }
  
  for (j in 1:q)	
  {
    fit2<- glmnet(X, Y[,j],lambda=lambda); 
    for (jj in 1:nlambda){
      beta[[jj]][,j]= fit2$beta[,jj]
    }
  }
  
  df=matrix(0,nlambda,1)
  loss=matrix(0,nlambda,1)
  BIC=matrix(0,nlambda,1)
  for (jj in 1:nlambda){
    df[jj]=sum(beta[[jj]]!=0)
    loss[jj]=sum((Y-X%*%beta[[jj]])^2)
  }
  BIC=log(loss)+log(n*q)*df/(n*q)
  
  beta_BIC=beta[[which.min(BIC)]]
  return(list(beta_BIC=beta_BIC,beta=beta,lambda=lambda,BIC=BIC,df=df))
}

#######################################################################################
#Reduced rank stochastic regression with a sparse singular value decomposition developed in Chen et al. (2013)
#input: X--n*p predictors ; Y--n*q response; nlambda--number of tuning parameters; lambda--sequence of tuning parameter; 
#  K--sequence of the rank
#######################################################################################
multi.rrsvd<-function(X,Y,nlambda=50,lambda.min.ratio=1e-2,K.max)
{
  n=dim(X)[1]
  p=dim(X)[2]
  q=dim(Y)[2]
  
  BIC_K=matrix(0,K.max,1)
  
  for (K in 1:K.max){
    
    fit<-rssvd(Y, X, nrank = K,control = list(nlambda=nlambda,minLambda=lambda.min.ratio))
    
    a<-fit$U.path
    
    beta=list()
    for (j in 1:nlambda){
      beta[[j]]=matrix(0,p,q)
    }
    
    df=matrix(0,nlambda,1)
    df_svd=matrix(0,nlambda,1)
    loss=matrix(1e+100,nlambda,1)
    BIC=matrix(1e100,nlambda,1)
    
    for (jj in 1:nlambda){
      temp=dim(fit$U.path)
      if ((length(temp)==3) & (temp[1]==p) & (temp[2]==q) & (temp[3]==K)){
        for (k in 1:K){
          beta[[jj]]= beta[[jj]]+matrix(fit$U.path[,jj,k]*fit$D.path[jj,k],p,1)%*%matrix(fit$V.path[,jj,k],1,q)
        }
        df[jj]=sum(beta[[jj]]!=0)
        df_svd[jj]=sum(fit$U.path[,jj,]!=0)+sum(fit$V.path[,jj,]!=0)-1
        loss[jj]=sum((Y-X%*%beta[[jj]])^2)
      }
    }
    
    BIC=log(loss)+log(n*q)*df_svd/(n*q)
    
    BIC_K[K]=min(BIC)
  }
  K=which.min(BIC_K)
  
  fit<-rssvd(Y, X, nrank = K,control = list(nlambda=nlambda,minLambda=lambda.min.ratio))
  a<-fit$U.path
  
  beta=list()
  for (j in 1:nlambda){
    beta[[j]]=matrix(0,p,q)
  }
  
  for (jj in 1:nlambda){
    for (k in 1:K){
      beta[[jj]]= beta[[jj]]+matrix(fit$U.path[,jj,k]*fit$D.path[jj,k],p,1)%*%matrix(fit$V.path[,jj,k],1,q)
    }
  }
  
  df=matrix(0,nlambda,1)
  df_svd=matrix(0,nlambda,1)
  loss=matrix(0,nlambda,1)
  BIC=matrix(0,nlambda,1)
  for (jj in 1:nlambda){
    df[jj]=sum(beta[[jj]]!=0)
    df_svd[jj]=sum(fit$U.path[,jj,]!=0)+sum(fit$V.path[,jj,]!=0)-1
    loss[jj]=sum((Y-X%*%beta[[jj]])^2)
  }
  
  BIC=log(loss)+log(n*q)*df_svd/(n*q)
  beta_BIC=beta[[which.min(BIC)]]
  
  return(list(beta_BIC=beta_BIC,beta=beta,BIC=BIC,df=df,K=K))
}

###############################################################
# Data   
###############################################################

A1000<-list()

A1000<-list(c(1:10),c(11:20),c(21:30),c(31:40),c(41:50),
            c(51:60),c(61:70),c(71:80),c(81:90),c(91:100),
           c(101:110),c(111:120),c(121:130),c(131:140),c(141:150),
           c(151:160),c(161:170),c(171:180),c(181:190),c(191:200),
           c(201:210),c(211:220),c(221:230),c(231:240),c(241:250),
           c(251:260),c(261:270),c(271:280),c(281:290),c(291:300),
           c(301:310),c(311:320),c(321:330),c(331:340),c(341:350),
           c(351:360),c(361:370),c(371:380),c(381:390),c(391:400),
           c(401:410),c(411:420),c(421:430),c(431:440),c(441:450),
           c(451:460),c(461:470),c(471:480),c(481:490),c(491:500),
           c(501:510),c(511:520),c(521:530),c(531:540),c(541:550),
           c(551:560),c(561:570),c(571:580),c(581:590),c(591:600),
           c(601:610),c(611:620),c(621:630),c(631:640),c(641:650),
           c(651:660),c(661:670),c(671:680),c(681:690),c(691:700),
           c(701:710),c(711:720),c(721:730),c(731:740),c(741:750),
           c(751:760),c(761:770),c(771:780),c(781:790),c(791:800),
           c(801:810),c(811:820),c(821:830),c(831:840),c(841:850),
           c(851:860),c(861:870),c(871:880),c(881:890),c(891:900),
           c(901:910),c(911:920),c(921:930),c(931:940),c(941:950),
           c(951:960),c(961:970),c(971:980),c(981:990),c(991:1000)
           )

nA1000<-100

A2000<-list()

A2000<-list(c(1:10),c(11:20),c(21:30),c(31:40),c(41:50),
            c(51:60),c(61:70),c(71:80),c(81:90),c(91:100),
            c(101:110),c(111:120),c(121:130),c(131:140),c(141:150),
            c(151:160),c(161:170),c(171:180),c(181:190),c(191:200),
            c(201:210),c(211:220),c(221:230),c(231:240),c(241:250),
            c(251:260),c(261:270),c(271:280),c(281:290),c(291:300),
            c(301:310),c(311:320),c(321:330),c(331:340),c(341:350),
            c(351:360),c(361:370),c(371:380),c(381:390),c(391:400),
            c(401:410),c(411:420),c(421:430),c(431:440),c(441:450),
            c(451:460),c(461:470),c(471:480),c(481:490),c(491:500),
            c(501:510),c(511:520),c(521:530),c(531:540),c(541:550),
            c(551:560),c(561:570),c(571:580),c(581:590),c(591:600),
            c(601:610),c(611:620),c(621:630),c(631:640),c(641:650),
            c(651:660),c(661:670),c(671:680),c(681:690),c(691:700),
            c(701:710),c(711:720),c(721:730),c(731:740),c(741:750),
            c(751:760),c(761:770),c(771:780),c(781:790),c(791:800),
            c(801:810),c(811:820),c(821:830),c(831:840),c(841:850),
            c(851:860),c(861:870),c(871:880),c(881:890),c(891:900),
            c(901:910),c(911:920),c(921:930),c(931:940),c(941:950),
            c(951:960),c(961:970),c(971:980),c(981:990),c(991:1000),
            c(1001:1010),c(1011:1020),c(1021:1030),c(1031:1040),c(1041:1050),
            c(1051:1060),c(1061:1070),c(1071:1080),c(1081:1090),c(1091:1100),
            c(1101:1110),c(1111:1120),c(1121:1130),c(1131:1140),c(1141:1150),
            c(1151:1160),c(1161:1170),c(1171:1180),c(1181:1190),c(1191:1200),
            c(1201:1210),c(1211:1220),c(1221:1230),c(1231:1240),c(1241:1250),
            c(1251:1260),c(1261:1270),c(1271:1280),c(1281:1290),c(1291:1300),
            c(1301:1310),c(1311:1320),c(1321:1330),c(1331:1340),c(1341:1350),
            c(1351:1360),c(1361:1370),c(1371:1380),c(1381:1390),c(1391:1400),
            c(1401:1410),c(1411:1420),c(1421:1430),c(1431:1440),c(1441:1450),
            c(1451:1460),c(1461:1470),c(1471:1480),c(1481:1490),c(1491:1500),
            c(1501:1510),c(1511:1520),c(1521:1530),c(1531:1540),c(1541:1550),
            c(1551:1560),c(1561:1570),c(1571:1580),c(1581:1590),c(1591:1600),
            c(1601:1610),c(1611:1620),c(1621:1630),c(1631:1640),c(1641:1650),
            c(1651:1660),c(1661:1670),c(1671:1680),c(1681:1690),c(1691:1700),
            c(1701:1710),c(1711:1720),c(1721:1730),c(1731:1740),c(1741:1750),
            c(1751:1760),c(1761:1770),c(1771:1780),c(1781:1790),c(1791:1800),
            c(1801:1810),c(1811:1820),c(1821:1830),c(1831:1840),c(1841:1850),
            c(1851:1860),c(1861:1870),c(1871:1880),c(1881:1890),c(1891:1900),
            c(1901:1910),c(1911:1920),c(1921:1930),c(1931:1940),c(1941:1950),
            c(1951:1960),c(1961:1970),c(1971:1980),c(1981:1990),c(1991:2000))

nA2000<-200


B100<-list()
B100<-list(c(1:10), c(11:20), c(21:30), c(31:40), (41:50), 
           c(51:60), c(61:70), c(71:80), c(81:90), c(91:100))

nB100<-10           

B200<-list()
B200<-list(c(1:10), c(11:20), c(21:30), c(31:40), (41:50),
           c(51:60), c(61:70), c(71:80), c(81:90), c(91:100), 
           c(101:110), c(111:120), c(121:130),(131:140),c(141:150), 
           c(151:160), c(161:170), c(171:180),c(181:190),c(191:200)) 

nB200<-20           
          
###############################        
# Loop of sample size n  
###############################
for (n in sample.size)
{
  Z<-matrix(rnorm(n*p,0,1),n,p)
  
  gram.orth<-function(U)
  {
    pp<-dim(U)[2]
    V<-matrix(0,pp,pp)
    PC<-matrix(0,pp,pp)
    for (s in 1:pp)
    {
      W<-PC%*%U[,s]
      V[,s]<-(U[,s]-W)/norm((U[,s]-W),"2")
      if (s<pp)
      {
        PC<-PC+V[,s] %*% t(V[,s])
      }
    }
    return(V)
  }
  
  sigma.x<-sample(c(1:10),p,replace=TRUE)
  
  U0<-matrix(rnorm(p*p,0,3),p,p)
  
  Q<-gram.orth(U0)
  
  X<-Z %*% Q %*% diag(sigma.x) %*% t(Q)
  X<-scale(X,center=TRUE,scale=TRUE)

    Coefficient<-matrix(0,p,q)
    
    Theta0<-matrix(runif(p1*q1,min=-theta.h, max=theta.h), p1, q1)
    
    activesetp<-c(1:p1)
    activesetq<-c(1:q1)
    
    Coefficient[activesetp,activesetq]<-Theta0
    
    # Assign initial 0 values 
    
    beta4<-list()
    beta3<-list()
    beta1<-list()
    beta0<-list()
    
    for (m in 1:nreplication)
    {
      beta4[[m]]<-matrix(0,p,q)
      beta3[[m]]<-matrix(0,p,q)
      beta1[[m]]<-matrix(0,p,q)
      beta0[[m]]<-matrix(0,p,q)
    }
    
    SNR<-matrix(0,1,nreplication)
    
    MS0<-matrix(0,1,nreplication)
    TP0<-matrix(0,1,nreplication)
    TN0<-matrix(0,1,nreplication)
    Prec0<-matrix(0,1,nreplication)
    F0<-matrix(0,1,nreplication)
    
    MS0p<-matrix(0,1,nreplication)
    TP0p<-matrix(0,1,nreplication)
    TN0p<-matrix(0,1,nreplication)
    Prec0p<-matrix(0,1,nreplication)
    F0p<-matrix(0,1,nreplication)
    
    MS1<-matrix(0,1,nreplication)
    TP1<-matrix(0,1,nreplication)
    TN1<-matrix(0,1,nreplication)
    Prec1<-matrix(0,1,nreplication)
    F1<-matrix(0,1,nreplication)
    
    MS1p<-matrix(0,1,nreplication)
    TP1p<-matrix(0,1,nreplication)
    TN1p<-matrix(0,1,nreplication)
    Prec1p<-matrix(0,1,nreplication)
    F1p<-matrix(0,1,nreplication)
    
    MS3<-matrix(0,1,nreplication)
    TP3<-matrix(0,1,nreplication)
    TN3<-matrix(0,1,nreplication)
    Prec3<-matrix(0,1,nreplication)
    F3<-matrix(0,1,nreplication)
    
    MS3p<-matrix(0,1,nreplication)
    TP3p<-matrix(0,1,nreplication)
    TN3p<-matrix(0,1,nreplication)
    Prec3p<-matrix(0,1,nreplication)
    F3p<-matrix(0,1,nreplication)
    
    MS4<-matrix(0,1,nreplication)
    TP4<-matrix(0,1,nreplication)
    TN4<-matrix(0,1,nreplication)
    Prec4<-matrix(0,1,nreplication)
    F4<-matrix(0,1,nreplication)
    
    MS4p<-matrix(0,1,nreplication)
    TP4p<-matrix(0,1,nreplication)
    TN4p<-matrix(0,1,nreplication)
    Prec4p<-matrix(0,1,nreplication)
    F4p<-matrix(0,1,nreplication)
    
    SNR<-matrix(0,1,nreplication)
    
    SS0<-norm(X %*% Coefficient,"F")^2
    
    sigma<-sample(c(1:15),q,replace=TRUE)
    
    U1<-matrix(rnorm(q*q,0,2),q,q)
    
    O1<-gram.orth(U1)
    
    Sigma05<-O1 %*% diag(sigma) %*% t(O1)
    
    # Loop of replications    
    for (m in 1:nreplication)
    {
      print("m=")
      print(m)
      
      if (normal==TRUE)
      {
        Error0<-matrix(rnorm(n*q,0,1),n,q) 
      }
      else{
        Error0<-matrix(rt(n*q,3,ncp=0),n,q)
      }
      
      Error<-Error0 %*% Sigma05
      
      SNR[m]<-SS0/norm(Error[,activesetq],"F")^2 
      
      Y<- X %*% Coefficient+Error
    
      inactivesetp<-setdiff(c(1:p),activesetp) 
      inactivesetq<-setdiff(c(1:q),activesetq)
     
      ##########################
# [1] LNC method
      # Selecting active resopnse variables with groups
      
      BS.group<-GRBS(Y, X, nA2000, A2000, 0.00125)
      gy.set<-BS.group$bestset
      
      print("y selection (0)")
      print(c(gy.set))
      
      MS0[m]<-length(gy.set)
      TP0[m]<-length(intersect(activesetq,gy.set))
      TN0[m]<-length(intersect(inactivesetq,setdiff(c(1:q),gy.set)))
      if (TP0[m]!=0)
      {
      Prec0[m]<-TP0[m]/MS0[m]
      F0[m]<-(2/q1)*Prec0[m]*TP0[m]/(Prec0[m]+TP0[m]/q1)
      }
      else{F0[m]<-0}
      
      Gg<-diag(c(BS.group$bestdecision))
      
    # Selecting the active predictors under final RBS with pattern 
      RP.group<-predictor.group(X, Y, Gg, nB100, B100, 200)
      gx.set<-RP.group$x.bestset
      
      print("x selection")
      print(c(gx.set))
      
      MS0p[m]<-length(gx.set)
      TP0p[m]<-length(intersect(activesetp,gx.set))
      TN0p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),gx.set)))
      if (TP0p[m]!=0)
      {Prec0p[m]<-TP0p[m]/MS0p[m]
      F0p[m]<-(2/p1)*Prec0p[m]*TP0p[m]/(Prec0p[m]+TP0p[m]/p1)}
      else{F0p[m]<-0}
      
      for (j in 1:q)
      {
        beta0[[m]][,j]<-RP.group$Beta[,j]
      }
  
     # [2] LNC method

      # To find the response best-subset selector

            BS.single<-RBS(Y,X,0.00125)
            y.set<-BS.single$bestset

            print("y selection")
            print(c(y.set))

            MS1[m]<-length(y.set)
            TP1[m]<-length(intersect(activesetq,y.set))
            TN1[m]<-length(intersect(inactivesetq,setdiff(c(1:q),y.set)))
            if (TP1[m]!=0)
            {Prec1[m]<-TP1[m]/MS1[m]
            F1[m]<-(2/q1)*Prec1[m]*TP1[m]/(Prec1[m]+TP1[m]/q1)
            }else{F1[m]<-0}

            Gs<-diag(c(BS.single$bestdecision))

      # selecting the active prodictor variables

            RP.single<-predictor.choice(X,Y,Gs,200)
            x.set<-RP.single$x.bestset

            print("x selection")
            print(c(x.set))

            MS1p[m]<-length(x.set)
            TP1p[m]<-length(intersect(activesetp,x.set))
            TN1p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),x.set)))
            if (TP1p[m]!=0)
            {
              Prec1p[m]<-TP1p[m]/MS1p[m]
              F1p[m]<-(2/p1)*Prec1p[m]*TP1p[m]/(Prec1p[m]+TP1p[m]/p1)
            }else{F1p[m]<-0}

            for (j in 1:q)
            {
              beta1[[m]][,j]<-RP.single$Beta[,j]
            }
      
      # [3] LNC-bootstrap method
     
      BS.bootstrap<-RBS_bootstrap(Y,X,0.00125)
        #GRBS_bootstrap(Y, X, nA2000, A2000, 0.00125)
      y_lasso<-BS.bootstrap$bestset
      
      print("y selection-bootstrap")
      print(y_lasso)
     
      
      MS3[m]<-length(y_lasso)
      TP3[m]<-length(intersect(activesetq,y_lasso))
      TN3[m]<-length(intersect(inactivesetq,setdiff(c(1:q),y_lasso)))
      if (TP3[m]!=0)
      {
        Prec3[m]<-TP3[m]/MS3[m]
        F3[m]<-(2/q1)*Prec3[m]*TP3[m]/(Prec3[m]+TP3[m]/q1)
      }
      else{F3[m]<-0}
      
      Gs<-diag(c(BS.bootstrap$bestdecision))
      
      RP.bootstrap<-predictor.choice(X,Y,Gs,200)
    #    predictor.group(X, Y, Gg, nB100, B100, 200)
      x_lasso<-RP.bootstrap$x.bestset
      
      print("x selection")
      print(x_lasso)
      
      MS3p[m]<-length(x_lasso)
      TP3p[m]<-length(intersect(activesetp,x_lasso))
      TN3p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),x_lasso)))
      if (TP3p[m]!=0)
      {Prec3p[m]<-TP3p[m]/MS3p[m]
      F3p[m]<-(2/p1)*Prec3p[m]*TP3p[m]/(Prec3p[m]+TP3p[m]/p1)
      }
      else{F3p[m]<-0}
      
      
# # [2] LNC method
# 
# # To find the response best-subset selector
#       
#       BS.single<-RBS(Y,X,0.00125)
#       y.set<-BS.single$bestset
#       
#       print("y selection")
#       print(c(y.set))
#       
#       MS1[m]<-length(y.set)
#       TP1[m]<-length(intersect(activesetq,y.set))
#       TN1[m]<-length(intersect(inactivesetq,setdiff(c(1:q),y.set)))
#       if (TP1[m]!=0)
#       {Prec1[m]<-TP1[m]/MS1[m]
#       F1[m]<-(2/q1)*Prec1[m]*TP1[m]/(Prec1[m]+TP1[m]/q1)
#       }else{F1[m]<-0}
#       
#       Gs<-diag(c(BS.single$bestdecision))
#       
# # selecting the active prodictor variables 
#       
#       RP.single<-predictor.choice(X,Y,Gs,200)
#       x.set<-RP.single$x.bestset
#       
#       print("x selection")
#       print(c(x.set))
#       
#       MS1p[m]<-length(x.set)
#       TP1p[m]<-length(intersect(activesetp,x.set))
#       TN1p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),x.set)))
#       if (TP1p[m]!=0)
#       {
#         Prec1p[m]<-TP1p[m]/MS1p[m]
#         F1p[m]<-(2/p1)*Prec1p[m]*TP1p[m]/(Prec1p[m]+TP1p[m]/p1)
#       }else{F1p[m]<-0}
#       
#       for (j in 1:q)
#       {
#         beta1[[m]][,j]<-RP.single$Beta[,j]
#       }
#       
# # [3] Lasso method
# 
#        fit3<-multi.Lasso(X,Y,nlambda=50,lambda.min.ratio=0.1,lambda=NULL)
#       
#         for (j in c(1:q))
#        {
#          beta3[[m]][,j]<-fit3$beta_BIC[,j]
#        }
#        y_lasso<-c()
#        for (j in 1:q) {
#          if (norm(beta3[[m]][,j],"2") !=0)
#          {y_lasso<-cbind(y_lasso,j)}
#        }
#        x_lasso<-c()
#        for (i in 1:p) {
#          if (norm(beta3[[m]][i,],"2") !=0)
#          {x_lasso<-cbind(x_lasso,i)}
#        }
#        
#       print("y selection")
#       print(y_lasso)
#       print("x selection")
#       print(x_lasso)
#       
#        MS3[m]<-length(y_lasso)
#        TP3[m]<-length(intersect(activesetq,y_lasso))
#        TN3[m]<-length(intersect(inactivesetq,setdiff(c(1:q),y_lasso)))
#        if (TP3[m]!=0)
#        {
#        Prec3[m]<-TP3[m]/MS3[m]
#        F3[m]<-(2/q1)*Prec3[m]*TP3[m]/(Prec3[m]+TP3[m]/q1)
#        }
#        else{F3[m]<-0}
# 
#       MS3p[m]<-length(x_lasso)
#       TP3p[m]<-length(intersect(activesetp,x_lasso))
#       TN3p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),x_lasso)))
#       if (TP3p[m]!=0)
#       {Prec3p[m]<-TP3p[m]/MS3p[m]
#       F3p[m]<-(2/p1)*Prec3p[m]*TP3p[m]/(Prec3p[m]+TP3p[m]/p1)
#       }
#       else{F3p[m]<-0}
#       
# # [4] RSSVD method

      fit4<-multi.rrsvd(X,Y,nlambda=50,lambda.min.ratio=1e-2,3)

      for (j in 1:q)
      {
        beta4[[m]][,j]<-fit4$beta_BIC[,j]
      }
      y_rssvd<-c()
      for (j in 1:q)
      {if (norm(beta4[[m]][,j],"2") !=0) {y_rssvd<-cbind(y_rssvd,j)}}

      print("y selection rssvd")
      print(y_rssvd)

      x_rssvd<-c()
      for (i in 1:p) {if (norm(beta4[[m]][i,],"2") !=0)
      {x_rssvd<-cbind(x_rssvd,i)}}

      print("x selection")
      print(x_rssvd)

      MS4[m]<-length(y_rssvd)
      TP4[m]<-length(intersect(activesetq,y_rssvd))
      TN4[m]<-length(intersect(inactivesetq,setdiff(c(1:q),y_rssvd)))
      if (TP4[m]!=0)
      {Prec4[m]<-TP4[m]/MS4[m]
      F4[m]<-(2/q1)*Prec4[m]*TP4[m]/(Prec4[m]+TP4[m]/q1)}else{F4[m]<-0}

      MS4p[m]<-length(x_rssvd)
      TP4p[m]<-length(intersect(activesetp,x_rssvd))
      TN4p[m]<-length(intersect(inactivesetp,setdiff(c(1:p),x_rssvd)))
      if (TP4p[m]!=0)
      {
        Prec4p[m]<-TP4p[m]/MS4p[m]
        F4p[m]<-(2/p1)*Prec4p[m]*TP4p[m]/(Prec4p[m]+TP4p[m]/p1)
      }else{F4p[m]<-0}

    }
    
  # Averaging absolute replication means and sd of all theta[i,j]
         
    Tbeta4.bias<-matrix(0,p,q)
    Tbeta3.bias<-matrix(0,p,q)
    Tbeta1.bias<-matrix(0,p,q)
    Tbeta0.bias<-matrix(0,p,q)
    
    beta4.bias<-matrix(0,p,q)
    beta3.bias<-matrix(0,p,q)
    beta1.bias<-matrix(0,p,q)
    beta0.bias<-matrix(0,p,q)
    
    beta4.sse<-matrix(0,p,q)
    beta3.sse<-matrix(0,p,q)
    beta1.sse<-matrix(0,p,q) 
    beta0.sse<-matrix(0,p,q)
    
    beta4.mse<-matrix(0,p,q)
    beta3.mse<-matrix(0,p,q)
    beta1.mse<-matrix(0,p,q)
    beta0.mse<-matrix(0,p,q)
     
    for (i in 1:p)
     {
      for (j in 1:q)
      {
       for (m in 1:nreplication)
       {
        Tbeta4.bias[i,j]<-Tbeta4.bias[i,j]+beta4[[m]][i,j]-Coefficient[i,j]
        Tbeta3.bias[i,j]<-Tbeta3.bias[i,j]+beta3[[m]][i,j]-Coefficient[i,j]
        Tbeta1.bias[i,j]<-Tbeta1.bias[i,j]+beta1[[m]][i,j]-Coefficient[i,j]
        Tbeta0.bias[i,j]<-Tbeta0.bias[i,j]+beta0[[m]][i,j]-Coefficient[i,j]
 
        beta4.sse[i,j]<-beta4.sse[i,j]+(beta4[[m]][i,j]-Coefficient[i,j])^2
        beta3.sse[i,j]<-beta3.sse[i,j]+(beta3[[m]][i,j]-Coefficient[i,j])^2
        beta1.sse[i,j]<-beta1.sse[i,j]+(beta1[[m]][i,j]-Coefficient[i,j])^2
        beta0.sse[i,j]<-beta0.sse[i,j]+(beta0[[m]][i,j]-Coefficient[i,j])^2
       }
        
        beta4.bias[i,j]<-Tbeta4.bias[i,j]/nreplication
        beta3.bias[i,j]<-Tbeta3.bias[i,j]/nreplication
        beta1.bias[i,j]<-Tbeta1.bias[i,j]/nreplication
        beta0.bias[i,j]<-Tbeta0.bias[i,j]/nreplication
        
        beta4.mse[i,j]<-(beta4.sse[i,j]/nreplication)
        beta3.mse[i,j]<-(beta3.sse[i,j]/nreplication)
        beta1.mse[i,j]<-(beta1.sse[i,j]/nreplication)
        beta0.mse[i,j]<-(beta0.sse[i,j]/nreplication)
      }
    }
    
  AbsBias_gLNC<-mean(abs(beta0.bias))
  AbsBias_lnc<-mean(abs(beta1.bias))
  AbsBias_lasso<-mean(abs(beta3.bias))
  AbsBias_rssvd<-mean(abs(beta4.bias))
  
  AveMSE_gLNC<-mean(beta0.mse)
  AveMSE_lnc<-mean(beta1.mse)
  AveMSE_lasso<-mean(beta3.mse)
  AveMSE_rssvd<-mean(beta4.mse)
  
# Means of replications
  
  gLNC.mean<-c(mean(TP0)/q1, mean(TN0)/(q-q1), mean(F0), mean(TP0p)/p1, 
              mean(TN0p)/(p-p1), mean(F0p), AbsBias_gLNC, AveMSE_gLNC)
  
  LNC.mean<-c(mean(TP1)/q1, mean(TN1)/(q-q1), mean(F1), mean(TP1p)/p1, 
                mean(TN1p)/(p-p1), mean(F1p), AbsBias_lnc, AveMSE_lnc)
  
 Lasso.mean<-c(mean(TP3)/q1, mean(TN3)/(q-q1), mean(F3),mean(TP3p)/p1, 
                mean(TN3p)/(p-p1), mean(F3p), AbsBias_lasso, AveMSE_lasso)
  
 rssvd.mean<-c(mean(TP4)/q1, mean(TN4)/(q-q1), mean(F4), mean(TP4p)/p1, 
                 mean(TN4p)/(p-p1), mean(F4p), AbsBias_rssvd, AveMSE_rssvd)
  
  SNR.mean<-mean(SNR)
  
  
  list_data<-list(c(n,q,q1,p,p1,theta.h), 
                  c(round(SNR.mean, digits=4)),
                  "gLNC", round(gLNC.mean, digits=4),
                  "LNC", round(LNC.mean, digits=4),
                  "Multi_Lasso", round(Lasso.mean, digits=4),
                  "SRRR", round(rssvd.mean, digits=4))
  
  print(list_data)
  newdata<-rbind(newdata,list_data)
  write.matrix(newdata,file="~/Documents/JMVA/nonoverlaping.csv", sep = "\n")
}
