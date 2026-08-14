#######################################################################################
# 1. Group response best subset selector: 
#######################################################################################
# Y----nxq observation matrix, 
# alpha----significance level,
# X----nxp predictor matrix, p<n is required,
# number.groups----the number of response groups,
# A---- Response group list,
# steplength.alpha----searching step length of alpha in (0,0.05], 0.00125 is recommended.
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
# 2. Predictor-variable selection with group patern:
#############################################################################

# Y----nxq observation matrix, 
# X----nxp predictor matrix, p<n is required,
# G=diag(d_1,...,d_q), d_i=0,1 indicator, d_i=1 means y_i selected,
# number.group-----the number of predictor groups,
# B----predictor group list,
# number.lambda----the number of lambda raking values,100 is reccomended.

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



#######################################################################################
# 3. Response best subset selector (RBS)
#######################################################################################
# Y----nxq observation matrix, 
# X----nxp predictor matrix, p<n is required,
# steplength.alpha----searching step length of alpha in (0,0.05], 0.00125 is recommended.
    
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




###################################
# 4. Predictor variable selection:
###################################
# X----nxp predictor matrix, 
# Y----nxq observation matrix, 
# G=diag(d_1,...,d_q), d_i=0,1 indicator, d_i=1 means y_i selected,
# number.lambda----the number of lambda taking values, for example, 100,

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
