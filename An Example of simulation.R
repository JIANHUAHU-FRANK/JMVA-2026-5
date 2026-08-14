#install.packages("glmnet")
#install.packages("rrpack") 
#install.packages("MASS") 
#install.packages("Matrix") 

library(Matrix)
library(glmnet)
library(rrpack)
library(MASS)

nreplication<-200

sample.size<-c(200)         

q<-2000     # The number of response variables
q1<-10     # The number of active response variables 
p<-100      # The number of predictor variables
p1<-10      # Then number of active predictor variables

theta.h<-6.5   # normal, n=200, 8.5(2.7551),7.5()
normal<-TRUE    # t(3), n=200, 10(1.5580), 8.5(0.8427)


newdata<-c()

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
    
    beta1<-list()
    beta0<-list()
    
    
    ###############################        
    # Loop for replication  
    ###############################
    
    
    for (m in 1:nreplication)
    {
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
     
      
      ################################
      # [1] gLNC method
      ################################
      
      # The 1st step: Selecting active resopnse variables with groups
      
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
      
    # The 2nd step: Selecting the active predictors under final RBS with pattern 
      
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
  
      
     ##########################
     # [2] LNC method
     ##########################
      
     # The 1st step: finding the response best-subset selector

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

      # The 2nd step: selecting the active predictor variables

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
    }
    
  # loop end
     
    
  # Averaging absolute replication means and sd of all theta[i,j]
         

    Tbeta1.bias<-matrix(0,p,q)
    Tbeta0.bias<-matrix(0,p,q)
    
    beta1.bias<-matrix(0,p,q)
    beta0.bias<-matrix(0,p,q)
    
  
    beta1.sse<-matrix(0,p,q) 
    beta0.sse<-matrix(0,p,q)
    
  
    beta1.mse<-matrix(0,p,q)
    beta0.mse<-matrix(0,p,q)
     
    for (i in 1:p)
     {
      for (j in 1:q)
      {
       for (m in 1:nreplication)
       {
        Tbeta1.bias[i,j]<-Tbeta1.bias[i,j]+beta1[[m]][i,j]-Coefficient[i,j]
        Tbeta0.bias[i,j]<-Tbeta0.bias[i,j]+beta0[[m]][i,j]-Coefficient[i,j]
 
        beta1.sse[i,j]<-beta1.sse[i,j]+(beta1[[m]][i,j]-Coefficient[i,j])^2
        beta0.sse[i,j]<-beta0.sse[i,j]+(beta0[[m]][i,j]-Coefficient[i,j])^2
       }
        
        beta1.bias[i,j]<-Tbeta1.bias[i,j]/nreplication
        beta0.bias[i,j]<-Tbeta0.bias[i,j]/nreplication
        
        beta1.mse[i,j]<-(beta1.sse[i,j]/nreplication)
        beta0.mse[i,j]<-(beta0.sse[i,j]/nreplication)
      }
    }
    
  AbsBias_gLNC<-mean(abs(beta0.bias))
  AbsBias_lnc<-mean(abs(beta1.bias))

  
  AveMSE_gLNC<-mean(beta0.mse)
  AveMSE_lnc<-mean(beta1.mse)
  
  
# Means of replications
  
  gLNC.mean<-c(mean(TP0)/q1, mean(TN0)/(q-q1), mean(F0), mean(TP0p)/p1, 
              mean(TN0p)/(p-p1), mean(F0p), AbsBias_gLNC, AveMSE_gLNC)
  
  LNC.mean<-c(mean(TP1)/q1, mean(TN1)/(q-q1), mean(F1), mean(TP1p)/p1, 
                mean(TN1p)/(p-p1), mean(F1p), AbsBias_lnc, AveMSE_lnc)
  
  SNR.mean<-mean(SNR)
  
  
  list_data<-list(c(n,q,q1,p,p1,theta.h), 
                  c(round(SNR.mean, digits=4)),
                  "gLNC", round(gLNC.mean, digits=4),
                  "LNC", round(LNC.mean, digits=4))
  
  print(list_data)
  newdata<-rbind(newdata,list_data)
  write.matrix(newdata,file="~/Documents/JMVA/simulationResults.csv", sep = "\n")

