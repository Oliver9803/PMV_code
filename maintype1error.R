#set.seed(1)
library(energy)
library(rrcov)

library(PMVANOVA)


n =50;
p =10;  
 

example=33
Times =1000;
BstrpTimes = 199;
pVal =matrix(0,Times,1);
pValue =matrix(0,3,1);
drpVal =matrix(0,Times,1);
drpValue =matrix(0,3,1);
dcortpVal =matrix(0,Times,1);
dcortpValue =matrix(0,3,1);
wilrankpVal =matrix(0,Times,1);
wilrankpValue =matrix(0,3,1);

for (ii in 1:Times){

if (example==11){
delta.j=0
x1 = matrix(rt(p*n,4,delta.j),ncol = p)     #group 1 
x2 = matrix(rt(p*n,4,0),ncol = p)           #group 2  
x3 = matrix(rt(p*n,4,0),ncol = p)           #group 3
x4 = matrix(rt(p*n,4,0),ncol = p)           #group 4
 }

if (example==22){
delta.j=0
x11 = matrix(rcauchy(p*n/2,-delta.j,1),ncol = p)     #group 1 
x12 = matrix(rcauchy(p*n/2,delta.j,1),ncol = p)     #group 1 
x1=rbind(x11,x12)

x2 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 2  
x3 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 3
x4 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 4
 }
 

if (example==33){
delta.j=1
x1 = matrix(rcauchy(p*n,0,delta.j),ncol = p)     #group 1 
x2 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 2  
x3 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 3
x4 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 4
 }



if (example==44){
delta.j=0
x11 = matrix(rnorm(p*n*delta.j,0,1),ncol = p)     #group 1 
x12 = matrix(rcauchy(p*n*(1-delta.j),0,1),ncol = p)     #group 1 
x1=rbind(x11,x12)

x2 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 2  
x3 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 3
x4 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 4
 }


X=rbind(x1,x2,x3,x4)
Y=c(rep(0,n),rep(1,n),rep(2,n),rep(3,n)) #group numbers
 
 
cat("iter = ", ii,"\n")

###PPMV: our proposed method 
ppmv=test_ppmv(X,Y,B=BstrpTimes)
      pVal[ii,1]=ppmv$pvalue
      pValue[1,] = apply(pVal[1:ii,,drop=FALSE] <=  0.01,2,mean);
      pValue[2,] = apply(pVal[1:ii,,drop=FALSE] <=  0.05,2,mean);
      pValue[3,] = apply(pVal[1:ii,,drop=FALSE] <=  0.1,2,mean);
      cat("ppmv.test = ", pValue,"\n")


###DISCO:  Rizzo and   Szekely (2010).    
dis=disco(X, factors=Y, R=BstrpTimes)
      drpVal[ii,1] =as.numeric(dis$p.value);
      drpValue[1,] = apply(drpVal[1:ii,,drop=FALSE] <=  0.01,2,mean);
      drpValue[2,] = apply(drpVal[1:ii,,drop=FALSE] <=  0.05,2,mean);
      drpValue[3,] = apply(drpVal[1:ii,,drop=FALSE] <=  0.1,2,mean);
      cat("disco.test = ", drpValue,"\n")

 
###  One-way MANOVA: Wilks Lambda 
fit<-manova(X~Y)
dcortpVal[ii,1]=summary(fit)$stats[1, "Pr(>F)"]

     dcortpValue[1,] = apply(dcortpVal[1:ii,,drop=FALSE] <= 0.01,2,mean);
     dcortpValue[2,] = apply(dcortpVal[1:ii,,drop=FALSE] <= 0.05,2,mean);
     dcortpValue[3,] = apply(dcortpVal[1:ii,,drop=FALSE] <= 0.1,2,mean);
     cat("wilks.test = ", dcortpValue,"\n")


###Robust One-way MANOVA: Wilks Lambda 
wilrank=Wilks.test(X, grouping=Y,  method="rank")
wilrankpVal[ii,1] =wilrank$p.value;
     wilrankpValue[1,] = apply(wilrankpVal[1:ii,,drop=FALSE] <= 0.01,2,mean);
     wilrankpValue[2,] = apply(wilrankpVal[1:ii,,drop=FALSE] <= 0.05,2,mean);
     wilrankpValue[3,] = apply(wilrankpVal[1:ii,,drop=FALSE] <= 0.1,2,mean);
     cat("rank.test = ", wilrankpValue,"\n")

}

 

err.tpe1=rbind(rbind(t(pValue),t(drpValue)),rbind(t(dcortpValue),t(wilrankpValue)))

err.tpe1


