rm(list=ls())
library(energy)
library(rrcov)
 

library(PMVANOVA)


n =30;
p =10;  
 
example=52
Times =1000;
BstrpTimes = 199;
pVal =matrix(0,Times,1);
pValue =matrix(0,3,1);
 
discpva1=matrix(0,Times,1);
discpva1.mean =matrix(0,3,1);
 
discpva2=matrix(0,Times,1);
discpva2.mean =matrix(0,3,1);
 
discpva3=matrix(0,Times,1);
discpva3.mean =matrix(0,3,1);
 
discpva4=matrix(0,Times,1);
discpva4.mean =matrix(0,3,1);
 
 

 
delta=c(0)
Num.na=matrix(0,1,4);


###col name: methods, 0.01,0.05,0.1
em.power=array(0,dim=c(5, 3, length(delta)))

for (j in 1:length(delta)){

delta.j=delta[j]
cat("delta.j= ", delta.j,"\n")
num=floor(n*delta.j)


for (ii in 1:Times){


if (example==51){
x11 = matrix(rnorm(p*n*delta.j,0,1),ncol = p)     #group 1 
x12 = matrix(rcauchy(p*n*(1-delta.j),0,1),ncol = p)     #group 1 
x1=rbind(x11,x12)

x2 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 2  
x3 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 3
x4 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 4
 }
 

if (example==52){

x11 = matrix(rcauchy(p*num,0,1),ncol = p)         #group 1 
x12 = matrix(rcauchy(p*(n-num),0,1),ncol = p)     #group 1 
x1=rbind( (x11),exp(x12))  

x2 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 2  
x3 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 3
x4 = matrix(rcauchy(p*n,0,1),ncol = p)           #group 4
 
x2 =exp(x2)
x3 =exp(x3)
x4 =exp(x4)
} 


X=rbind(x1,x2,x3,x4)
Y=c(rep(0,n),rep(1,n),rep(2,n),rep(3,n))      #group numbers
 
 
cat("iter = ", ii,"\n")

###PPMV: our proposed method 
ppmv=test_ppmv(X,Y,B=BstrpTimes)
      pVal[ii,1]=ppmv$pvalue
      pValue[1,] = apply(pVal[1:ii,,drop=FALSE] <=  0.01,2,mean);
      pValue[2,] = apply(pVal[1:ii,,drop=FALSE] <=  0.05,2,mean);
      pValue[3,] = apply(pVal[1:ii,,drop=FALSE] <=  0.1,2,mean);
      cat("ppmv.test = ", pValue,"\n")


###DISCO:  Rizzo and   Szekely (2010).    
dis=disco(X, factors=Y, index=1, R=BstrpTimes)
     discpva1[ii,1]=as.numeric(dis$p.value) 
     discpva1.mean[1,] = mean(discpva1[1:ii,,drop=FALSE] <= 0.01, na.rm = T);
     discpva1.mean[2,] = mean(discpva1[1:ii,,drop=FALSE] <= 0.05, na.rm = T);
     discpva1.mean[3,] = mean(discpva1[1:ii,,drop=FALSE] <= 0.1, na.rm = T);
     cat("disco_1 = ", discpva1.mean,"\n")

 

dis2<- disco(X, factors=Y, index=0.5, R=BstrpTimes)
discpva2[ii,1]=as.numeric(dis2$p.value) 
     discpva2.mean[1,] = mean(discpva2[1:ii,,drop=FALSE] <= 0.01, na.rm = T);
     discpva2.mean[2,] = mean(discpva2[1:ii,,drop=FALSE] <= 0.05, na.rm = T);
     discpva2.mean[3,] = mean(discpva2[1:ii,,drop=FALSE] <= 0.1, na.rm = T);
     cat("disco_0.5= ", discpva2.mean,"\n")

 
 
dis3<-disco(X, factors=Y, index=0.2, R=BstrpTimes)
discpva3[ii,1]=as.numeric(dis3$p.value) 
     discpva3.mean[1,] = mean(discpva3[1:ii,,drop=FALSE] <= 0.01, na.rm = T);
     discpva3.mean[2,] = mean(discpva3[1:ii,,drop=FALSE] <= 0.05, na.rm = T);
     discpva3.mean[3,] = mean(discpva3[1:ii,,drop=FALSE] <= 0.1, na.rm = T);
     cat("disco_0.2= ", discpva3.mean,"\n")


dis4<-disco(X, factors=Y, index=0.02, R=BstrpTimes)
discpva4[ii,1]=as.numeric(dis4$p.value) 
     discpva4.mean[1,] = mean(discpva4[1:ii,,drop=FALSE] <= 0.01, na.rm = T);
     discpva4.mean[2,] = mean(discpva4[1:ii,,drop=FALSE] <= 0.05, na.rm = T);
     discpva4.mean[3,] = mean(discpva4[1:ii,,drop=FALSE] <= 0.1, na.rm = T);
     cat("disco_0.02= ", discpva4.mean,"\n")

}


cat("########### ", "\n")
##col name: methods, 0.01,0.05,0.1
em.power[,,j]=rbind(rbind(t(pValue),t(discpva1.mean)),rbind(t(discpva2.mean),
   t(discpva3.mean)),t(discpva4.mean))
 

Num.na[j,1]=sum(is.na(discpva1[1:ii,,drop=FALSE]))/ii
Num.na[j,2]=sum(is.na(discpva2[1:ii,,drop=FALSE]))/ii
Num.na[j,3]=sum(is.na(discpva3[1:ii,,drop=FALSE]))/ii
Num.na[j,4]=sum(is.na(discpva4[1:ii,,drop=FALSE]))/ii


}
 
em.power

         