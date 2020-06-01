library(energy)
library(rrcov)

library(PMVANOVA)

  
load("michigandata.RData")
X=t(michigandata)
Y=c(rep(1,24),rep(2,62))
  
#####################
Y=factor(Y)
X=as.matrix(X)
n=length(Y)
BstrpTimes = 999;


####################################################################
######PPMV: our proposed method ######################################### 
####################################################################
ppmv=test_ppmv(as.matrix(X),Y,B=BstrpTimes)
       ppmv$pvalue
      
 

sum.factor=ppmv$ppmv*n^2
mean.factor=ppmv$ppmv*n^2/( length(levels(Y))-1)
sum.tot=TSS(X)*n^2
sum.bet=(sum.tot-sum.factor)
mean.bet=(sum.tot-sum.factor)/(n- length(levels(Y)))
F.ratio=mean.factor/mean.bet

 

####Result PMV test############
 
cat( "\n", "df.fac=", length(levels(Y))-1, " ", "sum.fa=", sum.factor, " ",
 "mean.fac=", mean.factor," ", "F.rat=", F.ratio," ", "pvalue=", ppmv$pvalue, "\n","\n",
 "df.bet=", n-length(levels(Y)), " ", "sum.bet=", sum.bet, " ", "mean.bet=", mean.bet, "\n","\n",
 "df.tot=", n-1, " ", "sum.tot=", sum.tot, "\n","\n")


 
####################################################################
###DISCO:  Rizzo and   Szekely (2010). 
####################################################################   
dis =disco(X, factors=Y, R=BstrpTimes)
 
dis 

  
 