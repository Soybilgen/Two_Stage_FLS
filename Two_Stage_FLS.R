# installing required packages
rm(list=ls(all=TRUE))
library("FKF")
library("foreach")
library("doParallel")
library('glmnet')

# loading required functions
source("Two_Stage_FLS_Functions.R")

# load the data file
RAW<-read.csv(file="data.csv",sep = ",", dec=".")

# FLS smoothing parameter
mu0.scale=100

set.seed(1)

# independent variables: output gap, inflation gap, real exchange rate gap
x=data.matrix(RAW[,c(3,4,5)], rownames.force = NA)
# dependent variable: real rate target
y=data.matrix(RAW$R, rownames.force = NA)

T=dim(x)[1]

# First stage lasso+OLS
xend=x

Z=cbind(lagmatrix(y,rbind(matrix(seq(2,6,1),nrow=5),9, 12))[13:T,],lagmatrix(x,rbind(matrix(seq(1,6,1),nrow=6),9, 12))[13:T,])
fitted.ls=matrix(rep(0,(3*(dim(Z)[1]))),nrow=dim(Z)[1])
indx.all={}
for (j in 1:dim(xend)[2]){
  obj.escv <- cv.glmnet(Z, xend[13:T,j],alpha=1)
  obj <- glmnet(Z, xend[13:T,j],alpha=1,lambda=obj.escv$lambda.min)
  indx=which(matrix(drop(obj$beta))!=0)
  indx.all=c(indx.all,indx)
  obj.ls=lm(xend[13:T,j]~Z[,indx])
  fitted.ls[,j]=obj.ls$fitted.values
}
indx.all=unique(indx.all)
Z.all=Z[,indx.all]

# Second stage Kalman
Yn=matrix(y[13:T],nrow=(T-12))
Yn.lag1=matrix(y[12:(T-1)],nrow=(T-12))
Xn=cbind(Yn.lag1,fitted.ls)
nvar=1
cx=dim(Xn)[2]
T=dim(Xn)[1]
ls.result=lm(Yn~Xn)

Xnf=cbind(1,Xn)
mu0=mu0.scale
bfls=fls(Xnf, Yn, mu=mu0)
bfls0=t(bfls)
res=Yn-bfls0[,1]-bfls0[,2]*Xnf[,2]-bfls0[,3]*Xnf[,3]-bfls0[,4]*Xnf[,4]-bfls0[,5]*Xnf[,5]
obj.ls.oir=lm(res~Z.all)
results<-bbfls(Xnf,Yn,mu0) 

# results[[1]]:bfls, results[[2]]: std(bfls), results[[3]]: 
# lower conf int, results[[4]]: upper confidence interval
coeff <- results[[1]]
coeff_ub <- results[[4]]
coeff_lb <- results[[3]]
