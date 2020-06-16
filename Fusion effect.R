fusion_effect=function(X,Y,id,lambda_grid){
library(MASS)
library(mclust)#adjustRI
library(nlme)
Sn=function(x,t){
xx=1-t/sqrt(sum(x^2))
if(xx>0){
s=xx*x
}else{
s=0*x
}
list(s=s)
}

Theta_MCP=function(u,lambda){
L=length(u)
Ta=rep(0,L)
for(k in 1:L){
uabs=abs(u[k])
if(uabs<=(v*lambda)){
Ta[k]=Sn(x=u[k],t=lambda/eta)$s/(1-1/(v*eta))
}else{
Ta[k]=u[k]
}
}
list(Ta=Ta)
}

Theta_SCAD=function(u,lambda){
L=length(u)
Ta=rep(0,L)
for(k in 1:L){
uabs=abs(u[k])
if(uabs<=(lambda+lambda/eta)){
Ta[k]=Sn(x=u[k],t=lambda/eta)$s
}else if(uabs>(lambda+lambda/eta)&uabs<=(v*lambda)){
Ta[k]=Sn(x=u[k],t=v*lambda/((v-1)*eta))$s/(1-1/((v-1)*eta))
}else{
Ta[k]=u[k]
}
}
list(Ta=Ta)
}

Fn_SCAD=function(lambda){
r=rep(0,(m-1)*m/2)
theta=theta0
upsilon=upsilon0
counter=1
repeat{
ahat=solve(t(Z)%*%Qx%*%Z+eta*t(Delta)%*%Delta)%*%(t(Z)%*%Qx%*%Y+eta*t(Delta)%*%(theta-1/eta*upsilon))
#ahat=as.vector(ahat)
betahat=solve(t(X)%*%X)%*%t(X)%*%(Y-Z%*%ahat)
for(i in 1:(m-1)){
Pi=ahat[i]-ahat[(1:m)>i]+1/eta*upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]
theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=Theta_SCAD(Pi,lambda)$Ta
upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]+eta*(ahat[i]-ahat[(1:m)>i]-theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)])
r[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=ahat[i]-ahat[(1:m)>i]-theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]
}
rnorm=sqrt(sum(r^2))
if(rnorm<0.00001){
break
}
if(counter>1000){
break
}
counter=counter+1
}
list(ahat=ahat,betahat=betahat)
}

Fn_MCP=function(lambda){
r=rep(0,(m-1)*m/2)
theta=theta0
upsilon=upsilon0
counter=1
repeat{
ahat=solve(t(Z)%*%Qx%*%Z+eta*t(Delta)%*%Delta)%*%(t(Z)%*%Qx%*%Y+eta*t(Delta)%*%(theta-1/eta*upsilon))
betahat=solve(t(X)%*%X)%*%t(X)%*%(Y-Z%*%ahat)
for(i in 1:(m-1)){
Pi=ahat[i]-ahat[(1:m)>i]+1/eta*upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]
theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=Theta_MCP(Pi,lambda)$Ta
upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=upsilon[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]+eta*(ahat[i]-ahat[(1:m)>i]-theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)])
r[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]=ahat[i]-ahat[(1:m)>i]-theta[((i-1)*(2*m-i)/2+1):((i-1)*(2*m-i)/2+m-i)]
}
rnorm=sqrt(sum(r^2))
if(rnorm<0.00001){
break
}
if(counter>1000){
break
}
counter=counter+1
}
list(ahat=ahat,betahat=betahat)
}

p=ifelse(is.matrix(X),dim(X)[2],1)
N=length(Y)
m=id[N]
ni=rep(0,m)
for(j in 1:m){
ni[j]=length(which(id==j))
}
Z=matrix(0,N,m)
for(k in 1:m){
if(k==1){
Z[(1:ni[1]),1]=rep(1,ni[1])
}else{
Z[(sum(ni[(1:m)<k])+1):(sum(ni[(1:m)<k])+ni[k]),k]=rep(1,ni[k])
}
}
eta=1
v=3

#Initialize theta and upsilon
li=c(1:m)
all=Z%*%li
data_ran=data.frame(Y=Y,all=all,X=X)
model=lme(Y~X,random=~1|all,data=data_ran)
a0=model$coef$random$all+model$coef$fixed[1]
up0=rep(0,m)
theta0=NULL
upsilon0=NULL
for(i in 1:(m-1)){
theta0=c(theta0,(a0[i]-a0[(1:m)>i]))
upsilon0=c(upsilon0,(up0[i]-up0[(1:m)>i]))
}
Im=diag(m)
Delta=NULL
for(i in 1:(m-1)){
Delta=cbind(Delta,(Im[,i]-Im[,(1:m)>i]))
}
Delta=t(Delta)
xx=solve(t(X)%*%X)
Qx=diag(N)-X%*%xx%*%t(X)

#the grid search of lambda
ngrid=length(lambda_grid)
#SCAD
Group=rep(0,ngrid)
bic=rep(0,ngrid)
b1=rep(0,ngrid)
b2=rep(0,ngrid)
ascad=list()
betascad=list()
for(i in 1:ngrid){
value=Fn_SCAD(lambda=lambda_grid[i])
ascad[[i]]=value$ahat
betascad[[i]]=value$betahat
z1=round(ascad[[i]],3)
ta1=as.data.frame(table(z1))
print(ta1)
Group[i]=dim(ta1)[1]
print(Group)
R=Y-Z%*%ascad[[i]]-X%*%betascad[[i]]
b1[i]=log(sum(R^2)/N)
Cn=5*log(log(N+p))
b2[i]=log(N)/N*(Group[i]+p)
bic[i]=b1[i]+Cn*b2[i]
}
plot(lambda_grid,bic)
I=which.min(bic)
aest=ascad[[I]]
betaest=betascad[[I]]
z=round(ascad[[I]],3)
ta=as.data.frame(table(z))
ngroup=dim(ta)[1]
print(ngroup)


#MCP
Group_mcp=rep(0,ngrid)
bic_mcp=rep(0,ngrid)
b3=rep(0,ngrid)
b4=rep(0,ngrid)
amcp=list()
betamcp=list()
for(i in 1:ngrid){
value_mcp=Fn_MCP(lambda=lambda_grid[i])
amcp[[i]]=value_mcp$ahat
betamcp[[i]]=value_mcp$betahat
z2=round(amcp[[i]],3)
ta2=as.data.frame(table(z2))
print(ta2)
Group_mcp[i]=dim(ta2)[1]
print(Group_mcp)
R_mcp=Y-Z%*%amcp[[i]]-X%*%betamcp[[i]]
b3[i]=log(sum(colSums(R_mcp^2))/N)
Cn=5*log(log(N+p))
b4[i]=log(N)/N*(Group_mcp[i]+p)
bic_mcp[i]=b3[i]+Cn*b4[i]
}
plot(lambda_grid,bic_mcp)
I1=which.min(bic_mcp)
aest_mcp=amcp[[I1]]
betaest_mcp=betamcp[[I1]]
z_mcp=round(amcp[[I1]],3)
ta_mcp=as.data.frame(table(z_mcp))
ngroup_mcp=dim(ta_mcp)[1]	
print(ngroup_mcp)
result=list(ngroup=ngroup,betaest=betaest,aest=aest,ngroup_mcp=ngroup_mcp,betaest_mcp=betaest_mcp,aest_mcp=aest_mcp)
return(result)
}


