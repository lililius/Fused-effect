#the number of subject
m=50 
# the number of observation of each subject
ni=sample(c(1,2),m,replace=TRUE,prob=c(1/2,1/2))
N=sum(ni) #the number of total observations
X=rnorm(N,0,1)
beta=2
eps=rnorm(N,0,0.4)
Z=matrix(0,N,m)
for(k in 1:m){
if(k==1){
Z[(1:ni[1]),1]=rep(1,ni[1])
}else{
Z[(sum(ni[(1:m)<k])+1):(sum(ni[(1:m)<k])+ni[k]),k]=rep(1,ni[k])
}
}
atrue=sample(c(-1.5,0,1.5),m,replace=TRUE,prob=c(1/3,1/3,1/3))
eps=rnorm(N,0,0.4)
Y=as.vector(Z%*%atrue)+as.vector(X*beta)+eps
id=NULL
for(k in 1:m){
id=c(id,rep(k,ni[k]))
}
data=as.data.frame(cbind(id,X,Y))
print(data)
Fused_effect(data$X,data$Y,data$id,seq(0.14,0.35,length=10))

