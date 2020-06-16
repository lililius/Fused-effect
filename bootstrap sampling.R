boot.cluster=function(x,id){
  boot.id=sample(unique(id),replace=T)
  out=lapply(1:length(boot.id), function(newid){cbind(x[which(id==boot.id[newid]),],newid)})
  return(do.call("rbind",out))
}
