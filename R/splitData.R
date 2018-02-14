splitData <-
function(response, array,fac){
  tarray<-t(array)
  tarray2<-as.data.frame(tarray[order(response),])
  response2<-response[order(response)]
  nspl<-split(tarray2,response2)
  
  ns = unlist(lapply(nspl,nrow))
  nr = round(ns*fac,0)
  fnres=lapply(1:length(nspl),function(x){
    tres=nspl[[x]]
    rsamp<-sample(c(1:ns[x]),nr[x])
    tr=tres[rsamp,]
    lr=tres[-rsamp,]
    list(tr=tr,lr=lr)
  })
  
  TR=do.call(rbind,lapply(fnres[1:length(nspl)],function(x){x['tr'][[1]]}))
  LR=do.call(rbind,lapply(fnres[1:length(nspl)],function(x){x['lr'][[1]]}))
  
  cid<-unique(response2)
  res<-unlist(lapply(1:length(nspl),function(x) {rep(cid[x],nr[x])}))
  res2<-unlist(lapply(1:length(nspl),function(x) {rep(cid[x],(ns[x]-nr[x]))}))
  TR2<-cbind(res,TR)
  LR2<-cbind(res2,LR)
  data<-list(TR2,LR2)
  names(data)<-c('Learning_set','Test_set')
  return(data)
}
