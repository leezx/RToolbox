# repeat the col or row for n times
# example: rep.df <- rep.col(all_tsne[,cluster], length(groups))
rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

