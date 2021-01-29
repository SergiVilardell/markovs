library(mixtools)
ids<-0:9
id<-0
x<-sort(read.table(paste("data/TEST",id,".txt",sep=""))[,1])
fitx = normalmixEM(x,k=3)
plot(fitx,which=2)
#lines(density(x), lty=2, lwd=2)
fitx$loglik

mus<-c()
sigmas<-c()
lambdas<-c()
best_fits <- list()
# pdf("mixN.pdf")
for(id in ids){
  x<-sort(read.table(paste("data/TEST",id,".txt",sep=""))[,1])
  bestll<--10000000
  for(int in 1:10){
    fitx<-normalmixEM(x,k=3)
    if(fitx$loglik>bestll){
      bestll<-fitx$loglik
      best<-fitx
    }
  }
  best_fits[[id+1]] <- best
  # par(mfrow=c(2,1))
  # plot(best,which=2)
  # z<-ecdf(x)
  # F<-function(x,lambda,mu,sigma){Fs<-0;for (i in 1:3) Fs<-Fs+lambda[i]*pnorm(x,mean=mu[i],sd=sigma[i]);Fs}
  # ys<-sapply(x,FUN=F,lambda=best$lambda,mu=best$mu,sigma=best$sigma)
  # plot(x,z(x),type="l")
  # lines(x,ys,col="red")
  # mus<-c(mus,best$mu)
  # sigmas<-c(sigmas,best$sigma)
  # lambdas<-c(lambdas,best$lambda)
}
saveRDS(best, file = "best_fit_mixture_data.RDS")
# dev.off()
#veiem el plot que només funcionava en alguns, p.e 0-3,7-9, pero depenia de la iteracio, s'ha de fer el fit millor, iterant i maximitzant loglik (fet)
#sembla clarament degut a la asimetria, provar amb gamma