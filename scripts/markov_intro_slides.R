xs<-seq(1,170,1)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,
        0.6*dnorm(xi,mean=5,sd=10)+
          0.3*dnorm(xi,mean=50,sd=10)+
          0.1*dnorm(xi,mean=100,sd=10)
  )
}
plot(xs,ys)

xs<-seq(1,170,1)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,1-(
    0.6*pnorm(xi,mean=5,sd=10)+
      0.3*pnorm(xi,mean=50,sd=10)+
      0.1*pnorm(xi,mean=100,sd=10)
  ))
}
cota<-(0.6*5+0.3*50+0.1*100)/xs
plot(xs,ys,type="l",log="y")
lines(xs,cota)
for (k in 1:100)
{
  mostra<-c(rnorm(6000000,mean=5,sd=10),rnorm(3000000,mean=50,sd=10),rnorm(1000000,mean=100,sd=10))
  cotak<-mean(mostra^k)/xs^k
  #cotak<-t.test(mostra^k)$conf.int[2]/xs^k
  #cotak<-(mean(mostra^k)+3*sd(mostra^k))/xs^k
  lines(xs,cotak,col="blue")
}
lines(xs,ys,col="red")

###########################
result<-c()
for (n in c(100,1000))
{
  level<-0.9999
  sigmas<-5
  risc<-1e-12
  for (sim in 1:2000)
  {
    mostra<-c(rnorm(0.6*n,mean=5,sd=10),rnorm(0.3*n,mean=50,sd=10),rnorm(0.1*n,mean=100,sd=10))
    taula<-c()
    for (i in seq(1,2,0.01))
    {
      a<-i*max(mostra)
      xi<-a
      true<-0.6*pnorm(xi,mean=5,sd=10)+0.3*pnorm(xi,mean=50,sd=10)+0.1*pnorm(xi,mean=100,sd=10)
      taula<-rbind(taula,c(a,1-true,i))
    }
    a<-taula[(abs(taula[,2]-risc))==min(abs(taula[,2]-risc)),1]
    i<-taula[(abs(taula[,2]-risc))==min(abs(taula[,2]-risc)),3]
    a
    i
    xi<-a
    true<-0.6*pnorm(xi,mean=5,sd=10)+0.3*pnorm(xi,mean=50,sd=10)+0.1*pnorm(xi,mean=100,sd=10)
    1-true
    mostra0<-mostra/a
    mean(mostra0)
    mm<-c();for (k in 1:140) mm<-c(mm,mean(mostra0^k))
    mm1<-c();for (k in 1:140) mm1<-c(mm1,mean(mostra0^k)+sigmas*sd(mostra0^k))#t.test(mostra0^k,conf.level=level)$conf.int[1])
    mm2<-c();for (k in 1:140) mm2<-c(mm2,mean(mostra0^k)-sigmas*sd(mostra0^k))#t.test(mostra0^k,conf.level=level)$conf.int[2])
    mm[mm<0]<-0
    mm1[mm1<0]<-0
    mm2[mm2<0]<-0
    plot(1:140,mm,ylim=c(min(mm,mm1,mm2)+1e-40,max(mm,mm1,mm2)),log="y",type="l")
    lines(1:140,mm1,lty=2)
    lines(1:140,mm2,lty=2)
    abline(h=1-true)
    #abline(v=max((1:140)[mm-risc>0])
    result<-rbind(result,c(n,i,a,mean(mostra),max(mostra),max((1:140)[mm-risc>0],na.rm=T)))
  }
}
plot(result[,4],result[,5])
plot(result[,4],result[,6])
plot(result[,5],result[,6])
result[,1]==1000
cvplot(-result[,3])
################################
# o sigui que amb k=40 per exemple, no infraestima per n=100
#### a veure
a<-167.06
true<-0.6*pnorm(xi,mean=5,sd=10)+0.3*pnorm(xi,mean=50,sd=10)+0.1*pnorm(xi,mean=100,sd=10)
risc<-1-true
n<-1000
k<-50
mostra<-c(rnorm(0.6*n,mean=5,sd=10),rnorm(0.3*n,mean=50,sd=10),rnorm(0.1*n,mean=100,sd=10))
as<-seq(max(mostra),2*max(mostra),0.01)
taula<-c()
for (ai in as)
{
  mostra0<-mostra/ai
  taula<-rbind(taula,c(ai,mean(mostra0^k)))
}
taula[abs(taula[,2]-risc)==min(abs(taula[,2]-risc)),1]
#194.9
thr<-thrselect(mostra,evi=0)$threshold
z<-fitpot(mostra,threshold=thr,evi=0)
ccdfplot(mostra,z)
qpot(1-1e-12,z$coeff)
#279
z<-fitpot(mostra,threshold=thr)
ccdfplot(mostra,z)
qpot(1-1e-12,z$coeff)
#141.4
###################################

n<-1000
mostra<-c(rnorm(0.6*n,mean=5,sd=10),rnorm(0.3*n,mean=50,sd=10),runif(0.1*n,100,130))
taula<-c()
for (i in seq(1,2,0.01))
{
  a<-i*max(mostra)
  xi<-a
  true<-0.6*pnorm(xi,mean=5,sd=10)+0.3*pnorm(xi,mean=50,sd=10)+0.1*punif(xi,100,130)
  taula<-rbind(taula,c(a,1-true,i))
}
a<-taula[(abs(taula[,2]-risc))==min(abs(taula[,2]-risc)),1]
i<-taula[(abs(taula[,2]-risc))==min(abs(taula[,2]-risc)),3]
a
i
xi<-a
true<-0.6*pnorm(xi,mean=5,sd=10)+0.3*pnorm(xi,mean=50,sd=10)+0.1*punif(xi,100,130)
1-true
mostra0<-mostra/a
mean(mostra0)
mm<-c();for (k in 1:140) mm<-c(mm,mean(mostra0^k))
mm1<-c();for (k in 1:140) mm1<-c(mm1,mean(mostra0^k)+sigmas*sd(mostra0^k))#t.test(mostra0^k,conf.level=level)$conf.int[1])
mm2<-c();for (k in 1:140) mm2<-c(mm2,mean(mostra0^k)-sigmas*sd(mostra0^k))#t.test(mostra0^k,conf.level=level)$conf.int[2])
mm[mm<0]<-0
mm1[mm1<0]<-0
mm2[mm2<0]<-0
plot(1:140,mm,ylim=c(min(mm,mm1,mm2)+1e-40,max(mm,mm1,mm2)),log="y",type="l")
lines(1:140,mm1,lty=2)
lines(1:140,mm2,lty=2)
abline(h=1-true)
#abline(v=max((1:140)[mm-risc>0])
result<-rbind(result,c(n,i,a,mean(mostra),max(mostra),max((1:140)[mm-risc>0],na.rm=T)))


##########################3
plot(xs,ys,type="l",log="y",col="grey")
n<-1000
level<-0.95
mostra<-c(rnorm(0.6*n,mean=5,sd=10),rnorm(0.3*n,mean=50,sd=10),rnorm(0.1*n,mean=100,sd=10))
cotak0s<-c()
cotak1s<-c()
cotak2s<-c()
for (k in 1:140)
{
  #cotak<-mean(mostra^k)/xs^k
  x0<-c();x1<-c();x2<-c()
  for (xi in xs)
  {
    mostra0<-mostra/xi
    x0<-c(x0,mean(mostra0^k))
    x1<-c(x1,t.test(mostra0^k,conf.level=level)$conf.int[1])
    x2<-c(x2,t.test(mostra0^k,conf.level=level)$conf.int[2])
  }
  x1[x1<0]<-NA;x1[x1>1]<-NA;x2[x2<0]<-NA;x2[x2>1]<-NA
  cotak0s<-cbind(cotak0s,x0)
  cotak1s<-cbind(cotak1s,x1)
  cotak2s<-cbind(cotak2s,x2)
}
lines(xs,apply(cotak2s,1,min,na.rm=T),col="red")
lines(xs,apply(cotak1s,1,min,na.rm=T))
lines(xs,sort(mostra),col="orange")
abline(h=1/n,lty=2)

image(cotak2s)
contour(xs,1:140,log(cotak2s[,1:140],10),log="y",levels=-10:-1)






cotak1s<-c()
cotak2s<-c()
for (k in 1:150)
{
  #cotak<-mean(mostra^k)/xs^k
  x<-t.test(mostra^k,conf.level=0.8)$conf.int[1]/xs^k
  x[x<0]<-NA
  x[x>1]<-NA
  cotak1s<-cbind(cotak1s,x)
  x<-t.test(mostra^k,conf.level=0.8)$conf.int[2]/xs^k
  x[x<0]<-NA
  x[x>1]<-NA
  cotak2s<-cbind(cotak2s,x)
  #cotak<-(mean(mostra^k)+3*sd(mostra^k))/xs^k
  #lines(xs,cotak,col="blue")
}
lines(xs,apply(cotak2s,1,min,na.rm=T),col="orange")
#lines(xs,apply(cotak1s,1,min,na.rm=T),col="blue")
abline(h=1/n,lty=2)
z<-ecdf(mostra)
lines(xs,1-z(xs))

########################


lines(xs,apply(cotak1s,1,max,na.rm=T),col="red")


for (k in 1:300)
  lines(xs,cotak2s[,k],col="red")



par(mfrow=c(1,2))
xs<-seq(0,170,0.01)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,
        0.6*dnorm(xi,mean=5,sd=10)+
          0.3*dnorm(xi,mean=50,sd=10)+
          0.1*dnorm(xi,mean=100,sd=10)
  )
}
plot(xs,ys)

par(mfrow=c(1,2))

xs<-seq(1,170,0.01)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,1-(
    0.6*pnorm(xi,mean=5,sd=10)+
      0.3*pnorm(xi,mean=50,sd=10)+
      0.1*pnorm(xi,mean=100,sd=10)
  ))
}
plot(xs,ys,type="l",log="y",ylim=c(1e-14,1))
cota<-(0.6*5+0.3*50+0.1*100)/xs
lines(xs,cota)

for (k in 1:300)
{
  mostra<-c(rnorm(600000,mean=5,sd=10),rnorm(300000,mean=50,sd=10),rnorm(100000,mean=100,sd=10))
  cotak<-median(mostra^k)/xs^k
  #cotak<-t.test(mostra^k)$conf.int[2]/xs^k
  lines(xs,cotak,col="blue")
}

lines(xs,ys,col="red")





mostra<-c(rnorm(60000,mean=5,sd=10),rnorm(30000,mean=50,sd=10),rnorm(10000,mean=100,sd=10))
z<-ecdf(mostra)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,1-z(xi))
}
plot(xs,ys,type="l",log="y",col="blue",ylim=c(1e-14,1))
cota<-mean(mostra)/xs
lines(xs,cota,col="blue")

k<-10
for (xi in xs)
{
  cotak<-mean((mostra/xi)^k)
  cotak1<-t.test((mostra/xi)^k)$conf.int[1]
  cotak2<-t.test((mostra/xi)^k)$conf.int[2]
  segments(xi,cotak1,xi,cotak2,col="blue")
}



lines(xs,ys,col="red")


for (k in c(1,10,100))#1:300)
{
  cotak<-mean(mostra^k)/xs^k
  cotak1<-t.test(mostra^k)$conf.int[1]/xs^k
  cotak2<-t.test(mostra^k)$conf.int[2]/xs^k
  segments(xs,cotak1,xs,cotak2,col="blue")
}
lines(xs,ys,col="red")



#########################################
require(lmom)

xs<-seq(0,150,0.01)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,
        0.6*dnorm(xi,mean=5,sd=10)+
          0.3*dnorm(xi,mean=50,sd=10)+
          0.1*dnorm(xi,mean=100,sd=10)
  )
}
plot(xs,ys)

xs<-seq(0,170,0.01)
ys<-c()
for (xi in xs)
{
  ys<-c(ys,1-(
    0.6*pnorm(xi,mean=5,sd=10)+
      0.3*pnorm(xi,mean=50,sd=10)+
      0.1*pnorm(xi,mean=100,sd=10)
  ))
}
EX<-(0.6*5+0.3*50+0.1*100)
EX2<-(0.6*(5^2+10^2)+0.3*(50^2+10^2)+0.1*(100^2+10^2))
EXk_normal<-function(kk,mu,sigma)
{
  k<-trunc(kk/2)
  if (kk-2*k==0) (VX^kk)*factorial(2*k)/((2^k)*factorial(k))+EX^kk
  else EX^kk
}
EXk<-function(k) (0.6*EXk_normal(k,5,10)+0.3*EXk_normal(k,50,10)+0.1*EXk_normal(k,100,10))

VX<-EX2-EX^2
plot(xs,ys,type="l",log="y")
lines(xs,cota)
for (k in 1:300)
{
  cotak_teo<-EXk(k)/xs^k
  mostra<-c(rnorm(60000,mean=5,sd=10),rnorm(30000,mean=50,sd=10),rnorm(10000,mean=100,sd=10))
  #cotak<-mean(mostra^k)/xs^k
  cotak<-t.test(mostra^k)$conf.int[2]/xs^k
  #cotak<-(mean(mostra^k)+3*sd(mostra^k))/xs^k
  lines(xs,cotak_teo,col="blue")
}
lines(xs,ys,col="red")

lmrnor(c(2,3), nmom=4)



