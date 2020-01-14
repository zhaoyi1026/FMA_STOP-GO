#################################################
# openfMRI ds030
# STOP/GO task
# Functional mediation

# STOP trial
# Bootstrap result summary
#################################################

# R cfma (causal functional mediation analysis) package
# install.packages("cfma")
library("cfma")

rm(list=ls())

# args<-commandArgs(trailingOnly=TRUE)
# Figure 6
args<-c("180319","180319","preSMA-post","M1-RH",11,11,500,10,3,2)
# Figure 7
# args<-c("180319","180319","preSMA-ant","M1-RH",11,11,500,4,15,3)

run.date<-args[1]
data.date<-args[2]

M.label<-args[3]
R.label<-args[4]

nbasis1<-as.numeric(args[5])
nbasis2<-as.numeric(args[6])

sims<-as.numeric(args[7])

# time invervals for M and Y
delta.grid1<-as.numeric(args[8])
delta.grid2<-as.numeric(args[9])
delta.grid3<-as.numeric(args[10])

dir.func<-getwd()
source(paste0(dir.func,"/functions.R"))

set.seed(100000)
cc<-colors()[c(24,136,564,500,469,81,200,460,17,2,652,90,8,146,464,52,2,sample(c(1:152,190:260,300:657),500))]
pp<-c(10,8,19,17,15,18,11)
pp.cex<-c(0.9,1.25,1.5,1.25,1.25,1.25,0.75)
ll.wd<-c(1,3,3,2,2,2,1,2,2)+2
rho.lty<-c(3,1,5)

cex.lab<-1.25
cex.axis<-1.25
pt.cex<-1.25
lab.cex<-1.25
legend.cex<-1.5

dir.fig0<-paste0(getwd(),"/fig")
if(file.exists(dir.fig0)==FALSE)
{
  dir.create(dir.fig0)
}
dir.fig1<-paste0(dir.fig0,"/STOP")
if(file.exists(dir.fig1)==FALSE)
{
  dir.create(dir.fig1)
}
dir.fig2<-paste0(dir.fig1,"/",M.label,"_",R.label)
if(file.exists(dir.fig2)==FALSE)
{
  dir.create(dir.fig2)
}

##############################################################################
# load data
load(paste0("Data_STOP_",M.label,"_",R.label,".RData"))
##############################################################################

##############################################################################
timegrids<-((1:ntp)-1)*TR
timeinv<-c(timegrids[1],timegrids[ntp])

basis1<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
basis2<-fourier.basis(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)

Ld2.basis1<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis1)
Ld2.basis2<-Ld2.fourier(timeinv=timeinv,ntp=ntp,nbasis=nbasis2)

# method parameter
intercept<-FALSE
lambda1<-10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1])
lambda2<-10^c(seq(-2,1,length.out=20),seq(1,3,length.out=11)[-1])

# nfolds cross-validation to choose lambda
nfolds<-5

# sims<-500
boot<-TRUE
boot.ci.type<-"perc"
conf.level<-0.95
verbose<-TRUE
##############################################################################

##############################################################################
# K-fold cross-validation choosing lambda
re.CV.FMA<-FMA.historical.CV(Z,M,Y,delta.grid1=delta.grid1,delta.grid2=delta.grid2,
                             delta.grid3=delta.grid3,intercept=intercept,
                             basis1=basis1,Ld2.basis1=Ld2.basis1,
                             basis2=basis2,Ld2.basis2=Ld2.basis2,
                             timeinv=timeinv,timegrids=timegrids,
                             lambda1=lambda1,lambda2=lambda2,nfolds=nfolds)

# FMA output
re.FMA<-FMA.historical(Z,M,Y,delta.grid1=delta.grid1,delta.grid2=delta.grid2,
                       delta.grid3=delta.grid3,intercept=intercept,
                       basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,
                       timeinv=timeinv,timegrids=timegrids,
                       lambda1.m=re.CV.FMA$M$lambda1,lambda2.m=re.CV.FMA$M$lambda2,
                       lambda1.y=re.CV.FMA$Y$lambda1,lambda2.y=re.CV.FMA$Y$lambda2)

# Bootstrap
IEsub.est=DEsub.est<-array(NA,c(N,ntp,sims))
alpha.z.int.est<-array(NA,c(N,ntp,sims))

for(b in 1:sims)
{
  set.seed(as.numeric(paste0("20",run.date))+b*100)
  idx.tmp<-sample(1:N,N,replace=TRUE)
  
  Ztmp<-Z[idx.tmp,]
  Mtmp<-M[idx.tmp,]
  Ytmp<-Y[idx.tmp,]
  
  re<-NULL
  try(re<-FMA.historical(Ztmp,Mtmp,Ytmp,delta.grid1=delta.grid1,delta.grid2=delta.grid2,delta.grid3=delta.grid3,
                         intercept=intercept,basis1=basis1,Ld2.basis1=Ld2.basis1,basis2=basis2,Ld2.basis2=Ld2.basis2,timeinv=timeinv,timegrids=timegrids,
                         lambda1.m=re.CV.FMA$M$lambda1,lambda2.m=re.CV.FMA$M$lambda2,lambda1.y=re.CV.FMA$Y$lambda1,lambda2.y=re.CV.FMA$Y$lambda2))
  
  for(i in 1:N)
  {
    for(j in 1:ntp)
    {
      rtmp1<-max(j-delta.grid1,1)
      rtmp2<-max(j-delta.grid2,1)
      rtmp3<-max(j-delta.grid3,1)
      
      alpha.z.int.est[i,j,b]<-int.func(Z[i,rtmp1:j]*re$M$curve$alpha[rtmp1:j,j],timeinv=c(timegrids[rtmp1],timegrids[j]),timegrids=timegrids[rtmp1:j])
      DEsub.est[i,j,b]<-int.func(Z[i,rtmp2:j]*re$Y$curve$gamma[rtmp2:j,j],timeinv=c(timegrids[rtmp2],timegrids[j]),timegrids=timegrids[rtmp2:j])
      IEsub.est[i,j,b]<-int.func(alpha.z.int.est[i,rtmp3:j,b]*re$Y$curve$beta[rtmp3:j,j],timeinv=c(timegrids[rtmp3],timegrids[j]),timegrids=timegrids[rtmp3:j])
    }
  }
  
  print(paste0("Bootstrap sample: ",b))
}
##############################################################################

save(list=c("re.CV.FMA","re.FMA","IEsub.est","DEsub.est"),file=paste0("FMA_STOP_boot_",M.label,"_",R.label,"_MZ",delta.grid1,"_YZ",delta.grid2,"_YM",delta.grid3,".RData"))

##############################################################################
# plot and summary

load(paste0("FMA_STOP_boot_",M.label,"_",R.label,"_MZ",delta.grid1,"_YZ",delta.grid2,"_YM",delta.grid3,".RData"))

# Bootstrap results
IE.sub.boot=DE.sub.boot<-array(NA,c(N,ntp,3))
IE.sub.boot[,,1]<-apply(IEsub.est,c(1,2),mean,na.rm=TRUE)
DE.sub.boot[,,1]<-apply(DEsub.est,c(1,2),mean,na.rm=TRUE)

if(boot.ci.type=="bca")
{
  ci.tmp<-apply(IEsub.est,c(1,2),BC.CI,sims=dim(IEsub.est)[3],conf.level=conf.level)
  IE.sub.boot[,,2]<-ci.tmp[1,,]
  IE.sub.boot[,,3]<-ci.tmp[2,,]
  
  ci.tmp<-apply(DEsub.est,c(1,2),BC.CI,sims=dim(DEsub.est)[3],conf.level=conf.level)
  DE.sub.boot[,,2]<-ci.tmp[1,,]
  DE.sub.boot[,,3]<-ci.tmp[2,,]
}
if(boot.ci.type=="perc")
{
  ci.tmp<-apply(IEsub.est,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2))
  IE.sub.boot[,,2]<-ci.tmp[1,,]
  IE.sub.boot[,,3]<-ci.tmp[2,,]
  
  ci.tmp<-apply(DEsub.est,c(1,2),quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2))
  DE.sub.boot[,,2]<-ci.tmp[1,,]
  DE.sub.boot[,,3]<-ci.tmp[2,,]
}

# subject 94 in Figure 6/7
i<-94
IE.range.tmp<-range(c(range(IE.sub.boot[i,,],na.rm=TRUE)))
DE.range.tmp<-range(c(range(DE.sub.boot[i,,],na.rm=TRUE)))

png(paste0(dir.fig,"/sub",i,"_ID",sub.ID[i],".png"),width=20,height=10,unit="in",res=300)
par(mfrow=c(2,1))
onset.idx<-which(onset[[i]][,1]==2)
# IE
plot(range(timegrids),IE.range.tmp,type="n",xlab="time (s)",ylab="IE",cex.lab=cex.lab,cex.axis=cex.axis)
abline(h=0,lty=2,lwd=2,col=8)
matplot(timegrids,IEsub.est[i,,],type="l",lty=1,lwd=1,col=8,add=TRUE)
abline(v=onset[[i]][onset.idx,2],lty=c(3,1)[onset[[i]][onset.idx,1]],col=cc[c(6,2)][onset[[i]][onset.idx,1]],lwd=1)
matplot(timegrids,IE.sub.boot[i,,],type="l",lty=c(1,3,3),lwd=c(2,2,2),col=cc[1],add=TRUE)
# DE
plot(range(timegrids),DE.range.tmp,type="n",xlab="time (s)",ylab="DE",cex.lab=cex.lab,cex.axis=cex.axis)
abline(h=0,lty=2,lwd=2,col=8)
matplot(timegrids,DEsub.est[i,,],type="l",lty=1,lwd=1,col=8,add=TRUE)
abline(v=onset[[i]][onset.idx,2],lty=c(3,1)[onset[[i]][onset.idx,1]],col=cc[c(6,2)][onset[[i]][onset.idx,1]],lwd=1)
matplot(timegrids,DE.sub.boot[i,,],type="l",lty=c(1,3,3),lwd=c(2,2,2),col=cc[1],add=TRUE)
dev.off()
##############################################################################