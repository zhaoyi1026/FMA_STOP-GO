#################################################
# openfMRI ds030
# STOP/GO task
# Functional mediation

# STOP trial
# Bootstrap result summary
#################################################

library("cfma")
library("plot3D")

rm(list=ls())

# args<-commandArgs(trailingOnly=TRUE)
# Figure 6
args<-c("180319","180319","preSMA-post","M1-RH",11,11,500,10,3,2)
# Figure 7
args<-c("180319","180319","preSMA-ant","M1-RH",11,11,500,4,15,3)

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

setwd(paste0("/Users/yizhao/Dropbox/MyFolder/Biostat-JHU/2017/Mediation/functional_mediation/R/",run.date,"/Rpkg"))

dir.func<-getwd()
source(paste0(dir.func,"/functions.R"))
##############################################################################
# load data
dir.data<-paste0("/Users/yizhao/Dropbox/MyFolder/Biostat-JHU/2017/Mediation/functional_mediation/R/",run.date,"/Data")

# read Z
env.tmp<-new.env()
load(paste0(dir.data,"/Data_Z.RData"),env.tmp)
Zraw<-env.tmp$Z
Z0<-env.tmp$Z0
Z1<-env.tmp$Z1
TR<-env.tmp$TR
rm("env.tmp")

# read M
env.tmp<-new.env()
load(paste0(dir.data,"/Data_",M.label,".RData"),env.tmp)
Mraw<-env.tmp$Y
rm("env.tmp")

# read R
env.tmp<-new.env()
load(paste0(dir.data,"/Data_",R.label,".RData"),env.tmp)
Yraw<-env.tmp$Y
rm("env.tmp")
##############################################################################

##############################################################################
sub.ID<-rownames(Yraw)
N<-nrow(Yraw)
ntp<-ncol(Zraw)

# Remove GO trial
M=Y<-matrix(NA,N,ntp)
rownames(M)=rownames(Y)<-sub.ID
for(i in 1:N)
{
  fit.m<-lm(Mraw[i,]~Z0[i,]+Z1[i,])
  M[i,]<-scale(Mraw[i,]-cbind(rep(1,ntp),Z0[i,])%*%fit.m$coefficients[1:2],center=TRUE,scale=TRUE)
  
  fit.y<-lm(Yraw[i,]~Z0[i,]+Z1[i,])
  Y[i,]<-scale(Yraw[i,]-cbind(rep(1,ntp),Z0[i,])%*%fit.y$coefficients[1:2],center=TRUE,scale=TRUE)
}
Z<-Z1

timegrids<-((1:ntp)-1)*TR
timeinv<-c(timegrids[1],timegrids[ntp])
##############################################################################

##############################################################################
# onset
dir.onset<-"/Users/yizhao/Dropbox/MyFolder/Biostat-Brown/2013_Fall/RA/Mediantion/Data/openfMRI/ds030/STOPGO_task/170306/EV"

onset<-list()
for(i in 1:length(sub.ID))
{
  file.onset.tmp<-paste0(dir.onset,"/",sub.ID[i],"_onset.RData")
  
  env.onset.tmp<-new.env()
  load(file.onset.tmp,env.onset.tmp)
  
  onset[[i]]<-env.onset.tmp$ev.sort
}
##############################################################################

save(list=c("sub.ID","N","ntp","TR","Z","M","Y","timegrids","timeinv","onset"),file=paste0("Data_STOP_",M.label,"_",R.label,".RData"))
