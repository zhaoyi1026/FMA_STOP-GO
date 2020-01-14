######################################################
# Fourier basis function
fourier.basis<-function(timeinv=c(0,1),ntp,nbasis=3)
{
  if(nbasis%%2==0)
  {
    nbasis<-nbasis-1
  }
  
  timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  
  r<-timeinv[2]-timeinv[1]
  
  basis<-matrix(NA,ntp,nbasis)
  basis[,1]<-rep(1/sqrt(r),ntp)
  for(j in 1:floor(nbasis/2))
  {
    basis[,j*2]<-sin(j*(2*pi*(timegrids-timeinv[1])/r))*sqrt(2/r)
    basis[,j*2+1]<-cos(j*(2*pi*(timegrids-timeinv[1])/r))*sqrt(2/r)
  }
  
  return(basis)
}
# second order derivative of Fourier basis fuction
Ld2.fourier<-function(timeinv=c(0,1),ntp,nbasis=3)
{
  if(nbasis%%2==0)
  {
    nbasis<-nbasis-1
  }
  
  timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  
  r<-timeinv[2]-timeinv[1]
  
  db<-matrix(NA,ntp,nbasis)
  db[,1]<-rep(0,ntp)
  for(j in 1:floor(nbasis/2))
  {
    db[,j*2]<--sqrt(2/r)*(2*pi*j/r)^2*sin(2*pi*j*(timegrids-timeinv[1])/r)
    db[,j*2+1]<--sqrt(2/r)*(2*pi*j/r)^2*cos(2*pi*j*(timegrids-timeinv[1])/r)
  }
  
  return(db)
}
######################################################

######################################################
# intergral of a function
int.func<-function(x,timeinv=c(0,1),timegrids=NULL)
{
  ntp<-length(x)
  
  if(is.null(timegrids))
  {
    timegrids<-seq(timeinv[1],timeinv[2],length.out=ntp)
  }
  
  if(timeinv[2]<=timeinv[1])
  {
    int<-0
  }else
  {
    if(length(timegrids)==ntp)
    {
      int<-0
      for(i in 1:(ntp-1))
      {
        int<-int+(x[i]+x[i+1])*(timegrids[i+1]-timegrids[i])/2
      }
    }else
    {
      stop("Error!")
    }
  }
  
  return(int)
}
######################################################
