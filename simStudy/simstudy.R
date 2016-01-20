### simulation study
library(geoR)
library(foreach)
library(doParallel)
library(asreml)
library(latticeExtra)
library(xtable)

Nlist=c(300,600,1000)  ## This determine the sample size. Large value = longer wait. For computer with <5 cores, consider a maximum of 600.
xlim=c(0,1000)
ylim=c(0,1000)
Nsim=500  ## The number of replication. Large value = longer wait. For computer with <5 cores, consider a max of 20.
covpars=c(0.5,200)
nugget=0.1
mu=5

Nblock=9
centers=expand.grid(seq(from=xlim[1],to=xlim[2],length=7)[c(2,4,6)],seq(from=xlim[1],to=xlim[2],length=7)[c(2,4,6)])
## to visualize the centers
plot(centers,type="n")
text(centers, labels=1:9)
pointsklpairing=list(c(1,2),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(2,6),c(3,5),c(3,6),c(4,5),c(4,7),c(4,8),c(5,6),c(5,7),c(5,8),c(5,9),c(6,8),c(6,9),c(7,8),c(8,9))
pointsjpair=plyr::alply(combn(20,2),.margins = 2)


cl=makeCluster(20)
registerDoParallel(cl)
t1=proc.time()
out=foreach(N=Nlist,.packages=c("geoR","Matrix","foreach","asreml"),.export=c("dist"),.errorhandling = "pass")%:%
  foreach(i=1:Nsim,.combine=rbind)%dopar%{
    simdata=geoR::grf(n=N,xlims = xlim,ylims = ylim,nsim = 1,cov.model = "exp",
                      cov.pars=covpars,nugget=nugget,mean=mu)
    v=variofit(variog(simdata,max.dist=500),ini.cov.pars = covpars*1.1,nugget=nugget*1.5,fix.nugget=FALSE)
    
    ## geoR
    gRt=system.time({
      gRfit=likfit(simdata,ini.cov.pars = v,fix.nugget = FALSE,nugget = v$nugget,lik.method="REML")
    })
    gRstat=data.frame(i=i,"usrTime"=gRt[1],"sysTime"=gRt[2],"elpsdTime"=gRt[3],"sigma2"=gRfit[[3]][1],"phi"=gRfit[[3]][2],"nugget"=gRfit$nugget,"beta"=gRfit$beta,method="geoR")
    
    
    df=as.data.frame(simdata)
    #   km=kmeans(df[,1:2],centers =centers,nstart = 10)
    #   pk=lapply(1:Nblock,function(x)which(km$cluster==x))
    #   pj=lapply(pointsklpairing,function(x)unlist(pk[x]))
    
    
    ## asreml
    asr.sv=asreml(data~1,random=~units,rcov = ~ieuc(x,y),data=df,start.values=TRUE)$gammas.table
    asr.sv$Value=c(nugget/covpars[1]*1.1,1,exp(-1/covpars[2]/1.1))
    asrT=system.time({asr=asreml(data~1,random=~units,rcov = ~ieuc(x,y),data=df,G.param=asr.sv,R.param=asr.sv)})
    asrcovpar=summary(asr)$varcomp$component
    asrstat=data.frame(i=i,"usrTime"=asrT[1],"sysTime"=asrT[2],"elpsdTime"=asrT[3],"sigma2"=asrcovpar[2]-asrcovpar[1],"phi"=-1/log(asrcovpar[3]),"nugget"=asrcovpar[1],"beta"=c(coef(asr)$fixed),method="asreml")
    #   
    ##CompReLik
    #   clsmfit=fitCLSM(form = data~1,data = df,dS = dSexp,coordsMat = df[,1:2],init = c(v$cov.pars,v$nugget),
    #           pointsj = pj,pointsjpair = pointsjpair,pointsk = NULL,
    #           cov.model = "exp",heter = FALSE,lower = c(0,1,0))
    #   clsmstat=data.frame(i=i,"usrTime"=clsmfit$time[1],"sysTime"=clsmfit$time[2],"sigma2"=clsmfit$par[1],"phi"=clsmfit$par[2],"nugget"=clsmfit$par[3],"beta"=clsmfit$betahat,method="clsm")
    
    #   out0=rbind(gRstat,asrstat,clsmstat)
    out0=rbind(gRstat,asrstat)
    rownames(out0)=NULL
    out0$N=N
    out0
  }
proc.time()-t1

# stopCluster(cl)

save(out,file="simresult4.rda")




### C-REML
t3=proc.time()
outclsm=foreach(N=Nlist,.packages=c("geoR","Matrix","foreach","asreml"),.export=c("dist"),.errorhandling = "pass")%:%
  foreach(i=1:Nsim,.combine=rbind)%do%{
    simdata=geoR::grf(n=N,xlims = xlim,ylims = ylim,nsim = 1,cov.model = "exp",
                      cov.pars=covpars,nugget=nugget,mean=mu)
    v=variofit(variog(simdata,max.dist=500),ini.cov.pars = covpars*1.1,nugget=nugget*1.5,fix.nugget=FALSE)
    
    
    df=as.data.frame(simdata)
    km=kmeans(df[,1:2],centers =centers,nstart = 10)
    pk=lapply(1:Nblock,function(x)which(km$cluster==x))
    pj=lapply(pointsklpairing,function(x)unlist(pk[x]))
    
    
    ##CompReLik
    clsmfit=fitCLSM(form = data~1,data = df,dS = dSexp,coordsMat = df[,1:2],init = c(v$cov.pars,v$nugget),
                    pointsj = pj,pointsjpair = pointsjpair,pointsk = NULL,
                    cov.model = "exp",heter = FALSE,lower = c(0,1,0))
    clsmstat=data.frame(i=i,"usrTime"=clsmfit$time[1],"sysTime"=clsmfit$time[2],"elpsdTime"=clsmfit$time[3],"sigma2"=clsmfit$par[1],"phi"=clsmfit$par[2],"nugget"=clsmfit$par[3],"beta"=clsmfit$betahat,method="clsm",N=N)
    
    out0=rbind(clsmstat)
    rownames(out0)=NULL
    out0
  }
proc.time()-t3

save(outclsm,file="simCLSM4.rda")

simout=rbind(Reduce(rbind,out),Reduce(rbind,outclsm))
## summarize results ------------------
## format data
simout$N=as.factor(simout$N)
bwplot(elpsdTime~method|N,data=simout,horizontal = FALSE)

xyplot(elpsdTime+usrTime~N,groups=method,
       data=aggregate(cbind(usrTime,elpsdTime)~method+N,data=simout,FUN=mean),
       ylab=list("User Time (s)",cex=2),xlab=list("N",cex=2),lwd=2,pch=1:3,
       type="b",grid=TRUE,cex=2,scales=list(x=list(cex=2),y=list(cex=2)),
       key=list(columns=3,text=list(c("REML (geoR)","AsReml","C-REML"),cex=2),
                points=list(pch=1:3,col=trellis.par.get("superpose.symbol")$col[1:3],cex=2),
                lines=list(lty=1:3,col=trellis.par.get("superpose.symbol")$col[1:3],cex=2,lwd=2)))


medianpar=aggregate(cbind(usrTime,elpsdTime,sigma2,phi,beta,nugget)~method+N,data=simout,FUN = median)
# meanpar=aggregate(cbind(usrTime,elpsdTime,sigma2,phi,beta,nugget)~method+N,data=simout,FUN = mean)
# aggregate(cbind(usrTime,sysTime,sigma2,phi,beta,nugget)~method+N,data=simout,FUN = quantile, probs=c(0.05,0.5,0.95))



xyplot(sigma2+phi+beta+nugget~N,groups=method,data=medianpar,grid=TRUE,
       ylab="",xlab=list("N",cex=2),par.strip.text=list(cex=2),pch=1:3,
       scales=list(x=list(cex=1.5,relation="free"),y=list(cex=1.5,relation="free",rot=0)),
       cex=2,panel=function(x,y,...){
         panel.xyplot(x,y,...)
         panel.abline(h=c(5,0.1,0.5,200),lty=2)
       },
       key=list(columns=3,text=list(c("REML (geoR)","AsReml","C-REML"),cex=2),
                points=list(pch=1:3,col=trellis.par.get("superpose.symbol")$col[1:3],cex=2)
       ))
