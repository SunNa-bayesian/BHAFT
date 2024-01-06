

library(brms)
library(rstan)
library(BhGLM)
library(survival)
library(survminer)
library(Hmisc)

set.seed(23456)
beta<-rnorm(11,0.7,0.05)
deta<-c(-1,1,1,-1,-1,
        1,-1,1,
        -1,1,1)
beta<-beta*deta  
beta<-c(beta,0.5)
nz=c(1,2,3,4,5,201,202,205,206,207,208)

beta




testsimsurv<- function(n=1, size=500, nsnp=200, nclin=5, corr=0.5, 
                       sigma,nz=nz,cr=0.3){
  simdatasum<-data.frame()
  for(j in 1:n){
    X=sim.x(n=size,m=nsnp,corr=corr)
    clin=sim.x(n=size,m=nclin,corr = corr)
    names(clin)<-c("x10001","x10002","x10003","x10004","x10005")
    var=cbind(X,clin)
    
    for(i in 1:nsnp){
      x10005x<-var[,i]*clin[,5]
      var<-cbind(var,x10005x)
    }
    aa=nsnp+nclin+1
    bb=nsnp+nclin+nsnp
    names(var)[aa:bb]<-as.list(paste("x10005x", c(1:nsnp), sep = ""))
    
    nn<-rep(j, times=size)
    ttime=stats::runif(size,min=0,max=1)
    epsilon=log(-log((1-ttime),exp(1)),exp(1))
    
    sigvar<-var[,nz]
    sigvarr<-cbind(sigvar,epsilon)
    sigvarr<-as.matrix(sigvarr)
    timeadj <- exp( sigvarr %*% beta )
    
    status=stats::rbinom(size,1,1-cr)
    clidata <- data.frame(nn,timeadj,status)
    
    clidata <- within(clidata,{
      indicate <- NA
      indicate[status==1] <- 1
      indicate[status==0] <- runif(size-sum(clidata$status), min=0, max=1)})
    clidata$time<-clidata$timeadj*clidata$indicate
    clidata <- clidata[,-c(2,4)]
    total<-cbind(clidata,var)
    simdatasum<-rbind(simdatasum,total)
  }
  return(simdatasum)
}


set.seed(12345)
data<-testsimsurv(n=100, size=500, nsnp=200, nclin=5, 
                  corr=0.5, sigma=1,
                  nz=nz,cr=0.3)

min=min(data$time)
max=max(data$time)
table(data$status)




set.seed(1234)
train<-testsimsurv(n=100, size=100, nsnp=200, nclin=5, 
                   corr=0.5, sigma=1,
                   nz=nz,cr=0.3)






setwd("C:\\Users\\Dell\\Desktop\\horseshoe\\")
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  horseshoe=set_prior("horseshoe(df=1,df_global=1)")
  #set.seed(20220902)
  f <- brm(time|cens(1-status) ~ .,
           data = datanew, 
           prior=horseshoe,
           chain=3,iter=3000, 
           #control = list(adapt_delta = 0.99,max_treedepth=12),
           family = weibull())
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  pre<-predict(f,newdata = trainnew)
  C_index=rcorr.cens(pre[,1],Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}












setwd("C:\\Users\\Dell\\Desktop\\ssehhorseshoe\\")
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  nsnp=200
  nvar=205
  nz=205
  stanvars=stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
    stanvar(x=nz,name="nz",scode="int nz;",block="data")
  horseshoe=set_prior("horseshoe(df=1,df_global=1)")
  
  sdata<-make_standata(time|cens(1-status) ~ .,
                       data = datanew,
                       family = weibull(),
                       prior=horseshoe,stanvars=stanvars)
  #set.seed(20220902)
  kf=rstan::stan(file ="C:\\Users\\Dell\\Desktop\\sshorseshoe_weibull20220706_intercept.stan",
                 data=sdata,
                 #control = list(adapt_delta = 0.99,max_treedepth=12),
                 chain=3,iter=3000)
  
  # feed the Stan model back into brms
  f <- brm(time|cens(1-status) ~ ., 
           data = datanew, 
           prior=horseshoe,stanvars=stanvars,
           family = weibull(),
           empty = TRUE)
  f$fit <- kf
  f <- rename_pars(f)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  pre<-predict(f,newdata = trainnew)
  C_index=rcorr.cens(pre[,1],Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}






setwd("C:\\Users\\Dell\\Desktop\\zzssde_hp2\\")
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  nvar=405
  s0=0.05
  ssde=set_prior("for (j in 1:nvar)
               target += log_sum_exp(log(1-gamma[j])+normal_lpdf(b[j]|0,s0*sqrttau),
                                     log(gamma[j])+normal_lpdf(b[j]|0,sqrttau))",check=F)+
    set_prior("target += inv_gamma_lpdf(tau | 3, 12)",check=F)+
    set_prior("for (j in 1:nvar)
               target += beta_lpdf(gamma[j]|1,1)",check=F)
  stanvars=stanvar(scode="real<lower=0> tau;",block="parameters")+
    stanvar(scode="real<lower=0> sqrttau=sqrt(tau);",block="tparameters")+
    stanvar(x=s0,name="s0",scode="real s0;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
    stanvar(scode="vector<lower=0,upper=1>[nvar] gamma;",block="parameters")
  
  #set.seed(20220902)
  f=brm(time|cens(1-status) ~ .,
        data = datanew, 
        family = weibull(),
        prior=ssde,stanvars=stanvars,
        chain=3,iter=5000)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  pre<-predict(f,newdata = trainnew)
  C_index=rcorr.cens(pre[,1],Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}






setwd("C:\\Users\\Dell\\Desktop\\zzhheredity_hp2\\")
n=100
for(j in 1:n){
  datanew<-data[which(data$nn==j),]
  datanew<-datanew[,-1]
  trainnew<-train[which(train$nn==j),]
  trainnew<-trainnew[,-1]
  
  fpath<-".\\f\\"
  sigvarpath<-".\\sigvar\\"
  Cindexpath<-".\\Cindex\\"
  
  nsnp=200 
  nvar=205
  s0=0.05
  
  eefssde=set_prior("for (j in 1:nsnp)
               target += log_sum_exp(log(1-gamma[j])+normal_lpdf(b[j]|0,s0*sqrttau),
                                     log(gamma[j])+normal_lpdf(b[j]|0,sqrttau))",check=F)+
    set_prior("for (j in (nsnp+1):nvar)
               target += normal_lpdf(b[j]|0,sqrttau)",check=F)+
    set_prior("for (j in (nvar+1):(nvar+nsnp))
               target += log_sum_exp(log(1-gamma[j-nvar])+normal_lpdf(b[j]|0,s0*sqrttau),
                                     log(gamma[j-nvar])+normal_lpdf(b[j]|0,sqrttau))",check=F)+
    set_prior("target += inv_gamma_lpdf(tau | 3, 12)",check=F)+
    set_prior("for (j in 1:nsnp)
               target += beta_lpdf(gamma[j]|1,1)",check=F)
  
  eefstanvars=stanvar(scode="vector<lower=0,upper=1>[nsnp] gamma;",block="parameters")+
    stanvar(scode="real<lower=0> tau;",block="parameters")+
    stanvar(scode="real<lower=0> sqrttau=sqrt(tau);",block="tparameters")+
    stanvar(x=s0,name="s0",scode="real s0;",block="data")+
    stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
    stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")
  
  #set.seed(20220902)
  f=brm(time|cens(1-status) ~ .,
        data = datanew, 
        family = weibull(),
        prior=eefssde,stanvars=eefstanvars,
        chain=3,iter=5000)
  
  result<-summary(f)
  resultt<-result[["fixed"]]
  sigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]
  
  pre<-predict(f,newdata = trainnew)
  C_index=rcorr.cens(pre[,1],Surv(trainnew$time,trainnew$status))
  
  save(f,file=paste(fpath,"f", j,".RData", sep = ""))
  write.csv(sigvar,paste(sigvarpath,"sigvar", j,".csv", sep = ""))
  write.csv(C_index,paste(Cindexpath,"Cindex", j,".csv", sep = ""))
  print(j)
}




