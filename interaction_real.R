

library(brms)
library(rstan)
library(BhGLM)
library(survival)
library(survminer)
load("C:\\Users\\Dell\\Desktop\\NEWTCGALUAD.Rdata")
load("C:\\Users\\Dell\\Desktop\\mRNA_clinical_LUAD.Rdata")
TCGALUAD<-TCGALUAD[,c(1,12,13,15)]
table(TCGALUAD$treatment)
table(TCGALUAD$race)
table(TCGALUAD$tobacco_smoking_history)
TCGALUAD<-TCGALUAD[-which(TCGALUAD$tobacco_smoking_history==5),]
TCGALUAD$tobacco_smoking_history<-ifelse(TCGALUAD$tobacco_smoking_history==4,
                                         3,
                                         TCGALUAD$tobacco_smoking_history)
TCGALUAD$treatment<-ifelse(TCGALUAD$treatment=="Pharmaceutical Therapy, NOS",
                           1,2)
TCGALUAD <- within(TCGALUAD,{
  race[race == "white"] <- 1
  race[race == "american indian or alaska native" ] <- 2
  race[race == "asian"] <- 2
  race[race == "black or african american"] <- 2
  race[race == "not reported"] <- 3})

clinical<-mRNA_clinical
clinical<-clinical[,-c(3)]

table(clinical$OS)
table(clinical$stage)
clinical$gender<-ifelse(clinical$gender=="male",
                        1,
                        0)

clinical <- within(clinical,{
  stage[stage == "not reported"] <- NA
  stage[stage == "stage i" ] <- 1
  stage[stage == "stage ii"] <- 2
  stage[stage == "stage iii"] <- 3
  stage[stage == "stage iv"] <- 4})

TCGALUAD$sample_id<-substr(TCGALUAD$id, 1, 16)
clinical<-merge(TCGALUAD[,c(5,4,3,2)],clinical,by="sample_id")



#clinical<-na.omit(clinical)
clinical<-clinical[which(clinical$Time>0),]
#clinical<-clinical[which(clinical$age>40),]
#clinical<-clinical[which(clinical$age<81),]


clinical$age<-scale(clinical$age)

#for(i in 7:2060){
#  clinical[,i]<-scale(clinical[,i])
#}
clinical[,11:2064]<-covariates(clinical[,11:2064])

data<-clinical
aa<-as.matrix(data[,11:2064])
dim(aa)
data$stage<-as.numeric(data$stage)
for(i in 1:2054){
  x10005x<-aa[,i]*data$stage
  data<-cbind(data,x10005x)
}
names(data)[2065:4118]<-as.character(paste("s",names(data)[11:2064], sep = ""))





colnames(data)<-gsub("-","",colnames(data))
#**********循环语句，做单因素cox**********
tolresult <- data.frame(matrix(ncol = 5, nrow = 0))
xnam <- colnames(data[,c(11:2064)])

nnn<-length(xnam)
for (i in 1:nnn){
  fmla <- as.formula(paste("Surv(Time,OS) ~ ", paste(xnam[i])))
  f<-coxph(fmla, data)
  resultt=summary(f)[["coefficients"]]
  resultt=as.data.frame(resultt)
  sigvar = resultt
  
  tolresult <- rbind(tolresult,sigvar)
  print(i)
}



#**********循环语句，做单因素cox**********
itolresult <- data.frame(matrix(ncol = 5, nrow = 0))
xnaminter <- colnames(data[,c(2065:4118)])

nnn<-length(xnaminter)
for (i in 1:nnn){
  fmla <- as.formula(paste("Surv(Time,OS) ~ ", paste(xnaminter[i])))
  f<-coxph(fmla, data)
  resultt=summary(f)[["coefficients"]]
  resultt=as.data.frame(resultt)
  sigvar = resultt
  
  itolresult <- rbind(itolresult,sigvar)
  print(i)
}







itolresult<-itolresult[order(abs(itolresult[,4])),]
itop<-itolresult[which(itolresult$`Pr(>|z|)`<0.05/2054),]
tolresult<-tolresult[order(abs(tolresult[,4])),]
top<-tolresult[which(tolresult$`Pr(>|z|)`<0.05/2054),]
#top <- top[rownames(top) %in% substr(rownames(itop),2,nchar(rownames(itop))),]
top<-c(rownames(top),substr(rownames(itop),2,nchar(rownames(itop))))


#aaa<-tolresult[2005:2054,]
#top<-rbind(top,aaa)

#data$smoking<-scale(data$smoking)
G <- data[,colnames(data) %in% top]
inter <- data[,colnames(data) %in% as.character(paste("s",top, sep = ""))]
time<-data[,5:6]
E<-data[,c(2:3,7:9,4)]
dumrace<-model.matrix(~race,E)
dumrace<-dumrace[,-1]
data<-data.frame(time,G,E[,1:5],inter)
data<-na.omit(data)
data$tobacco_smoking_history<-as.factor(data$tobacco_smoking_history)
dumsmoke<-model.matrix(~tobacco_smoking_history,data)
dumsmoke<-dumsmoke[,-1]
data<-data.frame(data[,1:113],dumsmoke,data[,c(115:229)])


table(data$OS)
data$gender<-as.numeric(data$gender)
data$stage<-as.numeric(data$stage)
#data$OS<-1-data$OS
table(data$OS)
data<-na.omit(data)






horseshoe=set_prior("horseshoe(df=1,df_global=1)")
hf <- brm(Time|cens(1-OS) ~ .,
          data = data,
          prior=horseshoe,
          chain=3,iter=6000, 
          #control = list(adapt_delta = 0.99,max_treedepth=12),
          family = weibull())
result<-summary(hf)
resultt<-result[["fixed"]]
hsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hkfold <- kfold(hf, chains = 3, K=5,
                save_fits = TRUE)
hreskfold<-hkfold[["estimates"]]

hkfp <- kfold_predict(hkfold)
hyhat<-hkfp[["yrep"]]
hyrep_mean <- colMeans(hyhat)
hcindex=rcorr.cens(hyrep_mean,Surv(data$Time,data$OS))




nsnp=111
nvar=117
nz=117
stanvars=stanvar(x=nsnp,name="nsnp",scode="int nsnp;",block="data")+
  stanvar(x=nvar,name="nvar",scode="int nvar;",block="data")+
  stanvar(x=nz,name="nz",scode="int nz;",block="data")
horseshoe=set_prior("horseshoe(df=1,df_global=1)")

sdata<-make_standata(Time|cens(1-OS) ~ .,
                     data = data, 
                     family = weibull(),
                     prior=horseshoe,stanvars=stanvars)

kf=rstan::stan(file ="C:\\Users\\Dell\\Desktop\\sshorseshoe_weibull20220706_intercept.stan",
               data=sdata,
               #control = list(adapt_delta = 0.999,max_treedepth=12),
               #control = list(adapt_delta = 0.99),
               chain=3,iter=5000)

# feed the Stan model back into brms
fit <- brm(Time|cens(1-OS) ~ .,
           data = data, 
           prior=horseshoe,stanvars=stanvars,
           family = weibull(),
           empty = TRUE)
fit$fit <- kf
fit <- rename_pars(fit)

result<-summary(fit)
resultt<-result[["fixed"]]
ehsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

ehkfold <- kfold(fit, chains = 3, K=5,
                 save_fits = TRUE)
ehreskfold<-ehkfold[["estimates"]]

ehkfp <- kfold_predict(ehkfold)
ehyhat<-ehkfp[["yrep"]]
ehyrep_mean <- colMeans(ehyhat)
ehcindex=rcorr.cens(ehyrep_mean,Surv(data$Time,data$OS))




library(Hmisc)
nvar=228
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


hpf=brm(Time|cens(1-OS) ~ .,
        data = data,family = weibull(),
        prior=ssde,stanvars=stanvars,
        chain=3,iter=8000)

result<-summary(hpf)
resultt<-result[["fixed"]]
hpsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hpkfold <- kfold(hpf, chains = 3, K=5,
                 save_fits = TRUE)
hpreskfold<-hpkfold[["estimates"]]

hpkfp <- kfold_predict(hpkfold)
hpyhat<-hpkfp[["yrep"]]
hpyrep_mean <- colMeans(hpyhat)
hpcindex=rcorr.cens(hpyrep_mean,Surv(data$Time,data$OS))







nsnp=111
nvar=117
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


hphf=brm(Time|cens(1-OS) ~ .,
         data = data,family = weibull(),
         prior=eefssde,stanvars=eefstanvars,
         chain=3,iter=8000)

result<-summary(hphf)
resultt<-result[["fixed"]]
hphsigvar<-resultt[which(resultt[,3]<0 & resultt[,4]<0 | resultt[,3]>0 & resultt[,4]>0),]

hphkfold <- kfold(hphf, chains = 3, K=5,
                  save_fits = TRUE)
hphreskfold<-hphkfold[["estimates"]]

hphkfp <- kfold_predict(hphkfold)
hphyhat<-hphkfp[["yrep"]]
hphyrep_mean <- colMeans(hphyhat)
hphcindex=rcorr.cens(hphyrep_mean,Surv(data$Time,data$OS))








hreskfold
ehreskfold
hpreskfold
hphreskfold

hcindex
ehcindex
hpcindex
hphcindex


#plot(f, variable = "b_sIFNE")
#plot(ff, variable = "b_sIFNE")
#plot(hf, variable = "b_sIFNE")
#plot(fit, variable = "b_sIFNE")


#plot(f, variable = "b_sSH2D5")
#plot(ff, variable = "b_sSH2D5")
#plot(hf, variable = "b_stage")
#plot(fit, variable = "b_sSH2D5")




