library(readxl)
library(forecast)
dt <- read_excel("parametros_compuestos.xlsx", 
                 col_types = c("text", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "numeric", "numeric", 
                               "numeric", "text", "numeric", "text"))
inter <- read_excel("parametros_compuestos.xlsx",sheet = "Hoja2")

dt$`Excretion - oct2`=NULL
dt$`Distribution - CNS`=NULL

inter=inter[-c(18,14),]
inter$MIN=as.numeric(inter$MIN)
inter$MAX=as.numeric(inter$MAX)
inter[15,2:3]=NA

dt$`Metabolism - CYP3A4`=factor(dt$`Metabolism - CYP3A4`)
levels(dt$`Metabolism - CYP3A4`)=c(-1,1)
dt$`Metabolism - CYP3A4`=as.integer(as.character(dt$`Metabolism - CYP3A4`))

M=as.matrix(data.frame(dt[,-1]))
f=c(1,2,2,2,2,2,2,3,3,3,3,4,4,5,5,6)
lam=apply(abs(M),2,function(x) BoxCox.lambda(abs(x),method="loglik"))
lam[c(15,5,12,8,10,1)]=NA

Mt=
mapply(function(x,y) {
  if(!is.na(y)) {
  abs(x)^y*sign(x)
  } else {
    x
  }
},as.data.frame(M),lam)

sds=apply(Mt,2,sd)

Mt=t(t(Mt)/sds)
interT=(inter[,2:3]^ifelse(is.na(lam),1,lam))/sds

Mt2=
mapply(function(x,y) {
  x=na.omit(x)
  if(length(x)==0) {
    v=matrix(y,ncol=1)
    v
  } else if(length(x)==1) {
    v=cbind((y-x)*(y<x),(y-x)*(y>=x))
    colnames(v)=c("L1","R1")
    v
  } else {
    v=do.call(cbind,lapply(x,function(xx) cbind((y-xx)*(y<xx),(y-xx)*(y>=xx)) ))
    colnames(v)=c("L1","R1","L2","R2")
    v
  }
},apply(interT,1,unique),
as.data.frame(Mt),SIMPLIFY = F)
names(Mt2)=inter$Compounds
f0=rep(1:length(f),sapply(Mt2,ncol))
f2=rep(f,sapply(Mt2,ncol))
Mt2=as.data.frame(data.frame(Mt2))
ix=which(apply(Mt2,2,sd)!=0)
nm=colnames(Mt2)
f3=f2
Mt3=Mt2

f2=f2[apply(Mt2,2,sd)!=0]
Mt2=Mt2[,apply(Mt2,2,sd)!=0]

library(fastICA)
s=
lapply(unique(f2),function(i) {
  if(sum(f2==i)>1) {
    set.seed(1994)
    s=fastICA(Mt2[,f2==i,drop=F],1)
    s$K=s$K %*% s$W;s$W=NULL
    s$K=s$K/sd(s$S)*sign(-cor(s$S,Mt2$Binding_Energy))[1,1]
    s$S=as.matrix(Mt2[,f2==i,drop=F]) %*% s$K
    s$X=NULL;s$A=NULL
    s$I=cbind(which(f2==i),rep(i,nrow(s$K)))
    s
  } else {
    s=list()
    s$S=Mt2[,f2==i,drop=F]
    sg=sign(-cor(s$S,Mt2$Binding_Energy))[1,1]
    s$S=s$S*sg
    colnames(s$S)=NULL
    s$K=matrix(1,ncol=1)*sg
    s$I=cbind(which(f2==i),rep(i,nrow(s$K)))
    s
  }
})
names(s)=unique(gsub(" -.{1,}","",inter$Compounds))
library(reshape2)
s1=data.frame(lapply(s,function(x) x$S))
sk=acast(do.call(rbind,lapply(s,function(x) data.frame(cbind(x$K,x$I)))),X2~X3,value.var = "X1",fill = 0)
colnames(sk)=colnames(s1)
rownames(sk)=colnames(Mt2)

set.seed(1994)
s2=fastICA(s1,n.comp = 1)
s2$K=s2$K %*% s2$W;s2$W=NULL
s2$K=s2$K/sd(s2$S)*sign(-cor(s2$S,Mt2$Binding_Energy))[1,1]
s2$S=as.matrix(s1) %*% s2$K
s2$X=NULL;s2$A=NULL
rownames(s2$K)=colnames(s1)
KT=sk %*% s2$K

KT=KT[match(1:43,ix),,drop=F]
rownames(KT)=nm
KT[is.na(KT)]=0
sk=sk[match(1:43,ix),]
rownames(sk)=nm
sk[is.na(sk)]=0
rm(Mt2,f2)
Mt3$s=s2$S[,1]

lim=unlist(apply(inter[,-1],1,function(x) {
  x=unique(x)
  rep(x,ifelse(is.na(x),1,2))
}))

dt$s=s2$S[,1]
tr_inter=data.frame(A=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,36),
                    B=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,15,16,15,17,18,20,19,20,19,21,22,24,23,24,23,25,26,28,27,28,27,29,30,32,31,32,31,33,34,36,35,36,35,37,38,40,39,40,39,41,42,43))
tr_inter=acast(tr_inter,A~B,fill = 0,fun.aggregate = length)
nm2=c("Binding_Energy","Physicochemical.characteristics -MW - <500","Physicochemical.characteristics -MW - >=500","Physicochemical.characteristics -HBA - <5","Physicochemical.characteristics -HBA - >=5","Physicochemical.characteristics -HBD - <10","Physicochemical.characteristics -HBD - >=10","Physicochemical.characteristics -Log.Po.w. - <5","Physicochemical.characteristics -Log.Po.w. - >=5","Physicochemical.characteristics -RotB - <10","Physicochemical.characteristics -RotB - >=10","Physicochemical.characteristics -PSA - <140","Physicochemical.characteristics -PSA - >=140","Absorption -LogS - <-6,5","Absorption -LogS - [-6,5;0,5)","Absorption -LogS - >=0,5","Absorption -.HOA - <0","Absorption -.HOA - [0;100)","Absorption -.HOA - >=100","Absorption -LogKp - <-8","Absorption -LogKp - [-8;-1)","Absorption -LogKp - >=-1","Absorption -Caco - <25","Absorption -Caco - [25;500)","Absorption -Caco - >=500","Distribution -LogBB - <-3","Distribution -LogBB - [-3;1,2)","Distribution -LogBB - >=1,2","Distribution -MDCK - <25","Distribution -MDCK - [25;500)","Distribution -MDCK - >=500","Metabolism -.Metab - <1","Metabolism -.Metab - [1;8)","Metabolism -.Metab - >=8","Metabolism -CYP3A4","Excretion -TC")
rownames(tr_inter)=nm2
colnames(tr_inter)=nm

tr_inter %*% KT
tr_inter %*% sk

write.table(tr_inter %*% KT,sep="\t",file = "clipboard")




