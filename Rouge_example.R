source("ext.stat.R")

####data example using Rouge Park invasion data.

one.t<-read.csv("Rouge data/data2013_one.t.csv")
rownames(one.t)<- one.t[,1]
one.t<- one.t[,2:ncol(one.t)]

two<-read.csv("Rouge data/data2019_two.csv")
rownames(two)<- two[,1]
two<- two[,2:ncol(two)]

DSV.change<-data.frame(site.plot=rownames(one.t),
                       DSV.change=(two$VIRO-one.t$VIRO),
                       DSV.2013=one.t$VIRO,DSV.2019=two$VIRO)

#DSV.change<-read.csv("Rouge data/DSV_change.csv")

#wrap in a loop
#need function to pull out number of extinctions and their ranks from before/after data

##data should include both time periods with common columns

out13<-list()

for (i in 1:nrow(one.t)){
  tmp<-as.numeric(one.t[i,])
  names(tmp)<-colnames(one.t)
  tmp<-tmp[names(tmp)!="VIRO"]#get rid of VIRO
  tmp<-tmp[tmp>0]
  out13[[i]]<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
}

names(out13)<-rownames(one.t)


#show that VIRO(DSV) was negatively correlated with richness in 2013
rich2013<-sapply(out13,nrow) 
DSV2013<-one.t$VIRO
DSV.change$Rich.2013<-rich2013

#get list of communities for 2019 as well.
out19<-list()

for (i in 1:nrow(two)){
  tmp<-as.numeric(two[i,])
  names(tmp)<-colnames(two)
  tmp<-tmp[names(tmp)!="VIRO"]#get rid of VIRO
  tmp<-tmp[tmp>0]
  out19[[i]]<-data.frame(Abun=tmp[order(tmp,decreasing=TRUE)],Rank=1:length(tmp))
}

names(out19)<-rownames(two)
rich2019<-sapply(out19,nrow) 
DSV.change$Rich.2019<-rich2019

#plots richness of 2013 vs 2019; and DSV richness
quartz()
par(mfrow=c(1,3))
l.rich<-lm(DSV.change$Rich.2019~DSV.change$Rich.2013)
summary(l.rich)
plot(DSV.change$Rich.2013,DSV.change$Rich.2019,pch=19,cex=1.3,
     cex.lab=1.4,ylab="Richness (2019)",
     xlab="Richness (2013)")
text(2,18,pos=4,expression(paste(beta," = ","0.58")))
text(2,16.8,pos=4,"P < 0.001")
abline(l.rich,lwd=4,col="Grey")

l.2013<-lm(DSV.change$Rich.2013~DSV.change$DSV.2013)
summary(l.2013)
plot(DSV.change$DSV.2013,DSV.change$Rich.2013,pch=19,cex=1.3,
     cex.lab=1.4,xlab="DSV abundance (2013)",
     ylab="Richness (2013)")
text(70,17,pos=4,expression(paste(beta," = ","-0.03")))
text(70,16,pos=4,"P < 0.01")
abline(l.2013,lwd=4,col="Grey")


l.2019<-lm(DSV.change$Rich.2019~DSV.change$DSV.2019)
summary(l.2019)
plot(DSV.change$DSV.2019,DSV.change$Rich.2019,pch=19,cex=1.3,
     cex.lab=1.4,xlab="DSV abundance (2019)",
     ylab="Richness (2019)")
text(120,18,pos=4,expression(paste(beta," = ","-0.04")))
text(120,16.8,pos=4,"P < 0.001")
abline(l.2019,lwd=4,col="Grey")


#focus on extinctions

exts<-NULL
ranks<-NULL

for (i in 1:length(out13)){
  a<-out13[[i]]
  b<-out19[[i]]
  tmp.e<-is.na(match(rownames(a),rownames(b)))
  exts[i]<-sum(tmp.e)
  ranks[i]<-mean(a$Rank[tmp.e])
}
ranks[is.nan(ranks)]<-NA

DSV.change$obs.ext<-exts
DSV.change$obs.rank.ext<-ranks

##do ses on each community

mod.out<-list()

#calculating only for plots with increasing DSV
for (i in 1:length(out13)){
  tmp<-out13[[i]]
  if (DSV.change$DSV.change[i] <= 0) mod.out[[i]]<-NA
  if (DSV.change$DSV.change[i] > 0){
    mod.out[[i]]<-ext.sum.stat(tmp,DSV.change$DSV.change[i], 999,
                DSV.change$obs.ext[i],DSV.change$obs.rank.ext[i])
  }
  
}

names(mod.out)<-names(out13)

###put into dataframe, maybe in a complex way
bar<-mod.out[!is.na(mod.out)]
final<-sapply(bar,rbind)
final<-as.data.frame(final)
final<-t(final)
final<-as.data.frame(final)
final<-apply(final,2,as.numeric)
final<-as.data.frame(final)
rownames(final)<-names(bar)

#write.csv(final,"Rouge data/final_output.csv")

ext.final<-read.csv("Rouge data/final_output.csv")
#DSV.change$Scenario<-ext.final$scenario[match(DSV.change$site.plot,ext.final$Plots)]
names(ext.final)[1]<-"Plots"

###just analyze plots with DSV increase and with > 4 species in 2013.
#add ext.z to DSV_change
DSV.change$ext.z<-ext.final$ext.z.value[match(DSV.change$site.plot,ext.final$Plots)]

scens.ext<-NULL
scens.rank<-NULL

for (i in 1:nrow(ext.final)){
  if (ext.final$ext.z.value[i] < -1.96 ) scens.ext[i]<-"S.neg"
  else if (ext.final$ext.z.value[i] > -1.96 && ext.final$ext.z.value[i] < 1.96) scens.ext[i]<-"N"
  else if (ext.final$ext.z.value[i] > 1.96 ) scens.ext[i]<-"S.pos"

  if (is.na(ext.final$rank.z.value[i])== FALSE) {
    if (ext.final$rank.z.value[i] < -1.96 ) scens.rank[i]<-"S.neg"
    else if (ext.final$rank.z.value[i] > -1.96 && ext.final$rank.z.value[i] < 1.96) scens.rank[i]<-"N"
    else if (ext.final$rank.z.value[i] > 1.96 ) scens.rank[i]<-"S.pos"
  }
  if (is.na(ext.final$rank.z.value[i])== TRUE) scens.rank[i]<-"N"
}

ext.final$Scenario[scens.ext=="N"&scens.rank=="S.neg"]<-3
ext.final$Scenario[scens.ext=="S.pos"&scens.rank=="S.neg"]<-4
ext.final$Scenario[scens.ext=="S.neg"&scens.rank=="S.neg"]<-3
ext.final$Scenario[scens.ext=="S.neg"&scens.rank=="N"]<-1

DSV.change$Scenario<-ext.final$Scenario[match(DSV.change$site.plot,ext.final$Plots)]

DSV.change<-DSV.change[DSV.change$Rich.2013 > 4,]
DSV.change<-DSV.change[DSV.change$DSV.change > 0,]
#write.csv(DSV.change,"Rouge data/DSV_change.csv")

#average ezt.Z
mean(DSV.change$ext.z)

scn<-as.factor(as.character(DSV.change$Scenario))
scn.d<-data.frame(Scenario=c("1","2","3","4"),Counts=c(sum(scn=="1",na.rm = TRUE),sum(scn=="2",na.rm = TRUE), sum(scn=="3",na.rm = TRUE),sum(scn=="4",na.rm = TRUE)))
quartz()
barplot(scn.d$Counts,names.arg =scn.d$Scenario,xlab="Scenario",ylab="Frequency")

#install.packages("nnet",dependencies=TRUE)
library(nnet)

m1<-multinom(DSV.change$Scenario~DSV.change$DSV.change)
summary(m1)
z1<-summary(m1)$coefficients/summary(m1)$standard.errors
p1<-(1-pnorm(abs(z1),0,1))*2

dp1<-as.data.frame(predict(m1,type="probs"))
dp1$DSV.change<-DSV.change$DSV.change[as.numeric(rownames(dp1))]

quartz()
plot(dp1$DSV.change,dp1[,3])


m2<-multinom(DSV.change$Scenario~DSV.change$DSV.2013)
summary(m2)
z2<-summary(m2)$coefficients/summary(m2)$standard.errors
p2<-(1-pnorm(abs(z2),0,1))*2
#p2 not sig

m3<-multinom(DSV.change$Scenario~DSV.change$DSV.2019)
summary(m3)
z3<-summary(m3)$coefficients/summary(m3)$standard.errors
p3<-(1-pnorm(abs(z3),0,1))*2
##p3 sig  -but DSV2019 highly correlated with DSV.change

m4<-multinom(DSV.change$Scenario~DSV.change$Rich.2013)
summary(m4)
z4<-summary(m4)$coefficients/summary(m4)$standard.errors
p4<-(1-pnorm(abs(z4),0,1))*2
#p4 sig -not cor with DSV.change
cor.test( DSV.change$DSV.change,DSV.change$Rich.2013)

m5<-lm(DSV.change$obs.ext~DSV.change$Rich.2013)
summary(m5)


m6<-lm(DSV.change$DSV.change~DSV.change$DSV.2013)
summary(m6)


###best model, Aic lower than either m1 or m4 alone
m10<-multinom(DSV.change$Scenario~DSV.change$DSV.change+DSV.change$Rich.2013)
summary(m10)
z10<-summary(m10)$coefficients/summary(m10)$standard.errors
p10<-(1-pnorm(abs(z10),0,1))*2

dp10<-as.data.frame(predict(m10,type="probs"))
dp10$DSV.change<-DSV.change$DSV.change[as.numeric(rownames(dp10))]
dp10$Rich.2013<-DSV.change$Rich.2013[as.numeric(rownames(dp10))]

quartz()
par(mfrow=c(2,2))
hist(DSV.change$DSV.change,xlab="Change in DSV cover",
     main=NULL,cex.lab=1.4)
hist(DSV.change$Rich.change,xlab="Change in species richness",
     main=NULL,cex.lab=1.4)
#plot(DSV.change$DSV.change,DSV.change$Rich.change,cex.lab=1.4,type="n",
  #   ylab="Change in species richness", xlab="Change in DSV cover")
#rect(0,-7.5,155,0,col=rgb(red=218,green=240,blue=227,maxColorValue = 255))
#points(DSV.change$DSV.change,DSV.change$Rich.change,pch=19,cex=1.3)

barplot(scn.d$Counts,names.arg =scn.d$Scenario,xlab="Scenario",
        ylab="Frequency",cex.lab=1.4)

plot(dp10$Rich.2013,dp10[,1],type="n",xlab="Resident richness (2013)",
     ylab="Probability",cex.lab=1.4,ylim =c(0,1))
lines(lowess(dp10$Rich.2013,dp10[,1],f=1),lwd=3,col="#EA4D1C")
lines(lowess(dp10$Rich.2013,dp10[,2],f=1),lwd=3,col="#CECE34")
lines(lowess(dp10$Rich.2013,dp10[,3],f=1),lwd=3,col="#CB4ECE")

text(5,0.99,pos=4,"Scenario 1",col="#EA4D1C",cex.lab=1.2)
text(5,0.92,pos=4,"Scenario 3",col="#CECE34",cex.lab=1.2)
text(5,0.85,pos=4,"Scenario 4",col="#CB4ECE",cex.lab=1.2)

plot(dp10$DSV.change,dp10[,1],type="n",xlab="Change in DSV cover",
     ylab="Probability",cex.lab=1.4,ylim =c(0,1))
lines(lowess(dp10$DSV.change,dp10[,1],f=1),lwd=3,col="#EA4D1C")

l.s3<-lowess(dp10$DSV.change,dp10[,2],f=1)
l.s3$y[l.s3$y>0.99]<-0.99

lines(l.s3,lwd=3,col="#CECE34")

l.s4<-lowess(dp10$DSV.change,dp10[,3],f=1)
l.s4$y[l.s4$y<0.01]<-0.01
lines(l.s4,lwd=3,col="#CB4ECE")






