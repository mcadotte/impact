source("ext.stat.R")

#if you want to create randomized log-normal rank-abundance curves
#ab<-abs(rlnorm(20,1,1))
#ab<-(ab/max(ab))*100

#single example
ab<-c(35.748064, 10.181998, 66.059484, 9.394792, 6.352681, 6.747668, 42.558520, 34.057773, 88.959203, 17.004277, 11.271972, 53.783547, 22.791321, 10.510519, 77.718661, 3.623766, 100.000000, 2.933351, 5.476969, 32.745390)

c1<-data.frame(Abun=ab[order(ab,decreasing=TRUE)],Rank=1:length(ab))

plot(c1$Rank,c1$Abun,type="h",xlab="Rank",ylab="Abundance",col="grey25",lwd=20,cex.lab=1.3,main="Community rank-abundance",cex.main=1.3)

###two scenarios of expected loss, inv = 200 and inv = 400

exp.ext(c1$Abun,200)
exp.ext(c1$Abun,400)

#example for 10 randomizations. useful for histograms or other analyses
exp.ext.rnd<-ext.rnd(c1$Abun,300,times=10)

#usage, scenario 1, same extinctions and rank not different
res1<-ext.sum.stat(c1$Abun,400,times=999,1,c(19))
res1<-as.data.frame(res1)

#usage, scenario 2, same extinctions and rank not different
res2<-ext.sum.stat(c1$Abun,400,times=999,3,c(18,19,20))
res2<-as.data.frame(res2)

#usage, scenario 3, same extinctions but rank different
res3<-ext.sum.stat(c1$Abun,400,times=999,2,c(8,14))
res3<-as.data.frame(res3)

#usage, scenario 4, more extinctions and rank different
res4<-ext.sum.stat(c1$Abun,400,times=999,7,c(2,4,6,8,10,14,15))
res4<-as.data.frame(res4)

res<-rbind(res1,res2,res3,res4)
rownames(res)<-c("Scenario 1","Scenario 2","Scenario 3","Scenario 4")
res

#write.csv(res,"~/Dropbox/2018WorkingFiles/Marc2018/Invasive species impact/R/simple-example.csv")