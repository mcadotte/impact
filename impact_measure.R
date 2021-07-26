###calculating impact from range, abundance and effect (sensu Parker et al. 1999)

impact.score<-function(R,A,E){
  return((R+A+E)/3)
}

####Effect is average of Z-values on extinctions
#DSV.change<-read.csv("Rouge data/DSV_change.csv")
EFF<-mean(DSV.change$ext.z)

####Range is from SES of occupancy in time t+1
two<-read.csv("Rouge data/data2019_two.csv")
rownames(two)<- two[,1]
two<- two[,2:ncol(two)]

#turn into pres/ab matrix
two.pa<-two
two.pa[two.pa > 0]<-1

#calculate SES of occupancy
DSV.occ.obs<-sum(two.pa$VIRO)

occ.ses<-function(pa.mat,times=999,spp.name="VIRO"){
  occ.obs<-sum(pa.mat[,colnames(pa.mat)==spp.name])
  tmp.occ<-NULL
  for (i in 1:times){
    tmp<-t(apply(pa.mat,1,sample))
    colnames(tmp)<-colnames(pa.mat)
    tmp.occ[i]<-sum(tmp[,colnames(tmp)==spp.name])
  }
  return((occ.obs-mean(tmp.occ))/sd(tmp.occ))
}

RGE<-occ.ses(two.pa)

####Abundance is from SES of randomized abundances for each plot
ab.ses<-function(mat,times=999,spp.name="VIRO"){
  ab.obs<-mean(mat[,colnames(mat)==spp.name],na.rm=TRUE)
  tmp.ab<-NULL
  for (i in 1:times){
    tmp<-t(apply(mat,1,sample))
    colnames(tmp)<-colnames(mat)
    tmp.ab[i]<-mean(tmp[,colnames(tmp)==spp.name],na.rm=TRUE)
  }
  return((ab.obs-mean(tmp.ab))/sd(tmp.ab))
}

AB<-ab.ses(two)


#funal calc
impact.score(EFF,RGE,AB)
