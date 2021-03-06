setwd("/home/bregnery/Documents/CaughlinRegnerySimulation/SeedGerminationTools") #set to whereever source code is stored

source("BoyerSeedCombFunctions.R")


#at each time step, animal goes to tree with BIGGEST tree.distance function
#spends a certain amount of time there
#then goes to the next biggest tree
# does not return to the same tree for 100 time steps
#pdf(file="test.pdf")
#seed.locs<-na.omit(Boyer(1000,adultR,10,100, ANG=pi))
#dev.off()
#################running

#news<-cbind(seed.locs[,1],seed.locs[,2]) #this will be where you put x and y coordinates of seeds from boyer model

#nbor.seeds<-smatADULTseeds(adultR,news,alpha.est,dis.est,alpha.sg,dis.sg) #news is a matrix of x and y cos from seeds
#a caveat here is that seeds only use adults over 20 m only, but seedlings *should* be taking all adults

#nborADULTS<-nbor.seeds$nddEST
#nbor.sgro<-nbor.seeds$sgro

#seedling vector

#seedlings<-subset(old,old[,4]==2)
seedlings<-Inits$seedlings
nborSEEDLINGS<-SEEDsSLINGSnhoodNDD(news,seedlings)



inht.mean<-(mu.inht)

newRECRUITSsurvive1<-rbinom(length(news[,1]),size=1,prob=plogis(mu.est+nborADULTS+beta.est*nborSEEDLINGS))



#write code to loop over the boyermodel with x and y coordinates of seeds as output, then apply this code
nsim<-60

stoch<-seq(from=0,to=pi,by=0.1)

stoch.use<-rep(stoch,times=nsim)

nruns<-length(stoch.use)

seed.vec<-rep(NA,times=nruns)
total.seeds<-rep(NA,times=nruns)

output<-vector("list",nruns)

pdf(file="gibbonruns.pdf")

for(i in 1:nruns) {
  seed.locs<-na.omit(Boyer(200,adultR,10,100, ANG=stoch.use[i]))

  #################running

  news<-na.omit(as.matrix(cbind(unlist(seed.locs[,1]),unlist(seed.locs[,2])))) #this will be where you put x and y coordinates of seeds from boyer model

  nbor.seeds<-smatADULTseeds(adultR,news,alpha.est,dis.est,alpha.sg,dis.sg) #news is a matrix of x and y cos from seeds
  #a caveat here is that seeds only use adults over 20 m only, but seedlings *should* be taking all adults

  nborADULTS<-nbor.seeds$nddEST
  nbor.sgro<-nbor.seeds$sgro

  #seedling vector

  #seedlings<-subset(old,old[,4]==2)
  seedlings<-Inits$seedlings
  nborSEEDLINGS<-SEEDsSLINGSnhoodNDD(news,seedlings)

  inht.mean<-(mu.inht)

  newRECRUITSsurvive1<-rbinom(length(news[,1]),size=1,prob=plogis(mu.est+nborADULTS+beta.est*nborSEEDLINGS))
  seed.vec[i]<-sum(newRECRUITSsurvive1)  
  total.seeds[i]<-length(newRECRUITSsurvive1)
  
  #it is storing seed.locs[,3] as a function (tree.locs?) rather than the actual distance
  #na.omit(unlist(seed.locs[,3])),
  data1<-data.frame(rep(i,times=length(news[,1])),rep(stoch.use[i],times=length(news[,1])))
  output[[i]]<-data1
}

#dev.off()

#Merge dataframes and name columns
outty<-do.call("rbind",output)
colnames(outty)<-c("Run","Stochasticity")
seedSurvive<-data.frame(cbind(seed.vec, total.seeds, stoch.use))
colnames(seedSurvive)<-c("Survival","total","Stochasticity")

#Create a Box Plot for number survived vs. Stochasticity

save(seedSurvive,file="seedSurvive.Rdata")

load(file="seedSurvive.Rdata")

plot(seedSurvive$Survival/seedSurvive$total~as.factor(seedSurvive$Stochasticity),cex=0.8)

plot(jitter(seedSurvive$Survival)~jitter(seedSurvive$Stochasticity),cex=0.8)

library("MASS")


m1<-glm.nb(seedSurvive$Survival~seedSurvive$Stochasticity)
m2<-glm(seedSurvive$Survival~seedSurvive$Stochasticity,family="poisson")
curve(exp(coef(m1)[1]+coef(m1)[2]*x),add=T,col="red",lwd=2)
curve(exp(coef(m2)[1]+coef(m2)[2]*x),add=T,col="blue",lwd=2)


m1<-glm.nb(seedSurvive$Survival~seedSurvive$Stochasticity+I(seedSurvive$Stochasticity^2))

curve(exp(coef(m1)[1]+coef(m1)[2]*x+coef(m1)[3]*x^2),add=T,col="red",lwd=2)

#Create a Box Plot for Seed Distance vs. Stocasticity
plot(outty$distance~as.factor(outty$Stochasticity),cex=0.8)
#dev.off()



curve(plogis(mu.est+((alpha.est/1)*((90)/(x)^(dis.est)))),from=0,to=20,ylim=c(0,0.03),cex.lab=1.4,cex.axis=1.2,xlab="Distance to neighbor tree (m)",
      ylab="Germination")

dev.off()

#