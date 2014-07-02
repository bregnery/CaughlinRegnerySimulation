setwd("/home/bregnery/Documents/CaughlinRegnerySimulation/SeedGerminationTools") #set to whereever source code is stored

source("BoyerSeedCombFunctions.R")


#at each time step, animal goes to tree with BIGGEST tree.distance function
#spends a certain amount of time there
#then goes to the next biggest tree
# does not return to the same tree for 100 time steps
seed.locs<-na.omit(Boyer(25,adultR,1000,ANG=0))

#################running

news<-cbind(seed.locs[,1],seed.locs[,2]) #this will be where you put x and y coordinates of seeds from boyer model

nbor.seeds<-smatADULTseeds(adultR,news,alpha.est,dis.est,alpha.sg,dis.sg) #news is a matrix of x and y cos from seeds
#a caveat here is that seeds only use adults over 20 m only, but seedlings *should* be taking all adults

nborADULTS<-nbor.seeds$nddEST
nbor.sgro<-nbor.seeds$sgro

#seedling vector

#seedlings<-subset(old,old[,4]==2)
seedlings<-data.frame(rep(1,times=10),runif(10,min=10,max=100),runif(10,min=10,max=100))
nborSEEDLINGS<-SEEDsSLINGSnhoodNDD(news,seedlings)



inht.mean<-(mu.inht)

newRECRUITSsurvive1<-rbinom(length(news[,1]),size=1,prob=plogis(mu.est+nborADULTS+beta.est*nborSEEDLINGS))



#write code to loop over the boyermodel with x and y coordinates of seeds as output, then apply this code
nsim<-10

stoch<-seq(from=0,to=2*pi,by=0.1)

stoch.use<-rep(stoch,times=nsim)

nruns<-length(stoch.use)

seed.vec<-rep(NA,times=nruns)
for(i in 1:nruns) {
  seed.locs<-na.omit(Boyer(25,adultR,1000,ANG=stoch.use[i]))

  #################running

  news<-na.omit(as.matrix(cbind(unlist(seed.locs[,1]),unlist(seed.locs[,2])))) #this will be where you put x and y coordinates of seeds from boyer model

  nbor.seeds<-smatADULTseeds(adultR,news,alpha.est,dis.est,alpha.sg,dis.sg) #news is a matrix of x and y cos from seeds
  #a caveat here is that seeds only use adults over 20 m only, but seedlings *should* be taking all adults

  nborADULTS<-nbor.seeds$nddEST
  nbor.sgro<-nbor.seeds$sgro

  #seedling vector

  #seedlings<-subset(old,old[,4]==2)
  seedlings<-data.frame(rep(1,times=10),runif(10,min=10,max=100),runif(10,min=10,max=100))
  nborSEEDLINGS<-SEEDsSLINGSnhoodNDD(news,seedlings)

  inht.mean<-(mu.inht)

  newRECRUITSsurvive1<-rbinom(length(news[,1]),size=1,prob=plogis(mu.est+nborADULTS+beta.est*nborSEEDLINGS))
  seed.vec[i]<-sum(newRECRUITSsurvive1)
  
  #it is storing seed.locs[,3] as a function (tree.locs?) rather than the actual distance
  #na.omit(unlist(seed.locs[,3])),
  data1<-data.frame(rep(i,times=length(news[,1])),rep(stoch.use[i],times=length(news[,1])))
  output[[i]]<-data1
}
#Merge dataframes and name columns
outty<-do.call("rbind",output)
colnames(outty)<-c("distance","Run","Stochasticity")
seedSurvive<-data.frame(cbind(seed.vec, stoch.use))
colnames(seedSurvive)<-c("Survival","Stochasticity")

#Create a Box Plot for number survived vs. Stochasticity
plot(seedSurvive$Survival~as.factor(seedSurvive$Stochasticity),cex=0.8)


plot(jitter(seedSurvive$Survival)~seedSurvive$Stochasticity,cex=0.8)

curve(exp(1.086479+0.009284*x),add=T,col="red")

#Create a Box Plot for Seed Distance vs. Stocasticity
plot(outty$distance~as.factor(outty$Stochasticity),cex=0.8)









