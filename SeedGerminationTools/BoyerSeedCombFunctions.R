setwd("/home/bregnery/Documents/CaughlinRegnerySimulation/SeedGerminationTools")

#install.packages("spam_0.29-3.tar.gz", repos=NULL, type="source")

################THIS READS IN THE PARAMETER VALUES FROM THE STATISTICAL ANALYSIS  
# IBMpars<-read.csv("IBMpars11.27.csv",header=T)
IBMpars<-read.csv("IBMpars12.16.csv",header=T)  
###############FOR THE TEST RUNS, I have been just looking at the median parameter values
med.par<-lapply(IBMpars,median)
attach(med.par)

load("Inits200m.Rdata")
initial.conditions.94<-Inits200m #should be a list

inits<-initial.conditions.94

#set the boundaries of the plot
L<-200

#initialize the dispersal kernel
L2<-L*2

nddist<-20

library("spam")

smatADULTseeds<-
  function(adults,seeds,alpha.est,dis.est,alpha.sg,dis.sg) #takes as input the entire matrix of seedlings
    ###need to change this so only adults >20
    
    #change to reflect neighborhood radius for seedlings?
    #change main code so that it takes the output from the list as a dollar sign
    
  {
    adults<-as.matrix(adults)
    size<-adults[,1]
    if(length(adults)==4) {
      x.off<-seeds[,1]
      y.off<-seeds[,2]
      dis1<-dis.fx.torus(adults[2],adults[3],x.off,y.off)
      surv<-ifelse(dis1<nddist,alpha.est*(adults[1]/dis1^dis.est),0)
      sgro<-ifelse(dis1<nddist,alpha.est*(adults[1]/dis1^dis.sg),0)
      return(list(nddEST=surv,sgro=sgro))}
    else {
      
      ads1<-adults[,2:3]
      snmat<-nearest.dist(seeds,ads1,method="toroidal",delta=nddist,L=L)
      sizemat<-snmat
      pos<-triplet(sizemat,tri=TRUE)
      sizemat@entries=size[pos$j]
      
      
      
      snmat@entries<-snmat@entries+0.0001
      snmat@entries<-snmat@entries^dis.est
      nddEST<-alpha.est*rowSums(sizemat/snmat,na.rm=T) #####DOUBLE CHECK THIS sum or rowSUMS more appropriate?
      
      snmat@entries<-snmat@entries^(1/dis.est)
      snmat@entries<-snmat@entries^dis.sg
      
      sgro<-alpha.sg*rowSums(sizemat/snmat,na.rm=T) #####DOUBLE CHECK THIS sum or rowSUMS more appropriate?
      
      
      return(list(nddEST=nddEST,sgro=sgro))
    }
  }

#maybe apply this to boyer model?
torus<-function(x,y,L) {
  size=1
  negL<-L*(-1)
  s1<-ifelse(x>2*L,0,size)
  s2<-ifelse(y>2*L,0,s1) #no long distance dispersal. assuming that if size goes to zero, then the plant is considered DEAD
  s3<-ifelse(x<negL,0,s2)
  s4<-ifelse(y<negL,0,s3)
  
  x1<-ifelse(x>L,x-L,x)
  y1<-ifelse(y>L,y-L,y)
  x3<-ifelse(x<0,x1+L,x1)
  y3<-ifelse(y<0,y1+L,y1)
  return(cbind(x3[which(s4==1)],y3[which(s4==1)]))
}




#########calculating neighborhoods around seeds
#########calculating neighborhoods around seeds
#########calculating neighborhoods around seeds

#adult neighborhoods
#this function calculates distance between two points in a torus
minitor<-function(x1,x2) {pmin(abs(x1-x2),L-abs(x1-x2))}

#calculates distance for a pair of two points
dis.fx.torus<-function(x1,y1,x2,y2) {sqrt(minitor(x1,x2)^2+minitor(y1,y2)^2)}

#creates the distance matrices
dmatSEEDS<-function(ads1,x.off,y.off)
{
  ads1=ads1
  x.off=x.off
  y.off=y.off
  
  dismat<-sqrt(outer(x.off,ads1[,1], "minitor")^2 + outer(y.off,ads1[,2],  #assuming that rows 2 and 3 are correct
                                                          "minitor")^2)
  
  return(dismat)}

#applies distance matrix to create vectors of survival and initial seedling establishment
#applies distance matrix to create vectors of survival and initial seedling establishment
#applies distance matrix to create vectors of survival and initial seedling establishment



#########calculating seedling neighborhoods around seeds
#########calculating seedling neighborhoods around seeds
#########calculating seedling neighborhoods around seeds
#########calculating seedling neighborhoods around seeds
SEEDsSLINGSnhoodNDD<-function(seeds,test)
{
  
  xvalz<-paste("(",floor(seeds[,1]),",",ceiling(seeds[,1]),"]",sep="")
  yvalz<-paste("(",floor(seeds[,2]),",",ceiling(seeds[,2]),"]",sep="")
  
  seedlocs<-paste(xvalz,yvalz,sep="+")
  
  
  
  
  xvalzSEEDLINGS<-paste("(",floor(test[which(test[,1]!=0),2]),",",ceiling(test[which(test[,1]!=0),2]),"]",sep="")
  yvalzSEEDLINGS<-paste("(",floor(test[which(test[,1]!=0),3]),",",ceiling(test[which(test[,1]!=0),3]),"]",sep="")
  
  
  seedlinglocs<-paste(xvalzSEEDLINGS,yvalzSEEDLINGS,sep="+")
  
  
  if(length(seedlocs) != 0) {
    testy<-rep(0,times=length(seedlocs))
    
    for(i in 1:length(seedlocs) ){
      testy[i]<-sum(seedlinglocs==seedlocs[i])
    }
    return(testy)
    
  } else{ return(0)}
}


library("plotrix")

treelocs<-function(adults) {
  
  
  trees1<-adults[,2]
  trees2<-adults[,3]
  
  trees<-cbind(trees1,trees2)
  
  plot(trees[,1]~trees[,2])
  return(trees)
}

tree.distance<-function(x,y,treeDAT) {
  
  treeDAT<-data.frame(treeDAT)
  
  n<-nrow(treeDAT)
  
  distance<-rep(NA,times=n)
  for(i in 1:n) {
    distance[i]<-(sqrt((x-treeDAT[i,1])^2+(y-treeDAT[i,2])^2))/treeDAT[i,3]
  }
  return(distance)
}

seed.distancef<-function(x,y,treeDAT) {
  
  treeDAT<-data.frame(treeDAT)
  
  n<-nrow(treeDAT)
  
  distance<-rep(NA,times=n)
  for(i in 1:n) {
    distance[i]<-(sqrt((x-treeDAT[i,1])^2+(y-treeDAT[i,2])^2))
  }
  return(distance)
}

Max.distance<-function(x,y,xp,yp) {
  distance<-(sqrt((x-xp)^2+(y-yp)^2))
  return(distance)
}

#returns angle with respect to the positive x axis
path.theta<-function(x,xp,y,yp,dist){
  theta<-(asin((yp-y)/dist))
  xdiff<-(xp-x)
  ydiff<-(yp-y)
  if(xdiff < 0 && ydiff > 0){
    theta = 3.141592654 - theta
  }
  else if(xdiff < 0 && ydiff < 0){
    theta = -3.141592654 - theta
  }
  else{
    theta = theta
  }
  return(theta)
}

moveTime<-function(velocity, x, y, i){
  
  movementTime<-(sqrt((x[i-1]-x[i])^2+(y[i-1]-y[i])^2))/velocity
  return(movementTime)
}

location<-function(times, x, y, i) {
  moveMagnitude = sqrt(x^2 + y^2)
  
  cosTheta = x/moveMagnitude
  sinTheta = y/moveMagnitude
  fruit<-matrix(runif(1*times,min=-10,max=10),nrow=times,ncol=1)
  moveMagnitude = moveMagnitude + fruit[i]
  x = moveMagnitude*cosTheta
  y = moveMagnitude*sinTheta
  coordinates<-c(x,y)
  return(coordinates)
}

adultR<-inits[[1]]

#if the movement distance is less than the step distance, there will be no randomness, the animal
#     will move directly between trees

Boyer<-function(times, adults, stepdist, ANG=0.0872,...) {
  tree1<-treelocs(adults)
  
  treeSIZE<-adults[,1]
  
  treeL<-vector("list",length=length(adults))
  
  treeL[[1]]<-cbind(tree1,treeSIZE)
  
  
  x=rep(NA,times)
  xp=rep(NA,times)
  y=rep(NA,times)
  yp=rep(NA,times)
  time.total = rep(0,1)
  time.near = rep(0,1)
  time.fraction = rep(0,1)
  time.movement = rep(0,1)
  time.move = rep(NA,times)
  time.current = rep(NA,times)
  time.current[1] = 0
  time.poop = rep(NA,100)
  seed.x=rep(NA,100)
  seed.y=rep(NA,100)
  seed.distance=rep(NA,100)
  theta.Stoc<-runif(times,min=-ANG,max=ANG)
  
  x[1]<-tree1[41,1]
  y[1]<-tree1[41,2]
  
  for (i in 2:times){
    treeM<-treeL[[i-1]]
    print(i)
    if(dim(t(data.frame(treeM)))[1]==1) {
      
      treeM<-t(data.frame(treeM))
      d=tree.distance(x[i-1],y[i-1],treeM)
      dis.sort<-sort(d)
      
      xp[i]=treeM[which(d==dis.sort[1]),1]
      yp[i]=treeM[which(d==dis.sort[1]),2]
      dist=Max.distance(x[i-1],y[i-1],xp[i],yp[i])
      if(dist > stepdist){
        theta=path.theta(x[i-1],xp[i],y[i-1],yp[i],dist)
        x[i]=x[i-1]+stepdist*cos(theta+theta.Stoc[i])
        y[i]=y[i-1]+stepdist*sin(theta+theta.Stoc[i])
      }
      else{
        x[i]=xp[i]
        y[i]=yp[i]
      }
      treeL[[i]]<-treeL[[1]][-which(d==dis.sort[1]),]
    }
    else{
      
      
      d=tree.distance(x[i-1],y[i-1],treeM)
      dis.sort<-sort(d)
      xp[i]=treeM[which(d==dis.sort[1]),1]
      yp[i]=treeM[which(d==dis.sort[1]),2]
      dist=Max.distance(x[i-1],y[i-1],xp[i],yp[i])
      if(dist > stepdist){
        theta=path.theta(x[i-1],xp[i],y[i-1],yp[i],dist)
        x[i]=x[i-1]+stepdist*cos(theta+theta.Stoc[i])
        y[i]=y[i-1]+stepdist*sin(theta+theta.Stoc[i])
      }
      else{
        x[i]=xp[i]
        y[i]=yp[i]
      }
      vist.dist<-Max.distance(x[i],y[i],xp[i],yp[i])
      if(vist.dist<10){
        treeL[[i]]<-treeM[-which(d==dis.sort[1]),]
      }
      else{
        treeL[[i]]<-treeM
      }
    }
    time.move[i]=moveTime(0.5,x,y,i)
    time.total= time.total + time.move[i] + 30
    time.current[i]=time.total
  }
  
  #Poop time
  time.poop0<-(rexp(100,rate=0.004))
  time.poop<-sort(ifelse(time.poop0>time.total,time.total,time.poop0))
  
  #Seed location
  for(i in 1:100){
    l=max(which(time.current<=time.poop[i]))
    j=min(which(time.current>=time.poop[i]))
    time.near = time.poop[i] - time.current[l]
    time.movement = time.current[j]-time.current[l]
    time.fraction = time.near/time.movement
    seed.x[i] = x[l]+time.fraction*(x[j]-x[l])
    seed.y[i] = y[l]+time.fraction*(y[j]-y[l])
  }
  
  #Histogram for seed distace from nearest tree
  for(i in 1:100){
    treeM<-treeL[[1]]
    treeM<-t(data.frame(treeM))
    seed.distance[i]=min(seed.distancef(seed.x[i],seed.y[i],treeM))
  }
  #hist(seed.distance,breaks=10,main="Seed Distance to the Nearest Tree")
  treeM1<-treeL[[1]]
  plot(treeM1[,1],treeM1[,2],pch="",xlim=c(-1,200),ylim=c(-1,200), main="Memory-Based", xlab="X(meters)", ylab="Y(meters)")
  lines(x,y,type="l") #,...)
  points(x,y,pch=19,cex=1.5)
  
  for(j in 1:length(tree1)/2) {
    draw.circle(tree1[j,1],tree1[j,2],radius=10,border="green",lwd=2,col=NA)
  }
  
  #plot(treeM[,1],treeM[,2],pch=21,cex=1,bg="red",main="Seed Dispersal Plot", xlab="X(meters)",ylab="Y(meters)")
  points(seed.x,seed.y,cex=0.9,pch=19,col="brown")
  
  #seed.locs[,1]=seed.x
  #seed.locs[,2]=seed.y
  
  #print(time.total)
  return(cbind(seed.x,seed.y,seed.distance))
  #return(cbind(x,y))
}
