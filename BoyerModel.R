#last bug remaining causes the animal to stay at the same tree every so many time steps
library("plotrix")

treelocs<-function(n) {
  
  
  trees<-matrix(runif(2*n,min=0,max=100),nrow=n,ncol=2)
  
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

seed.distance<-function(x,y,treeDAT) {
  
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

Boyer<-function(times, n, stepdist, ANG=0.0872,...) {
  tree1<-treelocs(n)
  
  treeSIZE<-runif(n,min=20,max=100)
  
  treeL<-vector("list",length=times)
  
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
  
  x[1]<-tree1[1,1]
  y[1]<-tree1[1,2]
  
  for (i in 2:times){
    treeM<-treeL[[i-1]]
  print(i)
    if(dim(t(data.frame(treeM)))[1]==1) {
      
      treeM<-t(data.frame(treeM))
      d=tree.distance(x[i-1],y[i-1],treeM)
      dis.sort<-sort(d)
      
      #Prevents animal from returning to the same tree when the loop is rerun
      #if(dis.sort[1]==0.){
      #  treeL[[i]]<-treeL[[1]][-which(d==dis.sort[1]),]
      #  d=tree.distance(x[i-1],y[i-1],treeM)
      #  dis.sort<-sort(d)
      #}
      
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
    seed.distance[i]=min(seed.distance(seed.x[i],seed.y[i],treeM))
    }
  #hist(seed.distance,breaks=10,main="Seed Distance to the Nearest Tree")
  treeM1<-treeL[[1]]
  plot(treeM1[,1],treeM1[,2],pch="",xlim=c(-1,100),ylim=c(-1,100), main="Memory-Based", xlab="X(meters)", ylab="Y(meters)")
  lines(x,y,type="l") #,...)
  points(x,y,pch=19,cex=1.5)
  
  for(j in 1:n) {
    draw.circle(tree1[j,1],tree1[j,2],radius=10,border="green",lwd=2,col=NA)
  }
  
  #plot(treeM[,1],treeM[,2],pch=21,cex=1,bg="red",main="Seed Dispersal Plot", xlab="X(meters)",ylab="Y(meters)")
  points(seed.x,seed.y,cex=0.9,pch=19,col="brown")
  
  print(time.total)
  #return(seed.distance)
  return(cbind(x,y))
}

#at each time step, animal goes to tree with BIGGEST tree.distance function
#spends a certain amount of time there
#then goes to the next biggest tree
# does not return to the same tree for 100 time steps
Boyer(25,10,1000,ANG=0)


seed.hist<-matrix(NA,ncol=30,nrow=100)

for(i in 1:30) {
  seed.hist[,i]<-Boyer(25,10,ANG=0)
  
}

hist(c(seed.hist),col="pink",
     xlab="Distance to nearest conspecific (m)",
     main="Memory-based movement model")

#