
  #to do: run initial states
  #to do: prevent huge seedlings. either have a mandatory fix to the switch function and have it so those large seedlings can't enter popdat as huge trees
  #(>100). Also, maybe repurpose basal area os that growth/reproduction ceases when >1
#to do: reparameterize seedlings so that exp(conspecific seedlings*seeds) is the norm
#keep in mind, reducing the size of the plot is a last-minute solution.
#I don't think basal area will fix things, considering that even 96 with the new basal area was 0.001


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


########################################
########################################
######UNDERLYING FUNCTIONS FOR IBM 7.29
########################################
########################################



########################################
########################################
######PRODUCTION OF NEW SEEDLINGS
########################################
########################################

#probability density function for dispersal distances
#input: u (parameterized dispersal parameter)
#x
PDFclark2dt<-function(x,u) {
p<-1
2*pi*x*(p/(3.141593*u))/(1+(x^2/u))^(p+1)
}

#fecundity: how many new seeds per tree
seedfx<-function(dbh,bf0,bf1,bf2) {
final<-rep(0,times=length(dbh))
for (i in 1:length(final)) {
final[i]<-ifelse(dbh[i]<20,0,
round(rpois(1, lambda=((bf0+1/(bf1+exp(bf2*dbh[i])))/9)*89)))
} #survival
return(final)}




nddist<-20

#########DISPERSAL:OUTPUT: WHERE DO SEEDS GO?
#########DISPERSAL:OUTPUT: WHERE DO SEEDS GO?
#########DISPERSAL:OUTPUT: WHERE DO SEEDS GO?
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



dispersal<-function(noffspring,x,y,adultsize,L,prob.distances) { #x and y of parent, number of offspring it makes
x.off<-rep(0,times=noffspring)
y.off<-rep(0,times=noffspring)  #this is incorrect, since it is sort of multiplying the dispersal distances by one another

a<-runif(noffspring, min =0, max = 2*pi)
dispersal.distance<-sample(dispersal.unif,size=noffspring,prob=prob.distances,replace=T)
dispersal.distanceF<-dispersal.distance+adultsize/200 #this is to correct for the size of the trunk
x.off<-x+cos(a)*dispersal.distanceF ###NOTE THAT THE WAY THE CODE IS HERE, it IS DOING a negative number.  so this could easily be negative
y.off<-y+sin(a)*dispersal.distanceF

noffs<-torus(x.off,y.off,L)
return(noffs)
}


disOUTLDD<-function(noffspring,adults,L,LDD,prob.distances) {
x<-adults[,2]
y<-adults[,3]
adultsize<-adults[,1]


listy<-vector("list",length(noffspring))

for(j in 1:length(noffspring)) {

listy[[j]]<-dispersal(noffspring[j],x[j],y[j],adultsize[j],L,prob.distances)
}
finality<-do.call("rbind",listy)
totalQ<-L*L
numberLDD<-rpois(totalQ,lambda=LDD)
nLDD<-sum(numberLDD)
if(nLDD<2) {return(finality)}
else{
LDDz<-cbind(runif(nLDD,min=0,max=L),runif(nLDD,min=0,max=L))
finality1<-rbind(LDDz,finality)
return(finality1)
}
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



  #embed this into the for loop for IBM: try looking at this with system time

        if(length(seedlocs) != 0) {
testy<-rep(0,times=length(seedlocs))

for(i in 1:length(seedlocs) ){
testy[i]<-sum(seedlinglocs==seedlocs[i])
}
return(testy)

    } else{ return(0)}
}


#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING
#########SEEDLING TO SEEDLING

##############Calculating adult neighborhoods around seedlings
 slingmatADULT<-function(adults,seedlings,alpha.sg,dis.sg) #takes as input the entire matrix of seedlings
###need to change this so only adults>20
{

if(length(as.matrix(seedlings))==4) {
x.off<-seedlings[2]
y.off<-seedlings[3]
seedling.size<-seedlings[1]}
else{
x.off<-seedlings[,2]
y.off<-seedlings[,3]
seedling.size<-seedlings[,1]
}
if(length(as.matrix(adults))==4) {
dis1<-dis.fx.torus(adults[2],adults[3],x.off,y.off)
stop(return(ifelse(c(dis1,dis1)<nddist,c((1*(adults[1]/dis1^1)),(alpha.sg*((adults[1]/seedling.size)/dis1^dis.sg))),0)))
}
else {

ads1<-adults[,2:3]
snmat<-nearest.dist(cbind(x.off,y.off),ads1,method="toroidal",delta=nddist,L=L)

size<-adults[,1]

  sizemat<-snmat
  pos<-triplet(sizemat,tri=TRUE)
  sizemat@entries=size[pos$j]
  snmat@entries<-snmat@entries+0.0001
  snmat@entries<-snmat@entries^dis.sg

nddEST<-exp((alpha.sg/seedling.size)*rowSums(sizemat/snmat,na.rm=T)) #####DOUBLE CHECK THIS sum or rowSUMS more appropriate?

return(nddEST)


}
}



##############Calculating seedling neighborhoods around seedlings
 match.func<-function(ind.val, vec.val) {
        vec.val[which(names(vec.val)==ind.val)]
 }


SEEDLINGnhoodNDD<-function(test)
{

        xvalz<-paste("(",floor(test[,2]),",",ceiling(test[,2]),"]",sep="")
        yvalz<-paste("(",floor(test[,3]),",",ceiling(test[,3]),"]",sep="")

        slocs<-paste(xvalz,yvalz,sep="+")
        sfac    <-unique(slocs)

  #embed this into the for loop for IBM: try looking at this with system time

        if(length(slocs) != 0) {

                testx <- table(slocs)
                turt    <- apply(as.matrix(slocs), 1, match.func, testx)
                return(unlist(turt)-1)

    } else{ return(0)}
}



##################SEEDLING SURVIVAL AND GROWTH KERNEL
pxSEEDLINGNDD<-function(ht,count,adz,mu.ss,beta.ss,size.ss,mu.sg,size.sg,beta.sg,shape.sg,scale.sg) { gt<-ifelse(ht==0,0,sxs(ht,count,mu.ss,beta.ss,size.ss)
*gxs(ht,count,adz,mu.sg,size.sg,beta.sg,shape.sg,scale.sg))
gt<-ifelse(gt<0,0,gt)
return(gt)}

sxs<-function(ht,count,mu.ss,beta.ss,size.ss) {return(rbinom(length(ht),size=1,prob=plogis(mu.ss + size.ss*ht+beta.ss*count)))}


gxs<-function(ht,count,adz,mu.sg,size.sg,beta.sg,shape.sg,scale.sg) {

mean.val<-(mu.sg+size.sg*ht)*exp(beta.sg*count)*exp(adz)
shape<-shape.sg #IS THIS CORRECT?
scale.val<-sqrt(scale.sg)



gro<-rsn(length(ht),location=mean.val-(scale.val*shape)/sqrt(1+shape^2)*sqrt(2/3.1415926),shape=shape,scale=scale.val)


final<-gro+ht
final1<-ifelse(final>230,230,final)
final2<-ifelse(final<0,0,final)
return(final2)
}


pxSEEDLINGNDD2<-function(ht,count,adz,mu.ss,beta.ss,size.ss,mu.sg,size.sg,beta.sg,shape.sg,scale.sg) { gt<-ifelse(ht==0,0,sxs2(ht,count,mu.ss,beta.ss,size.ss)
*gxs2(ht,count,adz,mu.sg,size.sg,beta.sg,shape.sg,scale.sg))
gt<-ifelse(gt<0,0,gt)
return(gt)}

sxs2<-function(ht,count,mu.ss,beta.ss,size.ss) {return(rbinom(length(ht),size=1,prob=plogis(mu.ss + size.ss*ht+beta.ss*count)^(8/12)))}


gxs2<-function(ht,count,adz,mu.sg,size.sg,beta.sg,shape.sg,scale.sg) {

mean.val<-(mu.sg+size.sg*ht)*exp(beta.sg*count)*exp(adz)*(8/12)
shape<-shape.sg #IS THIS CORRECT?
scale.val<-sqrt(scale.sg)



gro<-rsn(length(ht),location=mean.val-(scale.val*shape)/sqrt(1+shape^2)*sqrt(2/3.1415926),shape=shape,scale=scale.val)


final<-gro+ht
final1<-ifelse(final>230,230,final)
final2<-ifelse(final<0,0,final)
return(final2)
}


####SWITCH TO ADULTS. NOTE:  I am not considering this stochastic, because it doesn't reflect a biological process. sensitivity would be meaningless.
switch.fx<-function (slingline,aht,bht,sigma,a.switch,b.switch) {
ht<-slingline[1]
if(ht<250) {
return(slingline)
}
else {
switchorno<-rbinom(size=1,n=1,prob=plogis(a.switch+b.switch*(ht/100)))
if(switchorno>0) {
DBH<-rnorm(1,mean=aht+ht*bht,sd=sigma)
#DBH<-ifelse(DBH<1,1,DBH)
slingline[1]<-DBH
slingline[4]<-1
return(slingline)}
else{return(slingline)}
} #1 for adults
}



###############ADULTS

tordist<-function(x1,y1,x2,y2,X,Y) {
sqrt(min(abs(x1-x2),X-abs(x1-x2))^2+min(abs(y1-y2),Y-abs(y1-y2))^2)
}

minitor<-function(x1,x2) {pmin(abs(x1-x2),L-abs(x1-x2))}     #note: pmin takes the minimum of the vector

dist.fx<-function(x1,y1,x2,y2) {sqrt((x2-x1)^2+(y2-y1)^2)
}

#ok, this builds the distance matrix
dmat<-function(test)
{
nonzero<-which(test[,1]!=0)

dismat<-sqrt(outer(test[nonzero,2], test[nonzero,2], "minitor")^2 + outer(test[nonzero,3],  #assuming that rows 2 and 3 are correct
test[nonzero,3], "minitor")^2)

return(dismat)} #DOUBLECHECK

nddistADULTS<-25

omeg.dis<-10
nhoodNDD2<-function(adults,alpha.as,dis.as,alpha.ag,dis.ag)
{

ads1<-adults[,2:3]
dismat<-nearest.dist(ads1,ads1,method="toroidal",delta=nddistADULTS,L=L)

size<-adults[,1]

  sizemat<-dismat
  pos<-triplet(sizemat,tri=TRUE)
  sizemat@entries=size[pos$j]

  #not sure if this will work. also, you need to double check Nx. does it really not involve summation?
  omeg.mat<-dismat
  omeg.mat@entries<-ifelse(dismat@entries>omeg.dis,0,1)
  #it is necessary to have sum here, because
  Nx<-sum(omeg.mat@entries)

dismat@entries<-dismat@entries+0.0001
  dismat@entries<-dismat@entries^dis.as

diag(dismat)<-1
diag(sizemat)<-0

dbht<-size

finalS<-(alpha.as/(dbht))*rowSums(sizemat/dismat,na.rm=T)

  dismat@entries<-dismat@entries^(1/dis.as)
  dismat@entries<-dismat@entries^dis.ag


finalG<-(alpha.ag/(dbht))*rowSums(sizemat/dismat,na.rm=T)

gorp<-list(S=finalS,G=finalG,Nx=Nx)
#colnames(gorp)<-c("S","G","Nx")  #next, you need a way to save Nx in your for loop so that it can be saved in popdat. maybe assign()?

return(gorp)
}


dismatty<-function(ads1) {
dismat<-nearest.dist(ads1,ads1,method="toroidal",delta=nddistADULTS,L=L)

#size<-adults[,1]

#  sizemat<-dismat
#  pos<-triplet(sizemat,tri=TRUE)
#  sizemat@entries=size[pos$j]

  #not sure if this will work. also, you need to double check Nx. does it really not involve summation?
  omeg.mat<-dismat
  omeg.mat@entries<-ifelse(dismat@entries>omeg.dis,0,1)
  #it is necessary to have sum here, because
  Nx<-sum(omeg.mat@entries)
  return(Nx)
  }

##############SURVIVAL AND GROWTH KERNEL
hossfeld<-function(G,P,r,dbh) {
(G*r*dbh^(r-1))/(G+(dbh^r/P))^2}


plaw<-function(a,b,dbh) {a*dbh^b}


px<-function(x,SG,NG,G.g,P.g,G.s,P.s,a.sig,b.sig,a.shape,b.shape) { gt<-ifelse(x==0,0,sx(x,SG,G.s,P.s)*
gx(x,NG,G.g,P.g,a.sig,b.sig,a.shape,b.shape))
gt<-ifelse(gt>0 & gt<1,1,gt)
gt<-ifelse(gt<0,1,gt)
return(gt)}


#survival
sx<-function(x,NS,G.s,P.s) {
x<-x/10
rbinom(length(x),size=1,prob=plogis(hossfeld(exp(G.s),exp(P.s),2,x)+(NS/100))^(1/5))   #TURN THIS OFF IF YOU WANT THE PLAW AGAIN
}


#growth
gx<-function(dbh.old,NG,G.g,P.g,a.sig,b.sig,a.shape,b.shape) {
tree.mean<-hossfeld(G.g,exp(P.g),2,dbh.old)*exp(NG)
shape.ag=b.shape*exp(a.shape*dbh.old)
scale.ag=sqrt(a.sig*dbh.old^b.sig)


gro<-rsn(length(dbh.old),location=tree.mean-(scale.ag*shape.ag)/sqrt(1+shape.ag^2)*sqrt(2/3.1415926),
scale=scale.ag,shape=shape.ag)
final<-gro/5+dbh.old
final<-ifelse(final>100,100,final)
return(final)
}




bobles<-function(old,n) {

borders1<-c(0,L,L,0)
  borders2<-c(0,0,L,L)

  plot(borders1,borders2,main=as.character(paste("population at time",n)),cex=0.1)

filly<-old
  #filly<-subset(filler,filler[,2]==n)[,3:6]
  total<-cbind(borders1,borders2)
  polygon(total,col="white",lty=2)

  seedlings<-subset(filly,filly[,4]==2)

  points(seedlings[,2],seedlings[,3],pch=19,col="green",cex=0.02)

  adults<-subset(filly,filly[,4]==1)
  radius<-sqrt(adults[,1]/pi)
  symbols(adults[,2],adults[,3],circles=radius,inches=FALSE,fg="black",bg="red",add=T)
  }
  
  
bobNOW<-function(old,reading=T) {

borders1<-c(0,L,L,0)
  borders2<-c(0,0,L,L)

  plot(borders1,borders2,cex=0.1)


  if(reading==T) {filly<-old[,3:6]}
   else {filly<-old }
  
  
  total<-cbind(borders1,borders2)
  polygon(total,col="white",lty=2)

  seedlings<-subset(filly,filly[,4]==2)

  points(seedlings[,2],seedlings[,3],pch=19,col="green",cex=0.02)

  adults<-subset(filly,filly[,4]==1)
  radius<-sqrt(adults[,1]/pi)
  symbols(adults[,2],adults[,3],circles=radius,inches=FALSE,fg="black",bg="red",add=T)
  }
  
    
   popDATout<-function(filler,Nx)
  {
  seedlingCOUNT<-length(which(filler[,4]==2))
  adultCOUNT<-length(which(filler[,4]==1))
  total<-sum(seedlingCOUNT+adultCOUNT)
  basal.area<-sum(0.00007854*((filler[which(filler[,4]==1),1]/2)^2))/(L^2)
  Ax<-pi*omeg.dis^2*adultCOUNT
omega<-(Nx/Ax)/(adultCOUNT/(L*L))

  
  return(t(c(seedlingCOUNT,adultCOUNT,total,basal.area,omega)))
  }
  
