#setwd("C:/Users/Trevor/Dropbox/R files/Thailand/2010-2011/population stuff/1.28.2013")


  source("newFUNCTIONS 2.21.r")

  #
  #
  ##
  #
  #

  ###################
  #
  #
  #INITIALIZING DISPERSAL KERNEL
  #setwd("C:/Users/Trevor/Dropbox/R files/Thailand/2010-2011/population stuff/1.28.2013")
   #setwd("C:/Users/trevorcaughlin/Dropbox/R files/Thailand/2010-2011/population stuff/1.28.2013")
  ###########loading initial conditions
  load("Inits200m.Rdata")
  initial.conditions.94<-Inits200m

  #set the boundaries of the plot
  L<-200

  #initialize the dispersal kernel
  L2<-L*2
  dispersal.unif<-runif(90000,min=0,max=L2)
    ###############FOR LOOP
    #
    #
    #
    #
    #
    #
    #
    ############################

    IBM2<-function(run,initial,wantfull="Y",path,filePOPname,fileALLname,
    #############################SEED DISPERSAL
    u.dispersal,LDD,bf0,bf1,bf2,
    ############################SEED ESTABLISHMENT
    mu.est,alpha.est,dis.est,beta.est,
    ############################Initial seedling height
    mu.inht,scale.inht,shape.inht,
      ############################seedling survival
    mu.ss,beta.ss,size.ss,
      ############################seedling growth
    alpha.sg,dis.sg,mu.sg,size.sg,beta.sg,scale.sg,shape.sg,
    ############################switch from seedlings to adults
    aht,bht,sigma,a.switch,b.switch,
      ############################adult survival
      alpha.as,dis.as,P.s,G.s,
      ############################adult growth
       alpha.ag,dis.ag,P.g,G.g,a.sig,b.sig,a.shape,b.shape) {


 #dir.create(paste("C:/Users/Trevor/Dropbox/",path,sep=""))
# setwd(paste("C:/Users/Trevor/Dropbox/",path,sep=""))


   #########setting initial conditions
   #########setting initial conditions
   startz<-initial.conditions.94[[initial]]
   old<-cbind(startz,rep(1,times=length(startz[,1])))
   inits<-cbind(startz,rep(1,times=length(startz[,1])))

   Nx<-dismatty(inits)

  ###sets up vector of dispersal probabilities. this does not need to be run at each time step.
   prob.distances<-PDFclark2dt(dispersal.unif,u.dispersal)


    for (i in 2:run)
    {
     print(i)
  write.table(popDATout(old,Nx),filePOPname,append=T,sep=",",col.names=FALSE) #RECORD output for sensitivity JAKE: how do you want this out?
     ##############WHAT DO YOU WANT OUT? #do you want it to record the population every five years or not?
     if(wantfull=="Y") {

   if( any(i==c(seq(from=2,to=run,by=5),run))) {

    write.table(cbind(rep(i-1,times=length(old[,1])),old),file=paste(fileALLname,"-",i,".csv",sep=""),append=F,sep=",",col.names=FALSE) #RECORD location and size of every individual
   }
   else{ print(as.character(i))}
   }
   else{next}

    ####THIS SECTION GIVES SURVIVAL AND GROWTH OF EXISTING TREES
    if(sum(old[,1],na.rm=T)==0) { #if there are only zeroes, you only get zeroes
    return(cat("everything is dead")) #this should be faster than having it run through the entire for loop
  #again and again if there are only zeroes
    #filler[,1,i]<-rep(0,times=length(filler[,1,i])) # an alternative that fills it up with zero
    #again and again

    } else{



  filler<-matrix(nrow=length(old[,1]),ncol=length(old[1,]))

  new.mat<-filler

  old[,1]<-ifelse(old[,1]<=0,0,old[,1])


    ##################
    ################## Transition 1:
    ##################
    ################## ADULT to ADULT:
    if(i<8) {adultTRANS<-inits} else{
    if(length(which(old[,4]==1))==0) {adultTRANS<-c(0,0,0,0)}
    else {

    adults<-subset(old,old[,4]==1) #remember, 1=adults
    adults<-subset(adults,adults[,1]!=0)


   if(dim(adults)[1]==1){
   NDD<-data.frame(t(c(0,0)))
   colnames(NDD)<-c("S","G")
   }
   else{
    NDD<-nhoodNDD2(adults,alpha.as,dis.as,alpha.ag,dis.ag)
    }
    adultvec<-px(adults[,1],NDD$S,NDD$G,G.g,P.g,G.s,P.s,a.sig,b.sig,a.shape,b.shape)

    assign("Nx",NDD$Nx)

    if(length(adultvec)==1) {adultTRANS<-t(c(adultvec[1],adults[,2:4]))} else{
    adultTRANS<-cbind(adultvec,adults[,2:4])

    #adultTRANS<-adultTRANS[which(adultvec!=0),] #getting rid of all the zeroz

    }

    }
    }
    ##################
    ################## Transition 2:
    ##################
    ################## Seedling to Seedling:


    if(length(which(old[,4]==2))==0) {seedlingTRANSandNEWadults<-matrix(0,nrow=1,ncol=4)} else {

    seedlings<-subset(old,old[,4]==2) #remember, 2=seedlings

    seedlingCOUNT<-SEEDLINGnhoodNDD(seedlings)

    adults<-old[which(old[,4]==1),]

  if(sum(adults)==0) {adzNDD=0}
  else{

  adzNDD<-slingmatADULT(adults,seedlings,alpha.sg,dis.sg)
  }

  seedlingvec<-pxSEEDLINGNDD(seedlings[,1],seedlingCOUNT,adzNDD,mu.ss,beta.ss,size.ss,mu.sg,size.sg,beta.sg,scale.sg,shape.sg)

    if(length(seedlingvec)==1) {

    firstDOthis<-t(c(seedlingvec,seedlings[,2],seedlings[,3],2))
    seedlingTRANSandNEWadults<-switch.fx(firstDOthis,aht,bht,sigma,a.switch,b.switch)
    }
    else{
    seedlingTRANS<-cbind(seedlingvec,seedlings[,2:4])

    seedlingTRANS<-seedlingTRANS[which(seedlingvec!=0),] #removing dead ones
    seedlingTRANS<-seedlingTRANS[which(seedlingTRANS[,1]!=0),]
    ###########Transition 3: seedling to adult transition:
    if(sum(seedlingTRANS[,4])==2) {seedlingTRANSandNEWadults<-switch.fx(seedlingTRANS,aht,bht,sigma,a.switch,b.switch)}
    else{
    seedlingTRANSandNEWadults<-t(apply(seedlingTRANS,1,switch.fx,aht=aht,bht=bht,sigma=sigma,a.switch=a.switch,b.switch=b.switch)) #applying transition from seedlings to adults
    }
    }
    }


    ##################
    ################## Transition 4: Fecundity-production of new seedlings
    ##################


    if(length(which(old[,4]==1 & old[,1]>20))==0) {newRECRUITS<-c(0,0,0,0)}
    else {
    if(length(which(old[,1]>20))==0) {newRECRUITS<-c(0,0,0,0)}
      else{
      adultA<-subset(old,old[,4]==1)
    if(length(which(adultA[,1]>20))==0) {newRECRUITS<-c(0,0,0,0)} else {
     #remember, 1=adults
    adultA<-subset(old,old[,4]==1)
    adultR<-subset(adultA,adultA[,1]>20)
    recruits<-seedfx(adultR[,1],bf0,bf1,bf2) #number of new recruits per tree
  }
    if(sum(recruits)==0 | is.na(sum(recruits))==T) {newRECRUITS<-c(0,0,0,0)}

    else{

    #
    if(length(which(recruits!=0))==1 | length(adultR)==4) {adultREPRO<-t(adultR[which(recruits!=0),])} #annoyingly, if there is just one tree, the dimensions of the matrix are wrong and there are no columns
    else{
        adultREPRO<-as.matrix(adultR[which(recruits!=0),])
    }

    #adultREPRO<-as.matrix(adultR[which(recruits!=0),])

    news<-disOUTLDD(t(recruits[which(recruits!=0)]),adultREPRO,L,LDD,prob.distances)

    nbor.seeds<-smatADULTseeds(adultR,news,alpha.est,dis.est,alpha.sg,dis.sg)
    #a caveat here is that seeds only use adults over 20 m only, but seedlings *should* be taking all adults

    nborADULTS<-nbor.seeds$nddEST
    nbor.sgro<-nbor.seeds$sgro

    if(length(which(old[,4]==2))==0) {nborSEEDLINGS<-rep(0,times=length(news[,1]))}
    else {
    seedlings<-subset(old,old[,4]==2)
    nborSEEDLINGS<-SEEDsSLINGSnhoodNDD(news,seedlings)

    }

    inht.mean<-(mu.inht)

    newRECRUITSsurvive1<-rbinom(length(news[,1]),size=1,prob=plogis(mu.est+nborADULTS+beta.est*nborSEEDLINGS))
    newRECRUITSsize1<-rsn(length(newRECRUITSsurvive1),location=inht.mean-(scale.inht*shape.inht)/sqrt(1+shape.inht^2)*sqrt(2/3.1415926),
    scale=scale.inht,shape=shape.inht)

    newRECRUITSest<-newRECRUITSsurvive1*newRECRUITSsize1

    newRECRUITSsurvive<-pxSEEDLINGNDD2(newRECRUITSest,nborSEEDLINGS,nbor.sgro,mu.ss,beta.ss,size.ss,mu.sg,size.sg,beta.sg,shape.sg,scale.sg)
    newRECRUITSsize<-newRECRUITSsurvive

    #print(newRECRUITSsize)
    if(length(which(newRECRUITSsurvive1==1))==0) {newRECRUITS<-c(0,0,0,0)
    }
    else{
    if(length(which(newRECRUITSsurvive1==1))==1) {
    newRECRUITS<-t(c(newRECRUITSsize[which(newRECRUITSsurvive!=0)],news[which(newRECRUITSsurvive!=0),],
    rep(2,times=length(which(newRECRUITSsurvive!=0)))))
    }
    else{
    newRECRUITS<-cbind(newRECRUITSsize[which(newRECRUITSsurvive!=0)],news[which(newRECRUITSsurvive!=0),],
    rep(2,times=length(which(newRECRUITSsurvive!=0))))
    }
    }
    }
    }
    }


    if(sum(c(newRECRUITS,as.matrix(seedlingTRANSandNEWadults),as.matrix(adultTRANS)))==0) {
    return(cat("everything is dead"))}
    else{
    new.mat<-rbind(newRECRUITS,as.matrix(seedlingTRANSandNEWadults),adultTRANS)


  #outs<-c(which(new.mat[,1]<=0))

  ins<-c(which(new.mat[,1]>0))
  new.mat<-new.mat[ins,]
  if(length(new.mat)==4) {
    assign("old",t(new.mat))
    }
    else{ assign("old",new.mat)}
  }

    }
    }
    return(new.mat)}

   #return(old) }
    #####################

  #NEW<-IBM2(run=20,initial=1,path="YEPPERO",wantfull="Y",filePOPname="popstats2.csv",fileALLname="alldat",u.dispersal=u.dispersal,LDD=LDD,bf0=bf0,bf1=bf1,bf2=bf2,
  #  ############################SEED ESTABLISHMENT
  #  mu.est=mu.est,alpha.est=alpha.est,dis.est=dis.est,beta.est=beta.est,
  #  ############################Initial seedling height
  #  mu.inht=mu.inht,scale.inht=scale.inht,shape.inht=shape.inht,
  #    ############################seedling survival
  #  mu.ss=mu.ss,beta.ss=beta.ss,size.ss=size.ss,
  #    ############################seedling growth
  #  alpha.sg=alpha.sg,dis.sg=dis.sg,mu.sg=mu.sg,size.sg=size.sg,beta.sg=beta.sg,scale.sg=scale.sg,shape.sg=shape.sg,
  #      ############################switch from seedlings to adults
  #  aht=aht,bht=bht,sigma=sigma,a.switch=a.switch,b.switch=b.switch,
  #  ############################adult survival
  #    alpha.as=alpha.as,dis.as=exp(dis.as),P.s=P.s,G.s=G.s,
  #    ############################adult growth
  #     alpha.ag=alpha.ag,dis.ag=dis.ag,P.g=P.g,G.g=G.g,a.sig=a.sig,b.sig=b.sig,a.shape=a.shape,b.shape=b.shape)
