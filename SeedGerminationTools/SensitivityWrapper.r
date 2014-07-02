args=(commandArgs(TRUE))

for(i in 1:length(args)) {
	eval(parse(text=args[[i]])) 
# 	flag <- 1
# 	iter <- 1
}

library("abind")
library("lattice")
# library("sn")
# library("spam")
library('mnormt', lib="~/Rpckgs")
library('spam', lib="~/Rpckgs")
library('sn', lib='~/Rpckgs')
#setting spam up to give a larger nearestdistnnz so it doesn't crash:
# spam.options(nearestdistnnz=   c(1e9,400))
spam.options(nearestdistnnz=c(53248651, 400))
  
set.seed(iter)
source('IBM_2.25.r')

################THIS READS IN THE PARAMETER VALUES FROM THE STATISTICAL ANALYSIS  
# IBMpars<-read.csv("IBMpars11.27.csv",header=T)
IBMpars<-read.csv("IBMpars12.16.csv",header=T)  
###############FOR THE TEST RUNS, I have been just looking at the median parameter values
#     med.par<-lapply(IBMpars,median)
#     attach(med.par)
#     draw.par 				<- quantile(IBMpars[,flag], probs=iter/100)
# 	med.par				<- lapply(IBMpars, median)
# 	med.par[[flag]]	<- draw.par
#     attach(med.par)
	
ran.iters 				<- sample(1:(dim(IBMpars)[1]), size=dim(IBMpars)[2], replace=TRUE) 
IBMdraw 					<- diag(as.matrix(IBMpars[ran.iters,]))
names(IBMdraw)	  <- dimnames(IBMpars)[[2]]
IBMdraw 					<- as.data.frame(t(IBMdraw))

# IBMdraw  					  <- NULL
# IBMdraw$a.switch	  <- median(IBMpars$a.switch)
# IBMdraw$b.switch	  <- median(IBMpars$b.switch)
# IBMdraw$alpha.inht	<- median(IBMpars$alpha.inht)

attach(IBMdraw)

initial.state					<- sample(1:8, 1)
names(initial.state)	<- "initial"
draw.file <- "output/ParDraws.csv"

if (file.exists(draw.file)) {
	write.table(file="output/ParDraws.csv", t(c(iter, initial.state, unlist(IBMdraw))), append=T, row.names=F, col.names=F)
} else {
	write.table(file="output/ParDraws.csv", t(c(iter, initial.state, unlist(IBMdraw))), append=T, row.names=F, col.names=T)
}


popfile	<- paste('output/popfile', iter, '.csv', sep='')
allfile 	<- paste('output/allfile', iter, sep='')


##what does path do?
  Try1<-IBM2(run=100, initial=initial.state, wantfull="Y", filePOPname=popfile, fileALLname=allfile, u.dispersal=u.dispersal,LDD=LDD,bf0=bf0,bf1=bf1,bf2=bf2,
    ############################SEED ESTABLISHMENT
    mu.est=mu.est,alpha.est=alpha.est,dis.est=dis.est,beta.est=beta.est,
    ############################Initial seedling height
    mu.inht=mu.inht,scale.inht=scale.inht,shape.inht=shape.inht,
      ############################seedling survival
    mu.ss=mu.ss,beta.ss=beta.ss,size.ss=size.ss,
      ############################seedling growth
    alpha.sg=alpha.sg,dis.sg=dis.sg,mu.sg=mu.sg,size.sg=size.sg,beta.sg=beta.sg,scale.sg=scale.sg,shape.sg=shape.sg,
        ############################switch from seedlings to adults
    aht=aht,bht=bht,sigma=sigma,a.switch=a.switch,b.switch=b.switch,
    ############################adult survival
      alpha.as=alpha.as,dis.as=dis.as,P.s=P.s,G.s=G.s,
      ############################adult growth
       alpha.ag=alpha.ag,dis.ag=dis.ag,P.g=P.g,G.g=G.g,a.sig=a.sig,b.sig=b.sig,a.shape=a.shape,b.shape=b.shape)

# Try1	<-  IBM2(run=100, initial=initial.state, wantfull="Y", filePOPname=popfile, fileALLname=allfile, u.dispersal=u.dispersal,LDD=LDD,bf0=bf0,bf1=bf1,bf2=bf2, mu.est=mu.est, alpha.est=alpha.est, dis.est=dis.est, beta.est=beta.est, mu.inht=mu.inht, alpha.inht=alpha.inht, dis.inht=dis.inht, scale.inht=scale.inht, shape.inht=shape.inht, mu.ss=mu.ss, beta.ss=beta.ss, size.ss=size.ss, alpha.sg=alpha.sg, dis.sg=dis.sg, mu.sg=mu.sg, size.sg=size.sg, beta.sg=beta.sg, scale.sg=scale.sg, shape.sg=shape.sg, aht=aht,bht=bht, sigma=sigma, a.switch=a.switch, b.switch=b.switch, alpha.as=alpha.as, dis.as=dis.as, P.s=P.s, G.s=G.s, alpha.ag=alpha.ag, dis.ag=dis.ag, P.g=P.g, G.g=G.g, a.sig=a.sig, b.sig=b.sig, a.shape=a.shape, b.shape=b.shape)

detach(IBMdraw)

q('no')
	     
