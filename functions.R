##%#########################################################################%##
#                                                                             #
####       Functions for comparing ROC and AROC curves (same marker)       ####
#                                                                             #
##%#########################################################################%##

#%########################################%#
## Authors: Arís Fanjul Hevia
## Purpose: This script contains the functions that are needed for testing whether
##        the ROC and the AROC curves of the same diagnostic marker are the same or not.
## Data: 11/2024
#%########################################%#

# Support functions -----------------------------------------------------------

#--- Epanechnikov's density function 
kernel=function(u){(0.75*(1-u^2))*(u<1)*(u>-1)}

#--- Epanechnikov's distribution function
Kerneldf=function(u){ifelse(u > 1,1, 0.25*(3*u-u^3+2)*(u>-1))} 

#--- Nadaraya-Watson estimator
NW=function(x, X, Y,h=hbw.cv(X,Y)){
  # x=points, X=X data, Y= Ydata, h=bw ( kernel used is Epa)
  matk=outer(x,X, function(x,X,h){kernel((x-X)/h )},h=h)
  drop((matk%*%Y) / (matk%*%rep(1,length(X))))
}
#--- Bandwidth selector - Cross Validation
hbw.cv <- function(X, Y, h.grid = diff(range(X)) * (seq(0.1, 1, l = 100))^2) {
  obj <- sapply(h.grid, function(h){
    return(sum(((Y - NW(x = X, X = X, Y = Y, h = h))/(1 - kernel(0) / colSums(kernel(outer(X, X, "-") / h)) ) )^2))
  })
  h <- h.grid[which.min(obj)]
  return(h)
}

# Graphics -----------------------------------------------------------
#--- Plotting conditional densities of F and G 
plot.densities.FG=function(XF,YF,XG,YG,x.name="Covariate",y.name="",main="",nxmar=15,col1=viridis(5,alpha=0.5)[2], col2=viridis(4,alpha=0.5)[3],font=1){
  cdeF=cde(XF,YF,x.name=x.name, y.name="y.name",x.margin = seq(min(XF,XG),max(XF,XG),length.out = nxmar),nxmargin = nxmar)
  cdeG=cde(XG,YG,x.name=x.name, y.name="y.name",x.margin = seq(min(XF,XG),max(XF,XG),length.out = nxmar),nxmargin = nxmar)
  
  ylab <- y.name
  xlab<-  x.name
  
  oldpar <- par()
  par(mar=c(3,0,0,1)+.2)
  maxden <- max(cdeF$z,cdeG$z,na.rm=TRUE)
  zlim <- c(0,maxden)
  xlim <- range(unlist(cdeF$x),unlist(cdeG$x),na.rm=TRUE)
  ylim <- range(cdeF$y,cdeG$y,na.rm=TRUE)
  xrange <- xlim[2]-xlim[1]
  yrange <- ylim[2]-ylim[1]
  junk <- persp(xlim,ylim,cbind(rep(0,2),rep(maxden,2)), zlim=zlim,box=FALSE,axes=FALSE,
                theta=-75,phi=25,border=NA)
  mtext(main,3,cex=1,font=font,line=-4)
  
  # Function to enable things to be added to perspective plot
  perspp <- function(x,y,z, pmat)
  {
    tr <- cbind(x,y,z,1) %*% pmat
    list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])
  }
  
  # add axes
  lines(perspp(xlim,c(ylim[1],ylim[1]),c(0,0),junk))
  lines(perspp(c(xlim[1],xlim[1]),ylim,c(0,0),junk))
  lines(perspp(xlim,c(ylim[2],ylim[2]),c(0,0),junk))
  lines(perspp(c(xlim[2],xlim[2]),ylim,c(0,0),junk))
  # add tick marks and scale
  ylabels <- pretty(ylim)
  ylabels <- ylabels[ylabels<=max(ylim) & ylabels >= min(ylim)]
  for(i in ylabels){
    lines(perspp(c(xlim[1],xlim[1] - xrange/75),c(i,i),c(0,0),junk))
    lines(perspp(c(xlim[1],xlim[1] + xrange/75),c(i,i),c(0,0),junk))
    text(perspp(xlim[1]-xrange/20,i,0,junk),paste(round(i,4)),cex=0.8,font=font)
  }
  xlabels <- pretty(xlim)
  xlabels <- xlabels[xlabels<=max(xlim) & xlabels >= min(xlim)]
  for(i in xlabels){
    lines(perspp(c(i,i),c(ylim[1],ylim[1]-yrange/75),c(0,0),junk))
    lines(perspp(c(i,i),c(ylim[2],ylim[2]+yrange/75),c(0,0),junk))
    text(perspp(i,ylim[1]-yrange/20,0,junk),paste(round(i,4)),cex=0.8,font=font)
  }
  midx <- 0.5*(xlim[2]+xlim[1])
  midy <- 0.5*(ylim[2]+ylim[1])
  text(perspp(midx,ylim[1] - yrange/7, 0,junk),xlab,cex=0.8,srt=265,font=font)
  text(perspp(xlim[1] - xrange/7,midy,0,junk),ylab, cex=0.8,srt=352,font=font)
  
  #n <- length(cbind(cde1$x,cde2$x))
  n1 <- length(cdeF$x)
  n2<-length(cdeG$x)
  n=max(n1,n2)
  
  for(i in n:1){
    if(n1>=i){
      pol <- perspp(rep(cdeF$x[i],length(cdeF$y)+2),
                    c(cdeF$y[1],cdeF$y,cdeF$y[length(cdeF$y)]),
                    c(0,.65*cdeF$z[i,],0),junk)
      polygon(pol$x, pol$y, col=col1, border=TRUE)
    }
    if(n2>=i){
      pol <- perspp(rep(cdeG$x[i],length(cdeG$y)+2),
                    c(cdeG$y[1],cdeG$y,cdeG$y[length(cdeG$y)]),
                    c(0,.65*cdeG$z[i,],0),junk)
      polygon(pol$x, pol$y, col=col2, border=TRUE)
    }
    
  }
  
  par(mar=oldpar$mar)
  par(new=FALSE)
}

# Main function ----------------------------------------------------------------
Test.ROCAROC=function(sampleY, sampleX, nboots=200, nr = 2){
## Arguments:
#.  - sampleY: list(YF,YG)
#.  - sampleX: list(XF,XG)
#.  - nboots: number of bootstrap iterations to consider for the test.
#.  - nr: 1/nr is the proportion of the sample that will be used for estimating 
#.      the ROC curve (the rest is used for estimating the AROC curve).

  XF=sampleX[[1]]
  XG=sampleX[[2]]
  YF=sampleY[[1]]
  YG=sampleY[[2]]
  
  n=length(YF)
  m=length(YG)
  resolution=n+m
  
  # Division of the sample **********************
  idF.roc=sample(1:n,round(n/nr),replace = FALSE)
  idG.roc=sample(1:m,round(m/nr),replace = FALSE)
  
  YF.roc=sampleY[[1]][idF.roc]
  YG.roc=sampleY[[2]][idG.roc]
  n.roc=length(YF.roc)
  m.roc=length(YG.roc)
  
  XF.aroc=sampleX[[1]][-idF.roc]
  XG.aroc=sampleX[[2]][-idG.roc]
  YF.aroc=sampleY[[1]][-idF.roc]
  YG.aroc=sampleY[[2]][-idG.roc]
  n.aroc=length(YF.aroc)
  m.aroc=length(YG.aroc)
  
  # AROC **************************
  # bandwiths:
  gF = hbw.cv(XF.aroc,YF.aroc)
  gG = hbw.cv(XG.aroc,YG.aroc)
  
  # mu and sigma
  muF.est = NW(XF.aroc,XF.aroc,YF.aroc,gF)
  muG.est = NW(XG.aroc,XG.aroc,YG.aroc,gG)
  sigmaF.est = sqrt(NW(XF.aroc,XF.aroc,(YF.aroc-muF.est)^2,gF))
  sigmaG.est = sqrt(NW(XG.aroc,XG.aroc,(YG.aroc-muG.est)^2,gG))
  
  # residuals:
  residF= (YF.aroc-muF.est)/sigmaF.est
  residG= (YG.aroc-muG.est)/sigmaG.est
  
  # AROC curve estimate:
  p=seq(0.001,0.999,length=resolution)
  muG.est.XF = NW(XF.aroc,XG.aroc,YG.aroc,gG)
  sigmaG.est.XF = sqrt(NW(XF.aroc,XG.aroc,(YG.aroc-muG.est)^2,gG))
  sss=sort(residG)[floor(m.aroc*(1-p))+1]
  yyy=(YF.aroc-muG.est.XF)/sigmaG.est.XF
  AROC0=rowMeans(outer(sss,yyy,FUN = "<"))
  
  # ROC ******************************
  ROC0= 1- ecdf(YF.roc)(sort(YG.roc)[floor(m.roc*(1-p))+1])
  
  # Statistics ***********************
  T1.0 = mean(abs(ROC0-AROC0 )) 
  T2.0 = mean((ROC0-AROC0 )^2) 
  TKS.0= max(abs(ROC0-AROC0))
  
  # Bootstrap loop ***********************
  T1.b=T2.b=TKS.b=numeric(nboots)
  gF.boots = gF
  gG.boots = gG
  
  pb=txtProgressBar(style=3) # progress bar
   
  for(b in 1:nboots){
    
    #AROC bootstrap 
    ## New sample
    residF.boots = sample(residF,replace = TRUE)
    residG.boots = sample(residG,replace = TRUE)
    
    YF.b.aroc= muF.est + sigmaF.est*residF.boots
    YG.b.aroc= muG.est + sigmaG.est*residG.boots  
    
    # muF.est.boots = NW(XF.aroc,XF.aroc,YF.b.aroc,gF.boots)
    muG.est.boots = NW(XG.aroc,XG.aroc,YG.b.aroc,gG.boots)
    # sigmaF.est.boots = sqrt(NW(XF.aroc,XF.aroc,(YF.b.aroc-muF.est.boots)^2,gF.boots))
    sigmaG.est.boots = sqrt(NW(XG.aroc,XG.aroc,(YG.b.aroc-muG.est.boots)^2,gG.boots))
    
    # residF.boots = (YF.b.aroc-muF.est.boots)/sigmaF.est.boots 
    residG.boots = (YG.b.aroc-muG.est.boots)/sigmaG.est.boots 
    
    muG.est.XF.boots = NW(XF.aroc,XG.aroc,YG.b.aroc,gG.boots)
    sigmaG.est.XF.boots = sqrt(NW(XF.aroc,XG.aroc,(YG.b.aroc-muG.est.boots)^2,gG.boots))
    sss.b=sort(residG.boots)[floor(m.aroc*(1-p))+1]
    yyy.b=(YF.b.aroc-muG.est.XF.boots)/sigmaG.est.XF.boots
    AROCb=rowMeans(outer(sss.b,yyy.b,FUN = "<"))
    
    # ROC bootstrap 
    
    YF.b.roc = sample(YF.roc,replace = TRUE)
    YG.b.roc = sample(YG.roc,replace = TRUE)
    ROCb= 1- ecdf(YF.b.roc)(sort(YG.b.roc)[floor(m.roc*(1-p))+1])
    
    # Boots statistics 
    T1.b[b] = mean(abs(ROCb-ROC0+AROC0-AROCb ))
    T2.b[b] = mean((ROCb-ROC0+AROC0-AROCb)^2)
    TKS.b[b] = max(abs(ROCb-ROC0+AROC0-AROCb))
    
    # progress bar
    setTxtProgressBar(pb,b/(nboots))
  }
  
  ## Computing statistics and pvalues *********
  pv.1=mean(T1.b>T1.0)
  pv.2=mean(T2.b>T2.0)
  pv.KS=mean(TKS.b>TKS.0)
  pval=list(pv.1,pv.2,pv.KS)
  names(pval)=c("pvalue T1","pvalue T2","pvalue TKS")
  
  stat=list(T1.0,T2.0,TKS.0)
  names(stat)=c("T1","T2","TKS")
  
  ## Output ******************
  output=list()             # The output is a list
  output$ROC=ROC0           # First item: estimation of ROC (using 1/nr of the sample)
  output$AROC=AROC0         # Second item: estimation of AROC (sign the rest of the sample) 
  output$p=p                # Third item: vector of values in [0,1] for which the previous curves are computed
  output$statistics = data.frame(stat)  # Fourth item: value of the statistics
  output$pvalues=data.frame(pval)       # Fifth item: pvalues for the 3 considered functions
  return(output)
}
