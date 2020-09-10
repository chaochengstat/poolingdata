####################################################################################
#                                                                                  #
#                                R Code for                                        #
#     Statistical Methods for Analysis of Combined Biomarker Data from Multiple    #
#                        Nested Case-Control Studies (continuous biomarker)        #
#                                                                                  #
#                                                                                  #
#              For any question, feel free to contact Chao Cheng                   #                                                                                #    
#                         Email: c.cheng@yale.edu                                  #
#                                                                                  #
#                                                                                  #
####################################################################################


#################################################################################################
#
# There are two sections, the 1st section is the functions, the 2nd is an illustrative example
# Note: the case-control ratio for all matched set must be either 1:2 or 1:1.
#
#################################################################################################

#############################################################################################
#   
# SECTION 1: FUNCTIONS
#
#############################################################################################


##########################################################################
##########################################################################
# calculate the gradient
##########################################################################
# INPUT
# f: traget function
# x0: point where the gradient is to build
# heps: step size
##########################################################################
# RETURN
# the numerical gradient of target function
##########################################################################

mygrad=function (f, x0,heps = 1e-6, ...) {
  if (!is.numeric(x0)) 
    stop("Argument 'x0' must be a numeric value.")
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  p =length(f(x0))
  n <- length(x0)
  hh <- rep(0, n)
  gr <- matrix(0,nrow=n,ncol=p)
  for (i in 1:n) {
    hh[i] <- heps
    gr[i,] <- (f(x0 + hh) - f(x0 - hh))/(2 * heps)
    hh[i] <- 0
  }
  return(gr)
}

##########################################################################
##########################################################################
# logit function: f(x) = e^x/(1+e^x)
##########################################################################
# INPUT
# x: point where we need to calculate its logit value
##########################################################################
# RETURN
# the logit value of x
##########################################################################

logit=function(x) {exp(x)/(1+exp(x))}

##########################################################################
##########################################################################
# Gauss-Hermite Quadrature Functions (From GHQ R package)
##########################################################################

hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  r <- cbind(x * sqrt(2), w/sum(w))
  colnames(r) <- c("Points", "Weights")
  r
}

mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}

weight.f=function(n=5,mean,sigma) {
  #n=5;mean=rnorm(10);sigma=seq(1,10,by=1)
  N=length(mean)
  loc=matrix(0,ncol=n,nrow=N)
  gaus=statmod::gauss.quad.prob(n=n,dist="normal")
  x=gaus$node
  w=gaus$weights
  for (i in (1:n)) {
    loc[,i] = mean + sqrt(sigma) * x[i]
  }
  list(loc=loc,weights=w)
}


##############################################################################################
##############################################################################################
# Data generating functions (based on the simulation studies section in the main paper)
##############################################################################################

data.sim.match=function(OR=1.25,beta2=0,type=2,matchratio=2,errorType="normal",n1=200) {
  N = 15000
  M=10
  var_x = 1
  var_l = runif(M+1,0.15,0.35)
  alpha0_l = rnorm(M,3,sd=1)
  alpha1_l = 0.5
  gam=rnorm(M+1,0,sd=0.1)
  epsi = rnorm(M+1,mean=0,sd=0.5)
  mat=matrix(0,ncol = 5,nrow = M*N)
  colnames(mat)=c("W","X","HC","HL","group")
  for (i in (1:M)) {
    alpha0.now=alpha0_l[i]
    W   = rnorm(N,mean=0,sd=sqrt(2))
    if (errorType=="normal") {
      X   = rnorm(N,mean = alpha0.now + alpha1_l*W ,sd=sqrt(var_x))
    } else if(errorType=="uniform") {
      X   = alpha0.now+alpha1_l*W + runif(N,min=-sqrt(3*var_x),max=sqrt(3*var_x))
    } else {
      X   = alpha0.now+alpha1_l*W + fGarch::rsnorm(N, mean = 0, sd = sqrt(var_x), xi = 1.5)
    }
    H_C = epsi[1]+(gam[1]+1)*rnorm(N,mean=X,sd=sqrt(var_l[1]))
    H_L = epsi[i+1]+(gam[i+1]+1)*rnorm(N,mean=X,sd=sqrt(var_l[i+1]))
    mat[((i-1)*N+1):(i*N),1] = W
    mat[((i-1)*N+1):(i*N),2] = X
    mat[((i-1)*N+1):(i*N),3] = H_C
    mat[((i-1)*N+1):(i*N),4] = H_L
    mat[((i-1)*N+1):(i*N),5] = i
  }
  mu.beta=matrix(0,ncol=2,nrow=M)
  for (i in (1:M)) {mu.beta[i,] = c(0,alpha0_l[i] + alpha1_l * 0) }
  sigma=matrix(c(5,alpha1_l*5,alpha1_l*5,alpha1_l*alpha1_l*5+var_x),ncol=2,nrow=2,byrow=T)
  beta0=rep(0,M)
  prob.vector=rep(0,M*N)
  beta0= beta0.f(prevalence=0.5,beta.x=log(OR),mu=mean(mat[,2]),sd=sd(mat[,2])) 
  for (i in (1:M)) {
    prob.vector[((i-1)*N+1):(i*N)] = logit(
      beta0 + log(OR)* mat[((i-1)*N+1):(i*N),2] + beta2*mat[((i-1)*N+1):(i*N),1])
  }
  y <- as.numeric(runif(M*N) < prob.vector)
  casenum = as.integer(n1/(1+matchratio))
  num=casenum*(1+matchratio)
  data=matrix(0,ncol = dim(mat)[2]+2,nrow = M*num)
  for (i in (1:M)) {
    y.now   = y[((i-1)*N+1):(i*N)]
    mat.now = mat[((i-1)*N+1):(i*N),]
    loc1=sample(which(y.now==1),casenum,replace=F)
    loc2=sample(which(y.now==0),num-casenum,replace=F)
    id1=i*1000+1:length(loc1);id=rep(id1,matchratio+1)
    data[((i-1)*num+1):(i*num),-c(1:2)] = mat.now[c(loc1,loc2),]
    data[((i-1)*num+1):(i*num),1:2] = cbind(y.now[c(loc1,loc2)],id)
  }
  n.cali= as.integer(n1/10)
  is.calibration=rep(0,M*num)
  for (i in (1:M)) {
    if (type==1) {
      a = sample(((i-1)*num+1):(i*num),n.cali,replace=F)
      is.calibration[a]=1
    } else {
      pool=intersect(which(data[,1]==0), ((i-1)*num+1):(i*num))
      a = sample(pool,n.cali,replace=F)
      is.calibration[a]=1
    }
  }
  data=cbind(data,is.calibration)
  colnames(data)=c("Y","matchID","W","X","HC","HL","group","calibration")
  list(data=data,sigma=c(var_x,var_l),beta=c(beta0,log(OR),beta2),theta=c(alpha0_l,alpha1_l,epsi,gam),epsi=epsi,gam=gam)
}

beta0.f=function(prevalence=0.05,beta.x=0.2,mu=0,sd=1) {
  #prevalence=0.05;beta.x=0.2;mu=0;sd=1
  #log(prevalence/(1-prevalence)) - beta.x * mu
  cc <- pracma::gaussHermite(5)
  myf=function(a) {
    sum(cc$w*exp(a + beta.x*1.414*sd*cc$x+beta.x*mu)/(1+exp(a + beta.x*1.414*sd*cc$x+beta.x*mu))*1/sqrt(3.1415))-prevalence
  }
  tryCatch(uniroot(myf,interval=c(-100,100))$root,
           error=function(e) {log(prevalence/(1-prevalence)) - beta.x * mu})
}

##############################################################################################
##############################################################################################
# Data transferring step: transferring the original dataset to the datasets that can be used in 
# the paramater estimation steps
##############################################################################################

transdata=function(data) {
  n=dim(data)[1]
  data=cbind(data,loc=1:n)
  y = data[,-which(colnames(data)=="HC")]
  y = cbind(y,lab=y[,"group"])
  ystar = data[which(data[,"calibration"]==1),][,-which(colnames(data)=="HL")] 
  ystar = cbind(ystar,lab=0)
  colnames(ystar)[which(colnames(ystar)=="HC")]="HL"
  res=rbind(y,ystar)
}

data_lme=function(data=data,Y="Y",HL="HL",HC="HC",Wname="W",ID="matchID",group="group",calibration="calibration") {
  WZ = setdiff(colnames(data),c(Y,ID,HC,HL,group,calibration))
  data=cbind(data[,c(Y,ID,HC,HL,group,calibration)],data[,WZ])
  colnames(data)[1:6] = c("Y","matchID","HC","HL","group","calibration")
  
  mydata = transdata(data)
  mydata = as.data.frame(mydata)
  Du     = relevel(as.factor(mydata[,group]),ref="1")
  if (is.null(Wname)) {
    fu = as.formula(paste(HL,"~","Du",sep=""))
    fit    = fitted(lm(fu,data=mydata))
  } else {
    Wfor = paste(Wname,sep="",collapse ="+")
    fu = as.formula(paste(HL,"~","Du+",Wfor,sep=""))
    fit    = fitted(lm(fu,data=mydata))
  }
  
  list(data=cbind(mydata,Du,fit),data0=data)
}

##############################################################################################
##############################################################################################
# Calibration parameters estimation function
##############################################################################################

lme_mydata1=function(mydata,type=1,gam=rep(0,6),Wname="W") {
  M=max(mydata$lab)
  mydata$gepsi = 1
  if (type==1) {
    for (m in (0:M)) {
      mylm=lm(HL~fit,data=mydata[which(mydata$lab==m),])
      mydata$gepsi[which(mydata$lab==m)]= 1 
    }
  }
  if (type==2) {
    for (m in (1:M)) {
      mydata$gepsi[which(mydata$calibration==0 & mydata$group==m)]=1+gam[m+1]
      mydata$gepsi[which(mydata$calibration==1 & mydata$group==m)]=1+gam[1]
    }
  }
  if (is.null(Wname)) {
    form <- HL ~ Du  + (1 | lab) + (gepsi-1 | loc) + (fit-1|lab)
  } else {
    form <- as.formula(paste("HL~Du  + (1 | lab) + (gepsi-1 | loc) + (fit-1|lab)",paste(Wname,collapse ="+"),sep="+"))
  }
  lf <- lme4::lFormula(formula = form, data = mydata, REML = FALSE, weights = rep(1,dim(mydata)[1]))
  mylmer <- function(phi) {
    u=phi[1:M]
    n=dim(mydata)[1]
    myweight=function(w){
      ret=rep(1,n)
      for (i in (1:M)) {
        ret[which(mydata$lab==i)] = w[i]
      }
      ret
    }
    wei=myweight(u)
    lf$fr[,"(weights)"]=wei
    devf <- do.call(lme4::mkLmerDevfun, lf)
    devf(phi[-c(1:M)])
  }
  z = nlminb(start=c(rep(1,M),1,1,1),objective=mylmer,lower=c(rep(0.01,M),-Inf,-Inf,-Inf,-Inf),upper=c(rep(Inf,M),Inf,Inf,Inf,Inf))
  myweight=function(w){
    ret=rep(1,dim(mydata)[1])
    for (i in (1:M)) {
      ret[which(mydata$lab==i)] = w[i]
    }
    ret
  }
  wei=myweight(z$par[1:M])
  if (is.null(Wname)) {
    m2=lme4::lmer(HL ~ Du + (1| lab)+(gepsi-1|loc)+(fit-1|lab) ,data = mydata, REML = F, 
                  weights = wei,control=lme4::lmerControl(check.conv.singular = .makeCC(action = "message",  tol = 0)))
  } else {
    form <- as.formula(paste("HL~Du  + (1 | lab) + (gepsi-1 | loc) + (fit-1|lab)",paste(Wname,collapse ="+"),sep="+"))
    m2=lme4::lmer(form ,data = mydata, REML = F, 
                  weights = wei,control=lme4::lmerControl(check.conv.singular = .makeCC(action = "message",  tol = 0)))
  }
  
  
  hessinfo=solve(pracma::hessian(mylmer,x0=z$par))
  
  list(model=m2,parinfo=z$par,hessinfo=hessinfo)
}

lme_mydata2=function(Wname,Zname) {
  if (is.null(Wname)) {
    form <- as.formula(paste("HL~Du"))
  } else {
    form <- as.formula(paste("HL~Du",paste(Wname,collapse ="+"),sep="+"))
  }
  mydata$data[,"fit"]=fitted(lm(form,data=mydata$data))
  mydata<<-list(data=mydata$data,data0=mydata$data0)
  options(warn=-1)
  
  i=1
  stop.cert=1
  while( stop.cert>0.001 & i<=100) {
    if (i==1) {
      m1<<-lme_mydata1(mydata$data,type=1,Wname=Wname)
    } else {
      m1<<-lme_mydata1(mydata$data,type=2,gam=gam,Wname=Wname)
    }
    mydata$data[,"fit"]=as.vector(model.matrix(form,data=mydata$data) %*% fixef(m1$model))
    xigma.mat = nlme::random.effects(m1$model,condVar = FALSE)$lab
    xigma = c(xigma.mat[,2],xigma.mat[,1]+1)
    M=dim(xigma.mat)[1]-1
    gam=xigma[(M+2):(2*M+2)]-1
    if (i==1) {
      loglike0=summary(m1$model)$logLik[1]+10
      loglike1=summary(m1$model)$logLik[1]
    } else {
      loglike1=summary(m1$model)$logLik[1]
    }
    stop.cert=abs(loglike1-loglike0)
    i=i+1
    loglike0=loglike1
  }
  
  return(m1)
}


##############################################################################################
##############################################################################################
# OR estimation function
##############################################################################################

get.beta=function(mat=mydata$data0$data,model=m1,Wname="W",Zname=NULL) {
  M=max(mat[,"group"])
  data=mat
  m1=model
  xigma.mat = nlme::random.effects(m1,condVar = FALSE)$lab
  xigma = c(xigma.mat[,2],xigma.mat[,1]+1)
  fixed = nlme::fixed.effects(m1)
  alpha = fixed[1]+c(0,fixed[2:M])
  altau = c(alpha,fixed[-c(1:M)])
  N=dim(mat)[1]
  
  rsigma=as.data.frame(summary(m1)$varcor)[,4][c(1,4)]
  w.data=m1@frame
  weigh=rep(1,M+1) 
  for (i in (1:M)) {
    weigh[i+1]=w.data[which(w.data[,"lab"]==i),][1,"(weights)"]
  }
  sigma= c(rsigma[1],abs(rsigma[2] * (1/weigh) ))
  
  M=length(sigma)-2
  sigma=abs(sigma)
  xi=xigma[1:(M+1)]
  gam=xigma[-c(1:(M+1))]
  alpha=altau[1:M]
  tau=altau[-c(1:M)]
  N = dim(mat)[1]
  
  sigmar=sigma[-c(1:2)]/(gam[-1]^2)
  sigma0=sigma[2]/(gam[1]^2)
  w1 = sigma[1]/(sigma[1] + sigmar)
  w2 = sigma[1]/(sigma[1]+sigma0*sigmar/(sigma0+sigmar))
  w2.in = sigmar/(sigma0+sigmar)
  
  sigma.hat.X1 = sigma[1]*sigmar/(sigma[1] + sigmar)
  sigma.hat.X2 = w2 * sigma[2]*sigmar/(sigma0+sigmar)
  
  Xe = H.mean = rep(0,N)
  sigma.vec=rep(0,N)
  for (i in (1:M)) {
    loc1 = which(data[,"group"]==i & data[,"calibration"]==0)
    loc2 = which(data[,"group"]==i & data[,"calibration"]==1)
    if (is.null(Wname)) {
      Xe[loc1] = w1[i]*(data[loc1,"HL"] - xi[i+1])/ gam[i+1] + (1-w1[i])*( alpha[i])
      Xe[loc2] = w2[i]*(  (1-w2.in[i])*(data[loc2,"HC"]-xi[1])/gam[1] + w2.in[i]*(data[loc2,"HL"]-xi[i+1])/gam[i+1]  ) + (1-w2[i])*(alpha[i])
    } else {
      Xe[loc1] = w1[i]*(data[loc1,"HL"] - xi[i+1])/ gam[i+1] + (1-w1[i])*( alpha[i] + data[loc1,Wname] %*% matrix(tau,ncol=1))
      Xe[loc2] = w2[i]*(  (1-w2.in[i])*(data[loc2,"HC"]-xi[1])/gam[1] + w2.in[i]*(data[loc2,"HL"]-xi[i+1])/gam[i+1]  ) + (1-w2[i])*(alpha[i] + data[loc2,Wname] %*% matrix(tau,ncol=1))
    }
    
    H.mean[loc1] = data[loc1,"HL"]
    H.mean[loc2] = (data[loc2,"HL"] + data[loc2,"HC"])/2
    sigma.vec[loc1] = sigma.hat.X1[i]
    sigma.vec[loc2] = sigma.hat.X2[i]
  }
  
  data = cbind(data,Xe,H.mean,sigma.vec)
  id=unique(data[,"matchID"])
  mergedata=function(id,nam="H.mean") {
    a=data[which((data[,"Y"]==1) & (data[,"matchID"]==id)),nam]
    b=data[which((data[,"Y"]==0) & (data[,"matchID"]==id)),nam]
    c(a,b)
  }
  data.n=sapply(id, mergedata, nam="H.mean", simplify = FALSE)
  data.c=sapply(id, mergedata, nam="Xe", simplify = FALSE)
  data.sigma=sapply(id, mergedata, nam="sigma.vec", simplify = FALSE)
  ###### naive and approximate method
  if (is.null(Zname)) {
    condloglike=function(par0,data=data.n) {
      -sum(log(base::unlist(lapply(data,function(x) { 1/(1+sum(exp(par0*(x[-1]-x[1]))))}))))
    }
  } else {
    data.z1=sapply(id, mergedata, nam=Zname[1], simplify = F)
    condloglike=function(par0,data=data.n) {
      -sum(log(base::unlist(mapply(function(x,y) { 1/(1+sum(exp(par0[1]*(x[-1]-x[1]) + par0[2]*(y[-1]-y[1])  )))},data,data.z1))))
    }
  }
  if (is.null(Zname)) {
    par.n=nlminb(start=1, objective=condloglike,data=data.n)$par
    sd.n = as.vector(sqrt(1/pracma::hessian(condloglike,par.n)))
    par.c=nlminb(start=par.n, objective=condloglike,data=data.c)$par
    sd.c=as.vector(sqrt(1/pracma::hessian(condloglike,par.c)))
  } else {
    par.n=nlminb(start=rep(1,2), objective=condloglike,data=data.n)$par
    sd.n = as.vector(sqrt(diag(1/pracma::hessian(condloglike,par.n))))
    par.c=nlminb(start=par.n, objective=condloglike,data=data.c)$par
    sd.c=as.vector(sqrt(diag(1/pracma::hessian(condloglike,par.c))))
  }
  
  
  ## exact calibration method
  # (1) Monte Carlo Method
  seed.mat=matrix(rnorm(25*dim(data)[1]),ncol=25,nrow=dim(data)[1])
  X.mat.func=function(a) {
    #a=id[1]
    x=data.c[[which(id==a)]]
    y=data.sigma[[which(id==a)]]
    z=seed.mat[which(data[,"matchID"]==a),]
    x %*% t(rep(1,25)) + sqrt(y %*% t(rep(1,25)))*z
  }
  X.mat=sapply(id,X.mat.func,simplify = F)
  if (is.null(Zname)) {
    loglike.e.mc=function(par) {
      -sum( log(mapply(function(x){mean(exp(par*x[1,])/apply(exp(par*x),2,sum))},X.mat)) )
    }
  } else {
    loglike.e.mc=function(par) {
      -sum( log(mapply(function(x,y){mean(exp(par[1]*x[1,]+rep(par[2]*y[1],25))/apply(exp(par[1]*x+(par[2]*y) %*% t(rep(1,25))),2,sum))},X.mat,data.z1)) )
    }
  }
  par.e1=nlminb(start=par.c, objective=loglike.e.mc)$par
  sd.e1=as.vector(sqrt(diag(1/pracma::hessian(loglike.e.mc,par.e1))))
  ## Gauss Hermite Method
  ratio=base::unlist(lapply(data.c,length))
  ratio2= ratio==2;ratio3= ratio==3
  mysubset=function(x,y){
    if(sum(y)==0) {
      NULL
    } else {
      t(sapply(subset(x,y),function(x) x,simplify=T))
    }
  }
  data.c2=mysubset(data.c,ratio2);data.c3=mysubset(data.c,ratio3)
  if (is.null(Zname)==0) {
    data.z12=mysubset(data.z1,ratio2);data.z13=mysubset(data.z1,ratio3)
  }
  
  data.sigma2=mysubset(data.sigma,ratio2);data.sigma3=mysubset(data.sigma,ratio3)
  
  weights0=weight.f(n=10,mean=data.c2[,2]-data.c2[,1],sigma=data.sigma2[,2]+data.sigma2[,1])
  
  mean.list=cbind(data.c3[,2]-data.c3[,1],data.c3[,3]-data.c3[,1])
  if (dim(mean.list)[1] !=0) {
    sigma.list=cbind(data.sigma3[,2]+data.sigma3[,1],data.sigma3[,1],data.sigma3[,1],data.sigma3[,3]+data.sigma3[,1])
    p1=matrix(0,ncol=9,nrow=dim(mean.list)[1])
    p2=matrix(0,ncol=9,nrow=dim(mean.list)[1])
    we=matrix(0,ncol=9,nrow=dim(mean.list)[1])
    for (m in (1:dim(mean.list)[1])) {
      zz=mgauss.hermite(3, mu=mean.list[m,], sigma=matrix(sigma.list[m,],ncol=2))
      p1[m,]=zz$points[,1];p2[m,]=zz$points[,2];we[m,]=zz$weights
    }
  }
  
  
  if (is.null(Zname)) {
    loglike.e=function(par) {
      if (is.null(data.c3)) {
        xx=0
      } else {
        xx=apply(log( ( 1/(1+exp(par*p1)+exp(par*p2)) ) %*% apply(we,2,mean)),2,sum)
      }
      if (is.null(data.c2)) {
        yy=0
      } else {
        yy=apply(log( ( 1/(1+exp(par*weights0$loc)) ) %*% weights0$weights),2,sum)
      }
      -yy-xx
    }
  } else {
    loglike.e=function(par) {
      if (is.null(data.c2)) {
        xx=0
      } else {
        xx=apply(log( ( 1/(1+exp(par[1]*weights0$loc + (par[2]*(data.z12[,-1]-data.z12[,1])
        )  %*% t(rep(1,10)) )) ) %*% weights0$weights),2,sum)
      }
      if (is.null(data.c3)) {
        yy=0
      } else {
        yy=apply(log( ( 1/(1+exp(par[1]*p1+ (par[2]*(data.z13[,2]-data.z13[,1])
        ) %*% t(rep(1,9)))+exp(par[1]*p2+(par[2]*(data.z13[,2]-data.z13[,1])
        ) %*% t(rep(1,9)))) ) %*% apply(we,2,mean)),2,sum)
      }
      -yy-xx
    }
  }
  par.e=nlminb(start=par.c, objective=loglike.e)$par
  sd.e=as.vector(sqrt(diag(1/pracma::hessian(loglike.e,par.e))))
  
  #### full calibration method
  fc=rep(0,length(mat[,1]))
  for (i in (1:M)) {
    mylm=lm(HC~HL,data=as.data.frame(mat[which(mat[,"group"]==i & mat[,"calibration"]==1),]))
    fc[which(mat[,"group"]==i)]=predict.lm(mylm,newdata=data.frame(HL=mat[which(mat[,"group"]==i),"HL"]))
  }
  data = cbind(data,fc)
  data.fc=sapply(id, mergedata, nam="fc", simplify = F)
  if (is.null(Zname)) {
    par.fc =nlminb(start=1, objective=condloglike,data=data.fc)$par
  } else {
    par.fc =nlminb(start=rep(1,2), objective=condloglike,data=data.fc)$par
  }
  sd.fc  = as.vector(sqrt(diag(1/pracma::hessian(condloglike,par.fc))))
  
  list(par.fc=par.fc[1],par.c=par.c[1],par.n=par.n[1],par.e=par.e[1],par.e1=par.e1[1],sd.c=sd.c[1],sd.fc=sd.fc[1],sd.n=sd.n[1],sd.e=sd.e[1],sd.e1=sd.e1[1],seed.mat=seed.mat)
}

##############################################################################################
##############################################################################################
# Confidence Interval using a resampling approach
##############################################################################################

mycondVar=function(object) {
  condVar=TRUE
  ans <- object@pp$b(1)
  if (!is.null(object@flist)) {
    levs <- lapply(fl <- object@flist, levels)
    asgn <- attr(fl, "assign")
    cnms <- object@cnms
    nc <- lengths(cnms)
    nb <- diff(object@Gp)
    nbseq <- rep.int(seq_along(nb), nb)
    ml <- split(ans, nbseq)
    for (i in seq_along(ml)) ml[[i]] <- matrix(ml[[i]], ncol = nc[i], 
                                               byrow = TRUE, dimnames = list(NULL, cnms[[i]]))
    ans <- lapply(seq_along(fl), function(i) {
      m <- ml[asgn == i]
      b2 <- vapply(m, nrow, numeric(1))
      ub2 <- unique(b2)
      if (length(ub2) > 1) 
        stop("differing numbers of b per group")
      rnms <- if (ub2 == length(levs[[i]])) 
        levs[[i]]
      else seq(ub2)
      data.frame(do.call(cbind, m), row.names = rnms, check.names = FALSE)
    })
    names(ans) <- names(fl)
    whichel = names(ans)
    stopifnot(is(whichel, "character"))
    whchL <- names(ans) %in% whichel
    ans <- ans[whchL]
    if (condVar) {
      sigsqr <- lme4:::sigma.merMod(object)^2
      rp <- rePos$new(object)
      ############################################
      cv=lme4:::condVar(object, scaled = TRUE)
      rp <- rePos$new(object)
      trms <- rp$terms
      n <- diff(rp$offsets)
      cv2 <- lme4:::bdiag_to_mlist(cv, n)
      cv3=list()
      cv3[[1]]=lme4:::bdiag_to_array(cv2[[2]],1)
      cv3[[2]]=lme4:::bdiag_to_array(cv2[[3]],1)
      names(cv3) <- list(rp$cnms[[2]],rp$cnms[[3]])
      res <- list()
      tt <- trms[[2]]
      res[[1]] <- cv3[1:2]
      names(res) <- names(rp$flist)[2]
      vv=res
      ############################################
      attr(ans[[2]], "postVar") <- vv
      #ans=ans[[2]]
    }
    class(ans) <- "ranef.mer"
  }
  ans
}

optimwrap00 = function(fn,par=0,lower=rep(-Inf,3),upper=rep(Inf,3),control=list(pp=0),...) {
  res <- nlminb(start=control$pp, objective=fn, lower=lower,upper=upper,control=list(eval.max=0))
  with(res, list(par  = par,
                 fval = objective,
                 feval= evaluations[1],
                 conv = convergence,
                 message = message))
}
get.beta.var=function(OR=1.5,mat=mydata$data0$data,rep=15,model=m1,mymodel=m1,Wname="W",seed.mat,Zname=NULL) {
  #mat=mydata$data0;model=m1$model;Wname=Wname;Zname=Zname
  #OR=1.5;rep=15;mymodel=m1
  
  M=max(mat[,"group"])
  N=dim(mat)[1]
  data=mat
  allsigma.mat=abs(mvtnorm::rmvnorm(rep, mean=mymodel$parinfo,sigma=mymodel$hessinfo))
  m1=model
  xigma.mat = nlme::random.effects(m1,condVar = FALSE)$lab
  gam=xigma.mat[,1]
  mydata$data$gepsi=1
  for (m in (1:M)) {
    mydata$data$gepsi[which(mydata$data$calibration==0 & mydata$data$group==m)]=1+gam[m+1]
    mydata$data$gepsi[which(mydata$data$calibration==1 & mydata$data$group==m)]=1+gam[1]
  }
  
  res=matrix(0,ncol=3,nrow=rep)
  sd.mat=matrix(0,ncol=3,nrow=rep)
  optimwrap00 = function(fn,par=0,lower=rep(-Inf,3),upper=rep(Inf,3),control=list(pp=0),...) {
    res <- nlminb(start=control$pp, objective=fn, lower=lower,upper=upper,control=list(eval.max=0))
    with(res, list(par  = par,
                   fval = objective,
                   feval= evaluations[1],
                   conv = convergence,
                   message = message))
  }
  
  
  for (j in (1:rep)) {
    
    myweight=function(w){
      ret=rep(1,dim(mydata$data)[1])
      for (i in (1:M)) {
        ret[which(mydata$data$lab==i)] = w[i] 
      }
      ret
    }
    wei=myweight(allsigma.mat[j,c(1:M)])
    if (is.null(Wname)) {
      form <- as.formula("HL~Du  + (1 | lab) + (gepsi-1 | loc) + (fit-1|lab)")
      m2=lme4::lmer(form,data = mydata$data, REML = F, weights = wei,
                    control = lme4::lmerControl(optimizer = "optimwrap00",optCtrl=list(pp=allsigma.mat[j,-c(1:M)]),
                                                check.conv.singular = .makeCC(action = "message",  tol = 0)))
    } else {
      form <- as.formula(paste("HL~Du  + (1 | lab) + (gepsi-1 | loc) + (fit-1|lab)",paste(Wname,collapse ="+"),sep="+"))
      m2=lme4::lmer(form,data = mydata$data, REML = F, weights = wei,
                    control = lme4::lmerControl(optimizer = "optimwrap00",optCtrl=list(pp=allsigma.mat[j,-c(1:M)]),
                                                check.conv.singular = .makeCC(action = "message",  tol = 0)))
    }
    
    
    random.var=mycondVar(m2)[[2]]
    fix.var=summary(m2)$coefficients
    fix.mat=as.vector(mvtnorm::rmvnorm(1,mean=fix.var[,1],sigma=diag(fix.var[,2]^2)))
    random.simfit=as.vector(mvtnorm::rmvnorm(1,mean=random.var[,1],sigma=diag(attr(random.var,"postVar")$lab$fit[,,1:(M+1)])))
    random.simint=as.vector(mvtnorm::rmvnorm(1,mean=random.var[,2],sigma=diag(attr(random.var,"postVar")$lab$`(Intercept)`[,,1:(M+1)])))
    w.data=m2@frame
    weigh=rep(1,M+1) 
    for (i in (1:M)) {
      weigh[i+1]=w.data[which(w.data[,"lab"]==i),][1,"(weights)"]
    }
    xigma.mat=c(random.simint,random.simfit+1)
    altau.mat=fix.mat
    rsigma=as.data.frame(summary(m2)$varcor)[,4][c(1,4)]
    xigma=xigma.mat;altau=altau.mat
    sigma= c(rsigma[1],abs(rsigma[2] * (1/weigh) ))
    #sigma=mydata$data0$sigma
    M=length(sigma)-2
    sigma=abs(sigma)
    xi=xigma[1:(M+1)]
    gam=xigma[-c(1:(M+1))]
    alpha=altau[1:M]
    tau=altau[-c(1:M)]
    N = dim(mat)[1]
    
    sigmar=sigma[-c(1:2)]/(gam[-1]^2)
    sigma0=sigma[2]/(gam[1]^2)
    w1 = sigma[1]/(sigma[1] + sigmar)
    w2 = sigma[1]/(sigma[1]+sigma0*sigmar/(sigma0+sigmar))
    w2.in = sigmar/(sigma0+sigmar)
    
    sigma.hat.X1 = sigma[1]*sigmar/(sigma[1] + sigmar)
    sigma.hat.X2 = w2 * sigma[2]*sigmar/(sigma0+sigmar)
    
    Xe = H.mean = rep(0,N)
    sigma.vec=rep(0,N)
    for (i in (1:M)) {
      loc1 = which(data[,"group"]==i & data[,"calibration"]==0)
      loc2 = which(data[,"group"]==i & data[,"calibration"]==1)
      Xe[loc1] = w1[i]*(data[loc1,"HL"] - xi[i+1])/ gam[i+1] + (1-w1[i])*( alpha[i] + data[loc1,Wname] %*% matrix(tau,ncol=1))
      Xe[loc2] = w2[i]*(  (1-w2.in[i])*(data[loc2,"HC"]-xi[1])/gam[1] + w2.in[i]*(data[loc2,"HL"]-xi[i+1])/gam[i+1]  ) + (1-w2[i])*(alpha[i] + data[loc2,Wname] %*% matrix(tau,ncol=1))
      H.mean[loc1] = data[loc1,"HL"]
      H.mean[loc2] = (data[loc2,"HL"] + data[loc2,"HC"])/2
      sigma.vec[loc1] = sigma.hat.X1[i]
      sigma.vec[loc2] = sigma.hat.X2[i]
    }
    data0 = cbind(data,Xe,H.mean,sigma.vec)
    id=unique(data0[,"matchID"])
    mergedata=function(id,nam="H.mean") {
      a=data0[which((data0[,"Y"]==1) & (data0[,"matchID"]==id)),nam]
      b=data0[which((data0[,"Y"]==0) & (data0[,"matchID"]==id)),nam]
      c(a,b)
    }
    data.c=sapply(id, mergedata, nam="Xe", simplify = F)
    data.sigma=sapply(id, mergedata, nam="sigma.vec", simplify = F)
    if (is.null(Zname)) {
      condloglike=function(par0,data=data.c) {
        -sum(log(base::unlist(lapply(data,function(x) { 1/(1+sum(exp(par0*(x[-1]-x[1]))))}))))
      }
    } else {
      data.z1=sapply(id, mergedata, nam=Zname[1], simplify = F)
      condloglike=function(par0,data=data.c) {
        -sum(log(base::unlist(mapply(function(x,y) { 1/(1+sum(exp(par0[1]*(x[-1]-x[1]) + par0[2]*(y[-1]-y[1]) )))},data,data.z1))))
      }
    }
    if (is.null(Zname)) {
      par.c=nlminb(start=1, objective=condloglike,data=data.c)$par
      sd.c=as.vector(sqrt(1/pracma::hessian(condloglike,par.c)))
    } else {
      par.c=nlminb(start=c(1,1), objective=condloglike,data=data.c)$par
      sd.c=as.vector(sqrt(1/pracma::hessian(condloglike,par.c)))
    }
    
    ## exact calibration method
    # (1) Monte Carlo Method
    X.mat.func=function(a) {
      x=data.c[[which(id==a)]]
      y=data.sigma[[which(id==a)]]
      z=seed.mat[which(data[,"matchID"]==a),]
      x %*% t(rep(1,25)) + sqrt(y %*% t(rep(1,25)))*z
    }
    X.mat=sapply(id,X.mat.func,simplify = F)
    if (is.null(Zname)) {
      loglike.e.mc=function(par) {
        -sum( log(mapply(function(x){mean(exp(par*x[1,])/apply(exp(par*x),2,sum))},X.mat)) )
      }
    } else {
      loglike.e.mc=function(par) {
        -sum( log(mapply(function(x,y){mean(exp(par[1]*x[1,]+rep(par[2]*y[1],25))/apply(exp(par[1]*x+(par[2]*y) %*% t(rep(1,25))),2,sum))},X.mat,data.z1)) )
      }
    }
    par.e1=nlminb(start=par.c, objective=loglike.e.mc)$par
    sd.e1=as.vector(sqrt(1/pracma::hessian(loglike.e.mc,par.e1)))
    ## Gauss Hermite Method
    ratio=base::unlist(lapply(data.c,length))
    ratio2= ratio==2;ratio3= ratio==3
    mysubset=function(x,y){
      if(sum(y)==0) {
        NULL
      } else {
        t(sapply(subset(x,y),function(x) x,simplify=T))
      }
    }
    data.c2=mysubset(data.c,ratio2);data.c3=mysubset(data.c,ratio3)
    if (is.null(Zname)==0) {
      data.z12=mysubset(data.z1,ratio2);data.z13=mysubset(data.z1,ratio3)
    }
    data.sigma2=mysubset(data.sigma,ratio2);data.sigma3=mysubset(data.sigma,ratio3)
    
    weights0=weight.f(n=10,mean=data.c2[,2]-data.c2[,1],sigma=data.sigma2[,2]+data.sigma2[,1])
    
    mean.list=cbind(data.c3[,2]-data.c3[,1],data.c3[,3]-data.c3[,1])
    if (dim(mean.list)[1] !=0) {
      sigma.list=cbind(data.sigma3[,2]+data.sigma3[,1],data.sigma3[,1],data.sigma3[,1],data.sigma3[,3]+data.sigma3[,1])
      p1=matrix(0,ncol=9,nrow=dim(mean.list)[1])
      p2=matrix(0,ncol=9,nrow=dim(mean.list)[1])
      we=matrix(0,ncol=9,nrow=dim(mean.list)[1])
      for (m in (1:dim(mean.list)[1])) {
        zz=mgauss.hermite(3, mu=mean.list[m,], sigma=matrix(sigma.list[m,],ncol=2))
        p1[m,]=zz$points[,1];p2[m,]=zz$points[,2];we[m,]=zz$weights
      }
    }
    
    
    if (is.null(Zname)) {
      loglike.e=function(par) {
        if (is.null(data.c3)) {
          xx=0
        } else {
          xx=apply(log( ( 1/(1+exp(par*p1)+exp(par*p2)) ) %*% apply(we,2,mean)),2,sum)
        }
        if (is.null(data.c2)) {
          yy=0
        } else {
          yy=apply(log( ( 1/(1+exp(par*weights0$loc)) ) %*% weights0$weights),2,sum)
        }
        -yy-xx
      }
    } else {
      loglike.e=function(par) {
        if (is.null(data.c2)) {
          xx=0
        } else {
          xx=apply(log( ( 1/(1+exp(par[1]*weights0$loc + (par[2]*(data.z12[,-1]-data.z12[,1])
          )  %*% t(rep(1,10)) )) ) %*% weights0$weights),2,sum)
        }
        if (is.null(data.c3)) {
          yy=0
        } else {
          yy=apply(log( ( 1/(1+exp(par[1]*p1+ (par[2]*(data.z13[,2]-data.z13[,1])
          ) %*% t(rep(1,9)))+exp(par[1]*p2+(par[2]*(data.z13[,2]-data.z13[,1])
          ) %*% t(rep(1,9)))) ) %*% apply(we,2,mean)),2,sum)
        }
        -yy-xx
      }
    }
    
    par.e=nlminb(start=par.c, objective=loglike.e)$par
    sd.e=as.vector(sqrt(1/pracma::hessian(loglike.e,par.e)))
    
    
    res[j,]=c(par.c[1],par.e[1],par.e1[1])
    sd.mat[j,]=c(sd.c[1],sd.e[1],sd.e1[1])
  }
  out=sqrt(apply(res,2,var)+apply(sd.mat^2,2,mean))
  names(out) = c("ACM","ECM1","ECM2")
  rm("curWarnings",envir=as.environment(1))
  out
}



main_func=function(data=data,Wname="W",Zname=NULL,Y="Y",HL="HL",HC="HC",ID="matchID",group="group",calibration="calibration") {
  ## transferring the orginal datasets to the dataset that can be used in 
  ## the parameter estimation functions.

  mydata<<-data_lme(data=data,                  ## the dataset name
                    Y=Y,                      ## disease outcome
                    HL=HL,                    ## local lab measurements
                    HC=HC,                    ## central lab measurements 
                    Wname=Wname,                ## covariates in the model for X
                    ID=ID,               ## ID for the matched set or strata
                    group=group,              ## group or study number
                    calibration=calibration   ## whether this subject is in the calibration subset (Yes 1, No 0)
  )
  
  # calibration parameter estimation 
  m1=lme_mydata2(Wname=Wname,Zname=Zname)
  # log(OR) estimates
  beta.mat=get.beta(mat=mydata$data0,model=m1$model,Wname=Wname,Zname=Zname)
  # S.E. of log(OR); waiting for about 20 seconds
  sd.mat=get.beta.var(OR=1,mat=mydata$data0,model=m1$model,
                      Wname=Wname,seed.mat=beta.mat$seed.mat,Zname=Zname)
  
  ### (1) calibration parameters
  xigma.mat=nlme::random.effects(m1$model,condVar = FALSE)$lab
  M = dim(xigma.mat)[1]-1
  ## fixed effects
  fixed_effect = summary(m1$model)$coefficients
  rownames(fixed_effect) = c("(intercept)",paste("study",2:M,sep=""),Wname)
  ## random effect point estimates (EBLUP)
  xi.estimate = xigma.mat[,2]
  gamma.estimate = xigma.mat[,1]+1
  names(xi.estimate)=names(gamma.estimate)=paste("lab_",0:M,sep="")
  ## sigma_xi^2 and sigma_gamma^2
  sigma_r=as.data.frame(summary(m1$model)$varcor)[,4][c(3,2)]
  names(sigma_r) = c("sigma_xi^2","sigma_gamma^2")
  ## sigma_x^2
  sigma_x = as.data.frame(summary(m1$model)$varcor)[,4][1]
  names(sigma_x) = c("sigma_x^2")
  ## sigma_d^2
  rsigma=as.data.frame(summary(m1$model)$varcor)[,4][c(1,4)]
  w.data=m1$model@frame
  weigh=rep(1,M+1) 
  for (i in (1:M)) {
    weigh[i+1]=w.data[which(w.data[,"lab"]==i),][1,"(weights)"]
  }
  sigma_d= abs(rsigma[2] * (1/weigh) )
  names(sigma_d) = paste("sigma_",0:M,"^2",sep="")
  ## sigma_vec
  sigma_vec = c(sigma_r,sigma_x,sigma_d)
  
  ## (2) log Odds Ratio
  ### estimated log(OR) by ACM
  beta_ACM = beta.mat$par.c
  ### estimated log(OR) by ECM1
  beta_ECM1 = beta.mat$par.e
  ### estimated log(OR) by ECM2
  beta_ECM2=beta.mat$par.e1
  
  ### point estimates of log Odds ratio
  beta_point = c(beta_ACM,beta_ECM1,beta_ECM2)
  names(beta_point) = c("ACM","ECM1","ECM2")
  ### S.E. of ACM and ECMs
  beta_se=sd.mat
  
  output = list(fixed_effect=fixed_effect,gamma=gamma.estimate,xi=xi.estimate,
                sigma=sigma_vec,beta_point=beta_point,beta_se=beta_se)
  #output
  return(output)
}



