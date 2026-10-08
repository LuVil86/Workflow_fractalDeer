## ----setup, include=FALSE------------------------------------------------
require(migrateR)

## ---- fig.width = 7.5, fig.height = 3.5, echo=FALSE, warning=F-----------
pcol <- c("darkgreen","blue","purple","orange","red")
  plty <- c(2,1,3:5)
  ylabs <- c(expression("NSD "(Km^2)),"Elevation (m)")
  x1 <- matrix(1:366,ncol=1)
  parame <- pEst()
    parame$gamma <- c(0,300,400)
    parame$delta 	<- c(0,150,500)
    param <- parame[,-1]
	parame <- parame[,-c(7,8)]
	parame$delta <-	c(-400,-300,0)

  
		
  dday <- c(1:366)
  nsd <- matrix(nrow = 366, ncol = 5)
    nsd[,1] <- c(rep(1,100),seq(2,300,length.out=18),rep(300,130),seq(299,100,length.out=18),rep(100,100))+rnorm(366)
    nsd[,2] <-c(rep(1,110),seq(1,280,length.out=18),rep(280,110),seq(280,1,length.out=18),rep(1,110))
    nsd[,3] <- c(rep(1,120),seq(1,220,length.out=26),rep(220,220))+rnorm(366)	#rep(c(1,260),each=183)
    nsd[,4] <- seq(0,180,length.out=366)
    nsd[,5] <- c(1:16,rep(17,350))+rnorm(366)
  
  control1=nls.control(maxiter=150,warnOnly=T)  # sets nls controls | #tol=1e-10
  nsdm <- list()
      mm.p <- c("theta","phi","delta","rho","phi2","zeta")
      nsdm[[1]] <- nls(nsd[,1]~(
        delta/(1+exp((theta-dday)/phi)))+
        (-(delta*zeta)/(1+exp((theta+2*phi+2*phi2+rho-dday)/phi2))), 
        algorithm = "port", start = param["strt", mm.p],
        lower = param["lwr", mm.p], upper = param["upr", mm.p],
        control = control1
      )

      mig.p <- c("theta","phi","delta","rho","phi2") 
      nsdm[[2]] <- nls(nsd[,2]~(
        delta/(1+exp((theta-dday)/phi)))+
        (-delta/(1+exp((theta+2*phi+2*phi2+rho-dday)/phi2))),
        algorithm = "port", start = param["strt", mig.p],
        lower = param["lwr", mig.p], upper = param["upr", mig.p],
        control = control1
      )

      dis.p <- c("theta","phi","delta")
      nsdm[[3]] <- nls(nsd[,3]~(delta/(1+exp((theta-dday)/phi))),
        algorithm = "port", start = param["strt", dis.p],
        lower = param["lwr", dis.p], upper = param["upr", dis.p],
        control = control1
      )

      nsdm[[4]] = glm(nsd[,4]~-1+dday)		

      p.res <- c("delta","kappa")
      nsdm[[5]] <- nls(nsd[,5]~(delta*(1-exp(kappa*dday))),
        algorithm = "port",
        start = param["strt",p.res],
		lower = param["lwr",p.res], 
		upper = param["upr",p.res] ,
		control = control1
      )

      nsdf <- sapply(nsdm,fitted)
        nsdf[,1] <- nsdf[,1]+10
        nsdf[,2] <- nsdf[,2] +5

	
  # ELEV:
  elv <- matrix(nrow=366,ncol=3)
    elv[,1] <-c(rep(300,110),seq(300,10,length.out=18),rep(10,110),seq(1,300,length.out=18),rep(300,110))
    elv[,2] <- c(rep(290,100),seq(280,1,length.out=26),rep(1,240))+rnorm(366)
    elv[,3] <- rep(310,366)+rnorm(366)

  param <- parame
  
  elvm <- list()
    elvm[[1]] <- nls(elv[,1]~(gamma+delta/(1+exp((theta-dday)/phi)) - delta/(1+exp((theta+phi*2+phi2*2+rho-dday)/phi2))), 
      algorithm="port", upper=param[3,], lower=param[1,], start=param[2,], control=control1)
    elvm[[2]] <- nls(elv[,2]~(	gamma+delta/(1+exp((theta-dday)/phi))), 
			algorithm="port", upper=param[3,c(1:4)], lower=param[1,c(1:4)], start=param[2,c(1:4)],
			control=control1)
    elvm[[3]] <- glm(elv[,3]~1)


  ##----------------------------------------------------------------------------  
  par(mfrow=c(1,2), mar=c(1.3,1.3,1.3,1.3))

  plot(cbind(c(0,366),c(0,340)),typ="n",xlab="", ylab="", xaxt="n",yaxt="n", las=1)
    title(xlab="Days",ylab=ylabs[1],line=0)
    matplot(x1,nsdf,add=T,typ="l",lty=plty,lwd=2,col=pcol,xlim=c(0,365))	
  legend("topleft",legend=c("mixed migrant", "migrant","disperser","nomad","resident"), col=pcol,lty=plty,lwd=2,cex=.58)

  plot(cbind(c(0,366),c(0,340)),typ="n",xlab="", ylab="", xaxt="n",yaxt="n", las=1)
    title(xlab="Days",ylab=ylabs[2],line=0)
    matplot(x1,sapply(elvm,fitted),add=T,typ="l",lty=plty[c(2,3,5)],lwd=2,col=pcol[c(2,3,5)],xlim=c(0,365))	
  legend("bottomleft",legend=c("migrant","disperser","resident"), col=pcol[c(2,3,5)],lty=plty[c(2,3,5)],lwd=2,cex=.58)

## ------------------------------------------------------------------------
require(migrateR)
data(bighorn)
bighorn

## ---- warning=FALSE------------------------------------------------------
  bhs.nsd <- mvmtClass(bighorn, stdt = "10-31")

## ---- warning=FALSE------------------------------------------------------
  bhs.elev <- mvmtClass(bighorn, fam = "elev", stdt = "10-31")

## ----echo=FALSE----------------------------------------------------------
  bhs.nsd <- mvmtClass(bighorn, stdt = "10-31")

## ------------------------------------------------------------------------
!fullmvmt(bhs.nsd)

## ------------------------------------------------------------------------
fullmvmt(bhs.nsd, out = "numer")

## ------------------------------------------------------------------------
fullmvmt(bhs.nsd, out = "name")

## ------------------------------------------------------------------------
  pest.n2 <- pEst(s.d = 15)
  bhs.nsd2 <- mvmtClass(bighorn, stdt = "10-31", p.est = pest.n2)
  all(fullmvmt(bhs.nsd2))

## ------------------------------------------------------------------------
   fullmvmt(bhs.elev)
   pest.e2 <- pEst(s.d = -500)
   bhs.elev2 <- refine(bhs.elev, p.est = pest.e2)
   all(fullmvmt(bhs.elev2))

## ---- fig.height=4, fig.width=7, echo=FALSE------------------------------
plot(bhs.nsd[[3]])

## ---- fig.height=4, fig.width=7, echo=T----------------------------------
plot(bhs.elev[[3]])

## ---- fig.width = 7.5, fig.height = 3.5----------------------------------
   pest.e3 <- pEst(u.d = 0)
   bhs.elev3<- refine(bhs.elev, p.est = pest.e3)
   plot(bhs.elev3[[3]])

## ---- fig.width=7.5, fig.height=3.5--------------------------------------
spatmig(bighorn[2], bhs.nsd[2])


## ------------------------------------------------------------------------
bhs.nsd2

## ------------------------------------------------------------------------
top.bhs.nsd2 <- topmvmt(bhs.nsd2)
bhs.nsd2.behavior <- names(top.bhs.nsd2)
table(bhs.nsd2.behavior)

## ------------------------------------------------------------------------
summary(bhs.elev3)  # default classification shown under "topmod"

top.bhs.elev3 <- topmvmt(bhs.elev3, mrho = 21, mdelta = 500)
table(names(top.bhs.elev3))

## ---- fig.width=7.5, fig.height=3.5,echo=FALSE---------------------------
par(mfrow=c(1,2))
plot(bhs.elev3[[1]])
plot(bhs.elev3[[2]])

## ------------------------------------------------------------------------
p.bhs.nsd2 <- mvmt2df(top.bhs.nsd2)
p.bhs.nsd2

## ------------------------------------------------------------------------
t2 <- theta2(bhs.nsd2)
t2

## ----echo=FALSE----------------------------------------------------------
par(mfrow=c(1,1))

## ---- fig.width = 7.5, fig.height = 3.5----------------------------------
plot(bhs.nsd[[3]])
abline(v = t2[3,1], lty = 2, lw = 2)

## ------------------------------------------------------------------------
d2 <- delta2(bhs.nsd2)
d2

## ---- fig.width = 7.5, fig.height = 3.5----------------------------------
plot(bhs.nsd2[[2]])
d <- coef(bhs.nsd2[[2]]@models$mixmig)["delta"]
abline(h = d-d2[2,1], lty = 2, lw = 2)

## ------------------------------------------------------------------------
mvmt2dt(bhs.nsd, mod = "mixmig")

## ------------------------------------------------------------------------
bhs3.nsd <- mvmtClass(bighorn[3], stdt = "10-31", ecut = "6-15")

## ---- fig.width = 7.5, fig.height = 3.5, echo=FALSE----------------------
par(mfrow=c(1,2))
plot(bhs.nsd2[[3]])
plot(bhs3.nsd[[1]])

## ---- warning=FALSE------------------------------------------------------
rlocs <- findrloc(bighorn)
bhs.rnsd <- mvmtClass(bighorn, rloc = rlocs$rloc, stdt = "10-31", p.est = pest.n2)
bhs.rnsd2 <- refine(bhs.rnsd, pEst(s.t = 220))
fullmvmt(bhs.rnsd2,"name")


## ---- fig.width = 7.5, fig.height = 3.5, echo=FALSE----------------------
   par(mfrow=c(1,2))
   plot(bhs.nsd[[1]], ylim=c(0,45))
   plot(bhs.rnsd2[[1]], ylim=c(0,45))

## ---- echo=FALSE---------------------------------------------------------
bighorn

## ------------------------------------------------------------------------
bhs.elev0 <- mvmtClass(bighorn, fam = "elev", p.est = pest.e2)
bhs.elev0r1 <- refine(bhs.elev0, p.est = pest.e3)
fullmvmt(bhs.elev0r1, out = "name")

bhs.elev.stdt.mig <- topmvmt(bhs.elev3, omit = c("resident", "disperser"))
bhs.elev.mig <- topmvmt(bhs.elev0r1, omit = c("resident", "disperser"))

p.elev.stdt <- mvmt2df(bhs.elev.stdt.mig)
p.elev <- mvmt2df(bhs.elev.mig)

round(p.elev.stdt[[1]] - p.elev[[1]] ,2)

## ------------------------------------------------------------------------
  str(bhs.nsd[["s110 2010"]])

