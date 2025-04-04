#'
#'
#' This function ---
#'
#'
#' @import interval Icens stats
#' @export
####################################################################
ic_joint<-function(data2,data.id,xx,ttt,yyy,xxx) {


  TL<-data.id$TL; TR<-data.id$TR; delta=ifelse(TR==Inf,0,1)
  idd<-data2$ID
  id0<-unique(idd)
  n=length(id0)
  nn=nrow(data2)
  nobs=table(data2$ID)

  p=ncol(xx)

  y<-t<-matrix(0,n,max(nobs))

  print (c(length(yyy), length(ttt), nrow(xxx), length(data2$ID)))

  for(i in 1:n) {
    jjj=0
    for(j in 1:nn) {
      if(id0[i]==idd[j]) {jjj=jjj+1; y[i,jjj]=yyy[j]; t[i,jjj]=ttt[j]}
    }
  }

  ########################
  ####  EM algorithm  ####
  ########################
  IDD<-data2$ID

  fit.lme <-nlme::lme(yyy~ttt+xxx, random=~ttt|IDD,
                 control = list(opt = "optim", optCtrl = list(maxfun = 100000)))

  bb <- fit.lme$coef$random$ID
  beta<-fit.lme$coef$fixed

  xxx2=as.matrix(xxx)
  p=ncol(xxx2)
  var1<-stats::var(bb)
  sigmasq_e.hat<-fit.lme$sigma^2
  V<-V.hat <-sigma<-var1

  fit.interval <- interval::icfit(Surv(TL,TR, type='interval2')~1)

  nr=ncol(fit.interval$A)  #nr=30
  survt <- fit.interval$intmap[1,1:nr]
  AA<-fit.interval$A
  gamma=rep(0.001,p+1)
  lamb<-rep(1/nr,nr)
  clamb<-clam<-cumsum(lamb)

  ##################################################
  xx=as.matrix(xx)


  for(ii2 in 1:10){

    k <- 7
    gherm1 <- gherm2 <- statmod::gauss.quad(k, kind='hermite')
    sigma <- V.hat
    sigma.c <- chol(sigma)
    W1=sigma.c[1,1]; W2=sigma.c[2,2]; W3=sigma.c[1,2]

    CE1=CE2=CE11=CE12=CE22=CE4=rep(0,n)
    CE3=CE5=CE6=matrix(0,n,nr)

    pp=ee=rr=zmat=matrix(0,n,nr);

    for(i in 1:n){
      TOTT=TOTT1=TOTT2=TOTT11=TOTT22=TOTT12=TOTT4=0
      tot1=tot2=tot3=tot5=1.0

      for(k1 in 1:k){
        for(k2 in 1:k){

          VC0 <- sqrt(2)*gherm1$nodes[k1]*W1+bb[i,1]
          VC1 <- sqrt(2)*(gherm2$nodes[k2]*W2 + gherm1$nodes[k1]*W3)+bb[i,2]
          tot1 = 1
          pi = 3.14159265359
          for(j in 1:nobs[i]){
            x11=sum(c(1,t[i,j],xx[i,1:p])%*%beta)
            ss <- y[i,j]-(x11+ VC0 + VC1*t[i,j])
            tot1 <- tot1*1/sqrt(2*pi*sigmasq_e.hat)*exp(-ss^2/(2*sigmasq_e.hat))
          }

          tot3=0
          SL1=SR1=0
          nr2=nr-1

          for(j4 in 1:nr2){
            x0=c(1,survt[j4],xx[i, 1:p])
            QQQ=VC0+VC1*survt[j4]+sum(x0%*%beta)
            gamt=sum(c(xx[i,1:p],QQQ)*gamma)
            if(survt[1]>TL[i]){SL1=1}
            else{ SL1=ifelse(survt[j4]<=TL[i]&TL[i]<survt[j4+1],exp(-clam[j4]*exp(gamt)),SL1) }
          }

          for(j5 in 1:nr2){
            x0=c(1,survt[j5],xx[i,1:p])
            QQQ=VC0+VC1*survt[j5]+sum(x0%*%beta)
            gamt=sum(c(xx[i,1:p],QQQ)*gamma)
            if(survt[1]>TR[i]){SR1=0.999}
            else{ SR1=ifelse(survt[j5]<=TR[i]&TR[i]<survt[j5+1],exp(-clam[j5]*exp(gamt)),SR1) }
          }

          if(SL1==0){SL1=1}
          SL1=ifelse(SL1<=SR1,SL1+0.001,SL1)
          tot3=(SL1-SR1*delta[i])
          tot5 <-tot1*tot3
          TOTT = TOTT +     gherm1$weights[k1]*gherm2$weights[k2]*tot5
          TOTT1 = TOTT1 +   gherm1$weights[k1]*gherm2$weights[k2]*tot5*VC0
          TOTT2 = TOTT2 +   gherm1$weights[k1]*gherm2$weights[k2]*tot5*VC1
          TOTT11 = TOTT11 + gherm1$weights[k1]*gherm2$weights[k2]*tot5*VC0^2
          TOTT22 = TOTT22 + gherm1$weights[k1]*gherm2$weights[k2]*tot5*VC1^2
          TOTT12 = TOTT12 + gherm1$weights[k1]*gherm2$weights[k2]*tot5*VC0*VC1


        }
      }

      CE1[i]=TOTT1/TOTT;     CE2[i]=TOTT2/TOTT
      CE11[i]=TOTT11/TOTT;    CE12[i]=TOTT12/TOTT
      CE22[i]=TOTT22/TOTT;

      for (j in 1:nr) {
        zmat[i,j]=(CE1[i]+CE2[i]*survt[j]+sum(c(1,survt[j],xx[i,1:p])*beta))
        if(AA[i,j]==1){
          elinear=exp( sum(c(xx[i,1:p],zmat[i,j])*gamma) )
          pp[i,j]=lamb[j]*elinear*exp(-clamb[j]*elinear)
        }
      }
      if(delta[i]==1) ee[i,]=pp[i,]/sum(pp[i,])
      for (j in 1:nr){
        rr[i,j]=ifelse(delta[i]==1,sum(ee[i,j:nr]),ifelse(survt[j]<=TL[i],1,0))
      }

    }


    ####---- beta0.hat, beta1.hat, beta2.hat

    nmax=max(nobs)
    IV=matrix(0,nn,nn); V0=matrix(0,nmax,nmax)
    kkk=0
    for(i in 1:n) {
      for(j in 1:nobs[i]){
        for(j1 in j:nobs[i]){
          if(j==j1) sig=sigma[1,1]+t[i,j]^2*sigma[2,2]+2*sigma[1,2]*t[i,j]+sigmasq_e.hat
          if(j<j1)  sig=sigma[1,1]+t[i,j]*t[i,j1]*sigma[2,2]+sigma[1,2]*(t[i,j]+t[i,j1])
          V0[j,j1]=V0[j1,j]=sig;
        }}
      V1<-solve(V0[1:nobs[i], 1:nobs[i]])
      for(jj in 1:nobs[i]){
        for(jj2 in 1:nobs[i]){
          IV[kkk+jj,kkk+jj2]=V1[jj,jj2]
        }
      }
      kkk=kkk+nobs[i]
    }


    resid=rep(0,nn)
    for(i in 1:nn){
      for(j in 1:n) {
        if(id0[j]==data2$ID[i]){
          resid[i]=yyy[i]-CE2[j]*ttt[i]-CE1[j]
        } }}

    xc=cbind(rep(1,nn),ttt,xxx)
    beta.hat<-solve(t(xc)%*%IV%*%xc)%*%(t(xc)%*%IV%*%resid)

    err1=sum(abs(beta-beta.hat))


    beta[1:p+2]<-beta.hat[1:p+2]


    TOT4=TOT41=0
    for(i in 1:n){
      for(j in 1:nobs[i]){
        TOT4 =TOT4+(y[i,j]-sum(c(1,t[i,j],xx[i,])*beta)-CE1[i]-t[i,j]*CE2[i])^2
      } }

    sigmasq_e.hat <-TOT4/sum(nobs)

    ##############################################################################################
    #
    ###################################################################################################
    nmax=max(nobs)
    IV=matrix(0,nn,nn); V0=matrix(0,nmax,nmax)
    kkk=0
    for(i in 1:n) {
      for(j in 1:nobs[i]){
        for(j1 in j:nobs[i]){
          if(j==j1) sig=sigma[1,1]+t[i,j]^2*sigma[2,2]+2*sigma[1,2]*t[i,j]+sigmasq_e.hat
          if(j<j1)  sig=sigma[1,1]+t[i,j]*t[i,j1]*sigma[2,2]+sigma[1,2]*(t[i,j]+t[i,j1])
          V0[j,j1]=V0[j1,j]=sig;
        }}
      V1<-solve(V0[1:nobs[i], 1:nobs[i]])
      for(jj in 1:nobs[i]){
        for(jj2 in 1:nobs[i]){
          IV[kkk+jj,kkk+jj2]=V1[jj,jj2]
        }
      }
      kkk=kkk+nobs[i]
    }

    sigma_be<-matrix(0,p+2,p+2)
    #see=sqrt(diag(solve(t(xc)%*%(IV)%*%xc)))
    sigma_be<-solve(t(xc)%*%(IV)%*%xc)
    # sigma_be<--solve(t(xc)%*%xc)*sigmasq_e.hat


    BB=see=rep(0,p+2)
    TOT1=TOT2=TOT3=TOT4=TOT44=TOT5=0
    ntt=sum(nobs)

    for(i in 1:n){
      for(j in 1:nobs[i]){
        bet0=sum(c(1,t[i,j],xx[i,])*beta)
        TOT1=TOT1+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])*1
        TOT2=TOT2+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])*t[i,j]
        TOT3=TOT3+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])*xx[i,1]
        if(p>1){ TOT4=TOT4+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])*xx[i,2]
        TOT5=TOT5+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])*xx[i,3]   }

        TOT44=TOT44+(y[i,j]-bet0-CE1[i]-t[i,j]*CE2[i])^2
      }
    }
    BB[1]=TOT1; BB[2]=TOT2; BB[3]=TOT3;
    if(p>1) {BB[4]=TOT4; BB[5]=TOT5}

    D=ntt/2*sigmasq_e.hat^2-(TOT44)/sigmasq_e.hat^3

    S_INV1=sigma_be+sigma_be%*%BB%*%solve(D-t(BB)%*%sigma_be%*%BB)%*%t(BB)%*%sigma_be
    see[1]=sqrt(S_INV1[1,1]);    see[2]=sqrt(S_INV1[2,2]);   see[3]=sqrt(S_INV1[3,3]);
    if(p>1) {see[4]=sqrt(S_INV1[4,4]);  see[5]=sqrt(S_INV1[5,5])}

    # S_INV2=sqrt(solve(D-t(BB)%*%sigma_be%*%BB))
    beta0se=see[1];  beta1se=see[2]; beta2se=see[3];
    if(p>1) {beta3se=see[4]; beta4se=see[5]}

    b00<-CE1 ; b11<-CE2
    zmat <-matrix(0,n,nr)

    for(i in 1:n) {
      for(j in 1:nr){
        bet0=sum(c(1,survt[j],xx[i,1:p])*beta)
        zmat[i,j]=(b00[i]+b11[i]*survt[j]+bet0)
      }
    }

    p2=p+1

    barx=rep(0,nr); barx1=barx11=matrix(0,nr,p2)


    for(kk in 1:p2){
      for(j in 1:nr){
        tot=tot1=0
        for(i in 1:n){
          x=c(xx[i,],zmat[i,j])
          elinear1=exp(sum(x*gamma))
          tot=tot+elinear1*rr[i,j]
          tot1=tot1+elinear1*rr[i,j]*x[kk]
        }
        barx[j]=tot; barx1[j,kk]=tot1;
      }}

    barx11=array(0,dim=c(p2,p2,nr))

    for(kk in 1:p2){
      for(kk1 in kk:p2){
        for(j in 1:nr){
          tot11=0
          for(i in 1:n){
            x=c(xx[i,],zmat[i,j])
            elinear1=exp(sum(x*gamma))
            tot11=tot11+elinear1*rr[i,j]*x[kk]*x[kk1]
          }
          barx11[kk,kk1,j]=barx11[kk1,kk,j]=tot11
        }}}

    g=rep(0,p2)
    for(kk in 1:p2){
      tot1=0
      for(j in 1:nr){
        for(i in 1:n){
          if(delta[i]==1){
            x=c(xx[i,],zmat[i,j])
            tot1=tot1+ee[i,j]*(x[kk]-barx1[j,kk]/barx[j])
          }
        }
      }
      g[kk]=tot1
    }


    Jacob=matrix(0,p2,p2)
    for(kk in 1:p2){
      for(kk1 in kk:p2){
        tot11=0
        for(i in 1:n){
          if(delta[i]==1){
            for(j in 1:nr){
              tot11<-tot11+ee[i,j]*(barx11[kk,kk1,j]/barx[j]-barx1[j,kk]*barx1[j,kk1]/barx[j]^2)
            } }
        }
        Jacob[kk,kk1]=Jacob[kk1,kk]=tot11
      }}


    gamma.hat =gamma+solve(Jacob)%*%g
    gamma<-gamma.hat

    lamb<-rep(0,nr)

    for(j in 1:nr) {
      tot=tot1=0
      for(i in 1:n) {
        x=c(xx[i,],zmat[i,j])
        elinear=exp(sum(x*gamma))
        if(delta[i]==1) tot=tot+ee[i,j]
        tot1=tot1+rr[i,j]*elinear

      }
      lamb[j]=tot/tot1
    }

    err2=sum(abs(solve(Jacob)%*%g))
    clamb=clam=cumsum(lamb)

    print(c(ii2,err1/sum(abs(beta)),err2, gamma))
    gammase=sqrt(diag(solve(Jacob)))

    if(min(err1/sum(abs(beta)),err2/sum(abs(gamma)))<0.001) break   ## simulation
  }


  return(list("gamma"=gamma,"clamb"=clamb, "time"=survt, "ee"=ee, "beta"=beta, "zmat"=zmat, "se2"=gammase, "se1"=see, "Jac"=solve(Jacob),"Cov"=S_INV1,
              "CE1"=CE1,"CE2"=CE2))

}

