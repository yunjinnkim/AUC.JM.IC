#' Compute Brier Score for Survival Analysis
#'
#' This function calculates
#'
#' @param data.id A dataset for analysis
#' @param survt description
#' @param survp description
#' @param ee description
#' @param s description
#' @param th_v description
#'
#' @return A list containing:
#'   \item{bs}{}
#'
#'
#' @import survival
#' @export
brier<-function(data.id,survt,survp,ee,s,th_v) {

  nn=nrow(data.id)
  nr=ncol(ee)
  mt=length(th_v)

  TL1=data.id$TL; TR1=data.id$TR

  delta=ifelse(TR1==Inf,0,1)

  surv0=rep(0,nn)
  for(i in 1:nn) {
    if(s<survt[1]) surv0[i]=1
    for(j in 2:nr) {
      if(survt[1]<=s&s<=survt[j]) {surv0[i]=survp[i,j]; break}
    }
  }


  mark2=ff=matrix(0,nn,4)
  for(k in 1:4){
    th=th_v[k]+s
    for(i in 1:nn){
      for(j in 2:nr){
        if(survt[j]<th)   {mark2[i,k]=survp[i,j]; ff[i,k]=1-sum(ee[i,1:j])}
      }}
    if(th<survt[2])       mark2[1:nn,k]=survp[1:nn,1]
  }

  ## estimate of censoring
  CC=ifelse(delta==1,TR1,TL1); delc=1-delta
  fitc<-survival::survfit(Surv(CC,delc)~1)

  nc=length(fitc$time)-1
  indd<-max(ifelse(fitc$time>=s,1,0)*fitc$surv)

  condG<-matrix(0,nn,4)


  for(k in 1:4){
    indd2<-max(ifelse(fitc$time>=th_v[k],1,0)*fitc$surv)

    for(i in 1:nn){
      for(j in 1:nc){
        j2=j+1
        if(delta[i]==1&fitc$time[j]<=CC[i]&CC[i]<=fitc$time[j2]) condG[i,k]=fitc$surv[j]/indd
        else if(delta[i]==0)  condG[i,k]=indd2/indd
      }
    }
  }


  bs<-rep(0,4)

  for(k in 1:4){
    th=th_v[k]+s
    tot=0; tot1=tot2=tot3=tot4=tot5=0; ind1=ind2=ind3=ind4=ind5=0
    for(i in 1:nn){
      if(s<=TL1[i])                            {
        tot=tot+1
        if(th<=TL1[i])                                 {tot1=tot1+(1-mark2[i,k])^2; ind1=ind1+1}
        else if(TL1[i]<th&th<=TR1[i]&delta[i]==1)      {tot2=tot2+(ff[i,k]-mark2[i,k])^2; ind2=ind2+1}
        else if(TR1[i]<th)                             {tot3=tot3+(0-mark2[i,k])^2 ;ind3=ind3+1}
        else if(TL1[i]<th&delta[i]==0)                 {tot4=tot4+mark2[i,k]*(1-mark2[i,k])^2+(1-mark2[i,k])*(0-mark2[i,k])^2 ; ind4=ind4+1}
      }
      if(TL1[i]<s&th<TR1[i])                   {tot=tot+1; tot5=tot5+surv0[i]*mark2[i,k]*(1-mark2[i,k])^2; ind5=ind5+1}
    }
    bs[k]=(tot1+tot2+tot3+tot4+tot5)/tot
    # print(c(k,tot,ind1,ind2,ind3,ind4,ind5))
  }

  return(list("bs"=bs))

}
