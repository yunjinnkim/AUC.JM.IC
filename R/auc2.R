#' Compute AUC for Interval-Censored Data
#'
#' This function ---
#'
#' @param s description
#' @param th_v description
#' @param data3.id A dataset for analysis
#' @param survp description
#' @param survt description
#' @param ee description
#'
#' @return A list containing:
#'   \item{auc1}{}
#'
#' @import interval Icens
#' @export
AUC_ic2<-function(s,th_v,data3.id,survp,survt,ee) {

  nn=nrow(data3.id);
  nr=ncol(survp)
  nr0=nr-1
  surv0=mark2=rep(0,nn)
  TL1=data3.id$TL
  TR1=data3.id$TR

  delta=ifelse(TR1==Inf,0,1)


  surv0=rep(0,nn)

  for(i in 1:nn) {
    if(s<survt[1]) surv0[i]=1
    for(j in 2:nr) {
      if(survt[1]<s&s<=survt[j]) {surv0[i]=survp[i,j]; break}
    }
  }

  mark2=matrix(0,nn,4)
  for(k in 1:4){
    th=th_v[k]+s
    for(i in 1:nn){
      for(j in 2:nr){
        if(survt[j]<th)  mark2[i,k]=survp[i,j]
      }
    }
    if(th<survt[2])  mark2[1:nn,k]=survp[1:nn,1]
  }

  auc1<-rep(0,4)
  n1=nn-1
  for(k in 1:4){
    tot1=tot11=tot2=tot22=tot3=tot33=tot4=tot44=tot5=tot6=tot55=tot66=0
    th=th_v[k]+s
    for(i in 1:n1){
      i2=i+1
      for(j in i2:nn){
        if(s<=TL1[i]&TR1[i]<=th){
          if(th<=TL1[j])                                { tot1=tot1+1;           if(mark2[i,k]<mark2[j,k]) tot11=tot11+1 }
          if(s<=TL1[j]&TL1[j]<=th&th<TR1[j])            { tot2=tot2+mark2[j];    if(mark2[i,k]<mark2[j,k]) tot22=tot22+mark2[j,k]}
        }
        else if(s<TL1[i]&TL1[i]<th&th<TR1[i]) {
          if(th<=TL1[j])                               {tot3=tot3+(1-mark2[i,k]);            if(mark2[i,k]<mark2[j,k]) tot33=tot33+(1-mark2[i,k])}
          if(s<=TL1[j]&TL1[j]<=th&th<TR1[j])           {tot4=tot4+(1-mark2[i,k])*mark2[j,k]; if(mark2[i,k]<mark2[j,k]) tot44=tot44+(1-mark2[i,k])*mark2[j,k]}
        }
        else if(TL1[i]<s&TR1[i]<th){
          if(th<=TL1[j])                               {tot5=tot5+mark2[i,k];             if(mark2[i,k]<mark2[j,k]) tot55=tot55+mark2[i,k]}
          if(s<=TL1[j]&TL1[j]<=th&th<TR1[j])           {tot6=tot6+mark2[i,k]*mark2[j,k] ; if(mark2[i,k]<mark2[j,k]) tot66=tot66+mark2[i,k]*mark2[j,k]}
        }

      }
    }

    auc1[k]=(tot11+tot22+tot33+tot44+tot55+tot66)/(tot1+tot2+tot3+tot4+tot5+tot6)
  }

  return(list("auc1"=auc1))
}
###############################
