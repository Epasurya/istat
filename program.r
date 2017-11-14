#converted by suryantoepa@yahoo.co.id
#last update : 25 Sept 2017
options(digits=6)

gmm_sdm_sim=function (y1,y2,x1,x2,W,p)
{
  y1=as.matrix(y1)
  y2=as.matrix(y2)
  x1=as.matrix(x1)
  x2=as.matrix(x2)
  W=as.matrix(W)
  y3=cbind(y2,x1);
  n=nrow(y3)
  nvar=ncol(y3)
  
  #tahap pertama
  wy1=W%*%y1;
  wy2=W%*%y2;
  wx1=W%*%x1;
  wx2=W%*%x2;
  wwx2=W%*%W%*%x1;
  
  s=diag(n)-W*p;
  sinv=solve(s);
  g=W%*%W%*%sinv;
  gy2=g%*%y2;
  gx1=g%*%x1;
  gx2=g%*%x2;
  
  z1=cbind(y2,x1, wy2, wx1, wy1);
  z3=cbind(y2, x1, wy2, wx1);
  H=cbind(x1, x2 ,wx1, wx2)
  Hinv=solve(t(H)%*%H);
  PH=H%*%Hinv%*%t(H);
  wy1hat=PH%*%wy1;
  zhat=cbind(y2, x1 ,wy2, wx1, wy1hat)
  zhat1=solve(t(zhat)%*%zhat);
  deltahat=zhat1%*%t(zhat)%*%y1;	#parameter spasial lag witohut b0
  deltahat=t(deltahat);
  yhat=deltahat%*%t(z1);
  yhat=t(yhat);
  uhat=y1-yhat;	                    #residual
  uhat1=t(uhat)%*%uhat;
  var1=uhat1/n;
  
  q=cbind(y2, x1, x2, wy2, wx1, wx2, gy2, gx1, gx2);
  q1=t(q)%*%q;
  q1inv=solve(q1);
  A=q1inv*as.numeric(var1);
  
  q2=cbind(y2, x1, x2, wy2, wx1, wx2, gy2, gx1, gx2);
  x4=t(y3)%*%q2%*%A%*%t(q2)%*%y3;
  x4inv=solve(x4);
  x5=y3%*%x4inv%*%t(y3)%*%q2%*%A%*%t(q2);
  M=diag(n)-x5;
  M1=t(wy1)%*%M%*%q2%*%A%*%t(q2)%*%M%*%wy1;
  M1inv=solve(M1);
  rho=M1inv%*%t(wy1)%*%t(M)%*%q2%*%A%*%t(q2)%*%M%*%y1;
  
  rho1=diag(n)-W*as.numeric(rho);
  z5=t(z3)%*%q2%*%A%*%t(q2)%*%z3;
  z5inv=solve(z5);
  delta1hat=z5inv%*%t(z3)%*%q2%*%A%*%t(q2)%*%rho1%*%y1;
  delta2hat=cbind(t(delta1hat), rho);
  delta2hat=t(delta2hat);
  y2hat=t(delta2hat)%*%t(z1);
  y2hat=t(y2hat);
  mean_obs1=cbind(mean(y2), t(colMeans(x1)), mean(wy2), t(colMeans(wx1)), mean(wy1));
  
  #menghitung b0
  b01=mean(y1)-(t(delta2hat)%*%t(mean_obs1));
  satu=matrix(rep(1,n),ncol=1)
  obs1=cbind(satu, y2, x1, wy2, wx1, wy1);
  par1=rbind(b01,delta2hat);
  
  #parameter tahap pertama (x1,..xnvar,lambda,b0)
  
  y3hat=t(par1)%*%t(obs1);
  y3hat=t(y3hat);
  uhat2=y1-y3hat;	# residual
  uhat3=t(uhat2)%*%uhat2;
  uhat4=as.numeric(uhat3)*q1inv;
  
  #menentukan varian bi
  XZ=t(obs1)%*%q;
  ZX=t(q)%*%obs1;
  V0=XZ%*%uhat4%*%ZX;
  V=solve(V0);
  cii=diag(V);		
  se_bi=sqrt(cii);
  t0=par1/se_bi;
  t1=abs(t0);
  k=(nvar+1)*2;
  pval=1-pt(t1,n-k);
  pval1=2*pval;
  ttab=qt(0.9,n-1);
  
  #menghitung R-square
  sst0=y1-mean(y1);
  sst=t(sst0)%*%sst0;
  sse=uhat3;
  rsqr=1-(sse/sst);
  t1=t1;
  par1=par1;
  pval1=pval1;
  cat('*******************************************','\n')
  cat('Estimasi GMM model simultan spasial durbin','\n')
  cat('*******************************************','\n')
  Hasil_akhir=cbind(par1,se_bi,t0,pval1)
  print(Hasil_akhir)
  list(y1=y1,nobs=n,nvar=nvar,par1=par1,par1=par1,t1=t1,ttab=ttab,pval1=pval1,rsqr=rsqr,sse=sse,meth='SDM',resid=uhat2,yhat=y3hat)
}

s2sls_sdm=function (y1,y2,x1,x2,W)
{
  y1=as.matrix(y1)
  y2=as.matrix(y2)
  x1=as.matrix(x1)
  x2=as.matrix(x2)
  W=as.matrix(W)
  y3=cbind(y2,x1);
  n=nrow(y3)
  nvar=ncol(y3)
  
  #tahap pertama
  wy1=W%*%y1;
  wy2=W%*%y2;
  wx1=W%*%x1;
  wx2=W%*%x2;
  z1=cbind(y2 ,x1 ,wy2 ,wx1 ,wy1)
  H=cbind(y2 ,x1 ,x2 ,wy2 ,wx1 ,wx2)
  Hinv=solve(t(H)%*%H);
  PH=H%*%Hinv%*%t(H);
  wy1hat=PH%*%wy1;
  zhat=cbind(y2 ,x1 ,wy2 ,wx1 ,wy1hat)
  zhat1=solve(t(zhat)%*%zhat);
  deltahat=zhat1%*%t(zhat)%*%y1;	#parameter spasial lag witohut b0
  deltahat=t(deltahat);
  yhat=deltahat%*%t(z1);
  yhat=t(yhat);
  mean_obs1=cbind(mean(y2), t(colMeans(x1)), mean(wy2) ,t(colMeans(wx1)), mean(wy1));
  
  #menghitung b0
  b01=mean(y1)-(deltahat%*%t(mean_obs1));
  satu=matrix(rep(1,n),ncol=1)
  obs1=cbind(satu, y2 ,x1 ,wy2 ,wx1 ,wy1)
  par1=rbind(b01,t(deltahat))
  
  #parameter tahap pertama (x1,..xnvar,lambda,b0)
  y1hat=t(par1)%*%t(obs1);
  y1hat=t(yhat);
  uhat=y1-t(y1hat);	#residual
  uhat1=t(uhat)%*%uhat;
  var1=uhat1/n;
  
  #menentukan varian bi
  obs2=t(obs1)%*%obs1;      #X'*X
  obs2=solve(obs2);
  cii=diag(obs2);		#elemen diagonal X'*X
  var_bi=var1%*%cii;
  se_bi=sqrt(var_bi);
  t0=par1/t(se_bi);
  t1=abs(t0);
  k=(nvar+1)*2;
  pval=1-pt(t1,n-k);
  pval1=2*pval;
  ttab=qt(0.9,n-1);
  
  #menghitung R-square
  sst0=y1-mean(y1);
  sst=t(sst0)%*%sst0;
  sse=uhat1;
  rsqr=1-(sse/sst);
  t1=t1;
  par1=par1;
  pval1=pval1;
  cat('********************************************','\n')
  cat('Estimasi S2SLS model simultan spasial durbin','\n')
  cat('********************************************','\n')
  Hasil_akhir=cbind(par1 ,t(se_bi) ,t0 ,pval1)  
  print(Hasil_akhir)
  list(y1=y1,nobs=n,nvar=nvar,par1=par1,par1=par1,t1=t1,ttab=ttab,pval1=pval1,rsqr=rsqr,sse=sse,meth='S2SLS SDM',resid=uhat,yp=y1hat)
}

lm_error=function (y,x,W,alpha=0.05)
{
  if (nargs()==3) alpha=0.05
  if (nargs()<3) {stop('lm_error : Input Variabel Kurang')}
  y=as.matrix(y)
  x=as.matrix(x)
  W=as.matrix(W)
  n=nrow(x)
  k=ncol(x)
  l=ncol(W)
  m=nrow(W)
  if (l!=m) { stop('lm_error : matrix W bukan matrix bujur sangkar')
    stop	}
  
  z=t(x)%*%x;                     # Menghitung Invers Matrik x'*x
  xpxi=solve(z);                
  b = xpxi%*%(t(x)%*%y);            # Hitung nilai koefisien beta OLS
  M = diag(n) - x%*%xpxi%*%t(x);    # Hitung Nilai M 
  e = M%*%y;                   	  # Hitung nilai residual
  sighat = (t(e)%*%e)/n;          # Hitung nilai sigma hat
  
  T = sum(diag((W+t(W))%*%W));        # Hitung nilai penyebut
  lm1 = (t(e)%*%W%*%e)/sighat;# Hitung nilai pembilang
  print(lm1)
  lmerr = (lm1%*%lm1)%*%(1/T);      # Hasil LM error
  prob = 1-pchisq(lmerr,1);    # Nilai probabilitas LM error
  chi2_tabel=qchisq((1-alpha),1);
  result=(cbind(lmerr,chi2_tabel,prob))
  colnames(result)<-c('LM_Error','Chi-Square_Tabel','p-value')
  cat("Statistik Uji LM untuk spasial error","\n")  
  print(result)
  cat('Kesimpulan ','\n');
  if (lmerr<chi2_tabel){ cat('Gagal Tolak H0 ','\n')}
  else { cat('Tolak H0 ','\n')}
}

lm_lag=function (y,x,W,alpha=0.05)
{
  if (nargs()==3) alpha=0.05
  if (nargs()<3) {stop('lm_lag : Input Variabel Kurang')}
  y=as.matrix(y)
  x=as.matrix(x)
  W=as.matrix(W)
  n=nrow(x)
  k=ncol(x)
  l=ncol(W)
  m=nrow(W)
  if (l!=m) { stop('lm_lag : matrix W bukan matrix bujur sangkar')
    stop	}
  
  z=t(x)%*%x;                       # Menghitung Invers Matrik x'*x
  xpxi=solve(z);                
  b = xpxi%*%(t(x)%*%y);            # Hitung nilai koefisien beta OLS
  M = diag(n) - x%*%xpxi%*%t(x);    # Hitung Nilai M 
  e = M%*%y;                   	    # Hitung nilai residual
  sighat = (t(e)%*%e)/n;            # Hitung nilai sigma hat
  
  T = sum(diag((W+t(W))%*%W));      # Hitung nilai penyebut
  J = t(W%*%x%*%b)%*%M%*%(W%*%x%*%b)+(T%*%sighat);      
  lm1 = (t(e)%*%W%*%y)/sighat;      # Hitung nilai pembilang
  lmlag = (lm1*lm1)*(1/(J/sighat));      # Hasil LM lag
  prob = 1-pchisq(lmlag,1);           # Nilai probabilitas LM lag
  chi2_tabel=qchisq((1-alpha),1);
  result=(cbind(lmlag,chi2_tabel,prob))
  colnames(result)<-c('LM_Lag','Chi-Square_Tabel','p-value')
  cat("Statistik Uji LM untuk spasial Lag","\n")  
  print(result)
  cat('Kesimpulan ','\n');
  if (lmlag<chi2_tabel){ cat('Gagal Tolak H0 ','\n')}
  else { cat('Tolak H0 ','\n')}
}

lmerr_rob=function (y,x,W,alpha=0.05)
{
  if (nargs()==3) alpha=0.05
  if (nargs()<3) {stop('lm_error : Input Variabel Kurang')}
  y=as.matrix(y)
  x=as.matrix(x)
  W=as.matrix(W)
  n=nrow(x)
  k=ncol(x)
  l=ncol(W)
  m=nrow(W)	
  if (l!=m) { stop('lm_lag : matrix W bukan matrix bujur sangkar')
    stop	}
  
  z=t(x)%*%x;                       # Menghitung Invers Matrik x'*x
  xpxi=solve(z);                
  b = xpxi%*%(t(x)%*%y);            # Hitung nilai koefisien beta OLS
  M = diag(n) - x%*%xpxi%*%t(x);    # Hitung Nilai M 
  e = M%*%y;                   	    # Hitung nilai residual
  sighat = (t(e)%*%e)/n;            # Hitung nilai sigma hat  
  T = sum(diag((W+t(W))%*%W));      # Hitung nilai penyebut
  J = t(W%*%x%*%b)%*%M%*%(W%*%x%*%b)+(T%*%sighat)
  lm1 = (t(e)%*%W%*%e/sighat);      # Hitung nilai faktor koreksi
  lm2 = T%*%sighat%*%solve(J);
  lm3 = (t(e)%*%W%*%y/sighat);
  lmr1 = (lm1 - (lm2%*%lm3));
  lmr2 = lmr1%*%lmr1;
  den = T%*%(1-T%*%sighat%*%solve(J));
  lmerr_rob = lmr2/den;             # Hasil LM error robust
  prob = 1-pchisq(lmerr_rob,1);       # Nilai probabilitas LM error robust
  chi2_tabel=qchisq((1-alpha),1);
  result=(cbind(lmerr_rob,chi2_tabel,prob))
  colnames(result)<-c('LM_Error_Robust','Chi-Square_Tabel','p-value')
  cat("Statistik Uji LM untuk spasial error Robust","\n")  
  print(result)
  cat('Kesimpulan ','\n');
  if (lmerr_rob<chi2_tabel){ cat('Gagal Tolak H0 ','\n')}
  else { cat('Tolak H0 ','\n')}
}

lmlag_rob=function (y,x,W,alpha=0.05)
{
  if (nargs()==3) alpha=0.05
  if (nargs()<3) {stop('lm_error : Input Variabel Kurang')}
  y=as.matrix(y)
  x=as.matrix(x)
  W=as.matrix(W)
  n=nrow(x)
  k=ncol(x)
  l=ncol(W)
  m=nrow(W)	
  if (l!=m) { stop('lm_lag : matrix W bukan matrix bujur sangkar')
    stop	}
  
  z=t(x)%*%x;                       # Menghitung Invers Matrik x'*x
  xpxi=solve(z);                
  b = xpxi%*%(t(x)%*%y);            # Hitung nilai koefisien beta OLS
  M = diag(n) - x%*%xpxi%*%t(x);    # Hitung Nilai M 
  e = M%*%y;                   	    # Hitung nilai residual
  sighat = (t(e)%*%e)/n;            # Hitung nilai sigma hat  
  T = sum(diag((W+t(W))%*%W));      # Hitung nilai penyebut
  J = t(W%*%x%*%b)%*%M%*%(W%*%x%*%b)+(T%*%sighat)
  lm1 = (t(e)%*%W%*%y/sighat);      # Hitung nilai faktor koreksi
  lm2 = (t(e)%*%W%*%e/sighat);
  lmr1 = (lm1 - lm2);
  lmr2 = lmr1*lmr1;
  den  = (J/sighat) - T;
  lmlag_rob = lmr2/den;             # Hitung nilai LM lag robust
  prob = 1-pchisq(lmlag_rob,1);       # Nilai probabilitas LM error robust
  chi2_tabel=qchisq((1-alpha),1);
  result=(cbind(lmlag_rob,chi2_tabel,prob))
  colnames(result)<-c('LM_Lag_Robust','Chi-Square_Tabel','p-value')
  cat("Statistik Uji LM untuk spasial error Lag","\n")  
  print(result)
  cat('Kesimpulan ','\n');
  if (lmlag_rob<chi2_tabel){ cat('Gagal Tolak H0 ','\n')}
  else { cat('Tolak H0 ','\n')}
}

