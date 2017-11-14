function results=s2sls_sdm(y1,y2,x1,x2,W)
y3=[y2 x1];
[n nvar]=size(y3);
results.y1=y1;
results.nobs=n;
results.nvar=nvar;

%tahap pertama
wy1=W*y1;
wy2=W*y2;
wx1=W*x1;
wx2=W*x2;
z1=[y2 x1 wy2 wx1 wy1];
H=[y2 x1 x2 wy2 wx1 wx2];
Hinv=inv(H'*H);
PH=H*Hinv*H';
wy1hat=PH*wy1;
zhat=[y2 x1 wy2 wx1 wy1hat];
zhat1=inv(zhat'*zhat);
deltahat=zhat1*zhat'*y1;	%parameter spasial lag witohut b0
deltahat=deltahat';
yhat=deltahat*z1';
yhat=yhat';
mean_obs1=[mean(y2) mean(x1) mean(wy2) mean(wx1) mean(wy1)];

%menghitung b0
b01=mean(y1)-(deltahat*mean_obs1');
obs1=[ones(114,1) y2 x1 wy2 wx1 wy1];
par1=[b01;deltahat'];

%parameter tahap pertama (x1,..xnvar,lambda,b0)
results.par1=par1;
y1hat=par1'*obs1';
y1hat=yhat';
uhat=y1-y1hat';	%residual
uhat1=uhat'*uhat;
var1=uhat1/n;

%menentukan varian bi
obs2=obs1'*obs1;%X'*X
obs2=inv(obs2);
cii=diag(obs2);		%elemen diagonal X'*X
var_bi=var1*cii;
se_bi=sqrt(var_bi);
t0=par1./se_bi;
t1=abs(t0);
k=(nvar+1)*2;
pval=1-tcdf(t1,n-k);
pval1=2*pval;
ttab=tinv(0.9,n-1);

%menghitung R-square
sst0=y1-mean(y1);
sst=sst0'*sst0;
sse=uhat1;
rsqr=1-(sse/sst);
t1=t1;
par1=par1;
pval1=pval1;
fprintf('******************\n')
fprintf('Estimasi S2SLS model simultan spasial durbin\n')
fprintf('******************\n')
results.par1=par1;
results.t1=t1;
results.ttab=ttab;
results.pval1=pval1;
results.rsqr=rsqr;
results.sse=sse;
results.meth='S2SLS SDM'
results.resid=uhat;
results.yp=y1hat;
Hasil_akhir=[par1 se_bi t0 pval1]