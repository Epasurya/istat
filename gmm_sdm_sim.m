function results=gmm_sdm_sim(y1,y2,x1,x2,W,p)
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
wwx2=W*W*x1;

s=eye(n)-W*p;
sinv=inv(s);
g=W*W*sinv;
gy2=g*y2;
gx1=g*x1;
gx2=g*x2;

z1=[y2 x1 wy2 wx1 wy1];
z3=[y2 x1 wy2 wx1];
H=[x1 x2 wx1 wx2];
Hinv=inv(H'*H);
PH=H*Hinv*H';
wy1hat=PH*wy1;
zhat=[y2 x1 wy2 wx1 wy1hat];
zhat1=inv(zhat'*zhat);
deltahat=zhat1*zhat'*y1;	%parameter spasial lag witohut b0
deltahat=deltahat';
yhat=deltahat*z1';
yhat=yhat';
uhat=y1-yhat;	%residual
uhat1=uhat'*uhat;
var1=uhat1/n;

q=[y2 x1 x2 wy2 wx1 wx2 gy2 gx1 gx2];
q1=q'*q;
q1inv=inv(q1);
A=q1inv*var1;

q2=[y2 x1 x2 wy2 wx1 wx2 gy2 gx1 gx2];
x4=y3'*q2*A*q2'*y3;
x4inv=inv(x4);
x5=y3*x4inv*y3'*q2*A*q2';
M=eye(n)-x5;
M1=wy1'*M*q2*A*q2'*M*wy1;
M1inv=inv(M1);
rho=M1inv*wy1'*M'*q2*A*q2'*M*y1;

rho1=eye(n)-W*rho;
z5=z3'*q2*A*q2'*z3;
z5inv=inv(z5);
delta1hat=z5inv*z3'*q2*A*q2'*rho1*y1;
delta2hat=[delta1hat' rho];
delta2hat=delta2hat';
y2hat=delta2hat'*z1';
y2hat=y2hat';
mean_obs1=[mean(y2) mean(x1) mean(wy2) mean(wx1) mean(wy1)];

%menghitung b0
b01=mean(y1)-(delta2hat'*mean_obs1');
obs1=[ones(114,1) y2 x1 wy2 wx1 wy1];
par1=[b01;delta2hat];

%parameter tahap pertama (x1,..xnvar,lambda,b0)
results.par1=par1;
y3hat=par1'*obs1';
y3hat=y3hat';
uhat2=y1-y3hat;	%residual
uhat3=uhat2'*uhat2;
uhat4=uhat3*q1inv;

%menentukan varian bi
XZ=obs1'*q;
ZX=q'*obs1;
V0=XZ*uhat4*ZX;
V=inv(V0);
cii=diag(V);		
se_bi=sqrt(cii);
t0=par1./se_bi;
t1=abs(t0);
k=(nvar+1)*2;
pval=1-tcdf(t1,n-k);
pval1=2*pval;
ttab=tinv(0.9,n-1);

%menghitung R-square
sst0=y1-mean(y1);
sst=sst0'*sst0;
sse=uhat3;
rsqr=1-(sse/sst);
t1=t1;
par1=par1;
pval1=pval1;
fprintf('******************\n')
fprintf('Estimasi GMM model simultan spasial durbin\n')
fprintf('******************\n')
results.par1=par1;
results.t1=t1;
results.ttab=ttab;
results.pval1=pval1;
results.rsqr=rsqr;
results.sse=sse;
results.meth='SDM'
results.resid=uhat2;
results.yhat=y3hat;
Hasil_akhir=[par1 se_bi t0 pval1]