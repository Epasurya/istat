function Hasil=lm_lag(y,x,W,alpha)
[n k] = size(x); 
if nargin==3
    alpha=alpha;
end
if nargin<3
    error('lmerror: Input Variabel Kurang');
end
[l m] = size(W);
if l~=m
    error('lmerror: Matrix W bukan matrix bujursangkar');
end
z=x'*x;                   % Menghitung Invers Matrik x'*x
xpxi=inv(z);         
b = xpxi*(x'*y);          % Hitung nilai koefisien Beta OLS
M = eye(n) - x*xpxi*x';   % Hitung nilai M 
e = M*y;                  % Hitung nilai residual
sighat = (e'*e)/n;        % Hitung nilai sigma hat
 
T = trace((W+W')*W);           % Hitung nilai T
J = [(W*x*b)'*M*(W*x*b)+(T*sighat)];      
lm1 = (e'*W*y)/sighat;       % Hitung nilai pembilang
lmlag = (lm1*lm1)*(1/(J/sighat));   % Hasil LM lag
prob = 1-chi2cdf(lmlag,1);   % Nilai probabilitas LM error
chi2_tabel=chi2inv(1-alpha,1);
fprintf('Statistik Uji LM untuk spasial lag \n');
fprintf('LM Lag     Chi-Square Tabel       p-value \n');
[lmlag    chi2_tabel     prob]
fprintf('Kesimpulan \n');
if lmlag<chi2_tabel
    fprintf('Gagal Tolak H0 \n');
else
    fprintf('Tolak H0 \n');
end
