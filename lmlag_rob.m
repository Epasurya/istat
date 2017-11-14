function Hasil=lmlag_rob(y,x,W,alpha)
[n k] = size(x); 
if nargin==3
    alpha=0.05;
end
if nargin<3
    error('lmerror: Input Variabel Kurang');
end
[l m] = size(W);
if l~=m
    error('lmerror: Matrix W bukan matrix bujursangkar');
end
z=x'*x;                     % Menghitung Invers Matrik x'*x
xpxi=inv(z);                
b = xpxi*(x'*y);            % Hitung nilai koefisien bera OLS
M = eye(n) - x*xpxi*x';     % Hitung nilai M 
e = M*y;                    % Hitung nilai residual
sighat = (e'*e)/n;          % Hitung nilai sigma hat
T = trace((W+W')*W);        % Hitung nilai T
J = [(W*x*b)'*M*(W*x*b)+(T*sighat)];
lm1 = (e'*W*y/sighat);         % Htung nilai faktor koreksi
lm2 = (e'*W*e/sighat);
lmr1 = (lm1 - lm2);
lmr2 = lmr1*lmr1;
den  = (J/sighat) - T;
lmlag_rob = lmr2/den;          % Hitung nilai LM lag robust
prob = 1-chi2cdf(lmlag_rob,1); % Nilai probabilitas LM lag robust
chi2_tabel=chi2inv(1-alpha,1);
fprintf('Statistik Uji LM untuk spasial Lag Robust \n');
fprintf('LM Lag Robust  Chi-Square Tabel   p-value \n');
[lmlag_rob    chi2_tabel    prob]
fprintf('Kesimpulan \n');
if lmlag_rob<chi2_tabel
    fprintf('Gagal Tolak H0 \n');
else
    fprintf('Tolak H0 \n');
end
