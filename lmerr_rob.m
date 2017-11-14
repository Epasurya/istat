function Hasil=lmerr_rob(y,x,W,alpha)
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
b = xpxi*(x'*y);            % Hitung nilai koefisiean beta OLS
M = eye(n) - x*xpxi*x';     % Hitung nilai M 
e = M*y;                    % Hitung nilai residual 
sighat = (e'*e)/n;          % Hitung nilai sigma hat
T = trace((W+W')*W);        % Hitung nilai T
J = [(W*x*b)'*M*(W*x*b)+(T*sighat)];
lm1 = (e'*W*e/sighat);      % Hitung nilai faktor koreksi
lm2 = T*sighat*inv(J);
lm3 = (e'*W*y/sighat);
lmr1 = (lm1 - (lm2*lm3));
lmr2 = lmr1*lmr1;
den = T*(1-T*sighat*inv(J));
lmerr_rob = lmr2/den;            % Hasil LM error robust
prob = 1-chi2cdf(lmerr_rob,1); % Nilai probabilitas LM error robusut
chi2_tabel=chi2inv(1-alpha,1);
fprintf('Statistik Uji LM untuk spasial Error Robust \n');
fprintf('LM Error Robust  Chi-Square Tabel   p-value \n');
[lmerr_rob    chi2_tabel    prob]
fprintf('Kesimpulan \n');
if lmerr_rob<chi2_tabel
    fprintf('Gagal Tolak H0 \n');
else
    fprintf('Tolak H0 \n');
end
