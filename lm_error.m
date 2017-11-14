function Hasil=lm_error(y,x,W,alpha)
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
z=x'*x;                     % Menghitung Invers Matrik x'*x
xpxi=inv(z);                
b = xpxi*(x'*y);           % Hitung nilai koefisien beta OLS
M = eye(n) - x*xpxi*x';    % Hitung Nilai M 
e = M*y;                   % Hitung nilai residual
sighat = (e'*e)/n;         % Hitung nilai sigma hat
 
T = trace((W+W')*W);        % Hitung nilai penyebut
lm1 = (e'*W*e)/sighat;      % Hitung nilai pembilang
lmerr = (lm1*lm1)*(1/T);    % Hasil LM error
prob = 1-chi2cdf(lmerr,1);  % Nilai probabilitas LM error
chi2_tabel=chi2inv(1-alpha,1);
fprintf('Statistik Uji LM untuk spasial error \n');
fprintf('LM Error     Chi-Square Tabel       p-value \n');
[lmerr    chi2_tabel     prob]
fprintf('Kesimpulan \n');
if lmerr<chi2_tabel
    fprintf('Gagal Tolak H0 \n');
else
    fprintf('Tolak H0 \n');
end
