global M
global N
N = 250;
M = N+1;

smrtnost = 1.0;
dt = 0.001;
K = 150;

function P = verjetnost(sprememba, povprecje)
  if povprecje == 0
    P = 0;
  else
    P = poisspdf(abs(sprememba), povprecje);
  endif
endfunction

function A = povprecje(S)
  global N
  A = sum((N:-1:0)' .* S);
endfunction

function D = odmik(S, A)
  global N
  MS = sum((N:-1:0)' .* (N:-1:0)' .* S);
  D = sqrt(MS - A^2);
endfunction

W = zeros(M);
for j=1:M
  s = dt * smrtnost * (M-j);
  for i=1:M
    if i > j
      W(i,j) = verjetnost(i-j, s);
    else
      W(i,j) = 0;
    endif
  endfor
  W(j,j) = 1 - sum(W(1:M,j));
endfor

W = W^K;

S = [1; zeros(N,1)];
R = [0, 0, 1, 0];
v = 0;
T = dt * K;
for i=1:100
  S = W * S;
  A = povprecje(S);
  R = [R; i * T, S(M)-v, A/N, odmik(S,A)/sqrt(N)];
  v = S(M);
endfor

Ostanek = 1-v
Normalizacija = sum(S)

save "g_matrike_exp_250.dat" R


