N = 5;
M = N+1;

smrtnost = 1.0;
dt = 0.01;

function P = verjetnost(sprememba, povprecje)
  if povprecje == 0
    P = 0;
  else
    P = poisspdf(abs(sprememba), povprecje);
  endif
endfunction

W = zeros(M);
for j=1:M
  s = dt * smrtnost * j;
  for i=1:M
    if i > j
      W(i,j) = verjetnost(i-j, s);
    elseif i < j
      W(i,j) = 0;
    else
      W(i,j) = 0;
    endif
  endfor
  W(j,j) = 1 - sum(W(1:M,j));
endfor

W = W^100

S = [1; zeros(N,1)];
for i=1:10
  S = W * S;
  S'
endfor


