function [m, M, data] = normalize(D)
    m = min(D);
    M = max(D);
    if (m == M)
        M = m + 1;
    endif
    data = (D - m)/(M-m);
endfunction

function [c, a, cov, k] = svdfit(S, y, wcrit)
  [U, W, V] = svd(S);

  [m n] = size(S);
  M = zeros(n, m);
  k = 0;
  for i = 1:n
    if (W(i,i) > wcrit)
      M(i,i) = 1/W(i,i);
      k = k + 1;
    else
      W(i,i) = 0;
    endif
  endfor

  a = V * M * U' * y;
  cov = V * M(1:n,1:n) * V';
  c = sumsq(S*a - y) / (m-k);
endfunction

function [ChiSqRed, a, k] = prevodnost ( wcrit )
  ## Preberemo podatke in jih reskaliramo, da bodo vsi med 0 in 1
  Data = dlmread("prevodnost.dat");
  [t0, t1, T] = normalize(Data(:,1));
  [p0, p1, P] = normalize(Data(:,2));
  Y = Data(:,3);
  [y0, y1, y] = normalize(Y);
  Sigma = Data(:,4) / ( y1 - y0 );
  C = diag(Sigma.^-1);

  ## Ustvarimo matriko S da bodo nasi parametri resitve sistema S*a = y
  N = 3;
  n = N^2;
  m = size(T,1);
  S = zeros(m,n);
  for i = 0:N-1
      for j = 0:N-1
	  k = j*N+i+1;
	  S(:,k) = T.^i .* P.^j;
      endfor
  endfor

  [ChiSqRed, a, k] = svdfit(C*S, C*y, wcrit);
endfunction


function [c, a, k, cov] = kristali(wcrit)
  E = dlmread('wrboost_eight.dat');
  [e0, e1, E(:,2)] = normalize(E(:,2));
  O = dlmread('wrboost_omega.dat');
  [o0, o1, O(:,2)] = normalize(O(:,2));
  T = dlmread('wrboost_theta.dat');
  [t0, t1, T(:,2)] = normalize(O(:,2));

  m = size(E,1);
  D = E(:,2);
  n = 15;
  S = zeros(m, n);
  for i = 0:n-1
      S(:,i+1) = D.^i;
  endfor

  Sigma = ones(m,1) * 2 / (e1 - e0);
  C = diag(Sigma .^ -1);

  [c, a, k, cov] = svdfit(C*S, C*D, wcrit);
endfunction

[c, a, cov, k] = kristali(0.04);
c, k, [a diag(cov)]
