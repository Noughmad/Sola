global alpha

function Q = cena(chi, k)
  global alpha
  Q = (chi + 1/chi) + k*alpha;
endfunction

function [m, M, data] = normalize(D)
    m = min(D);
    M = max(D);
    if (m == M)
        M = m + 1;
    endif
    data = (2*D - M - m)/(M-m);
endfunction

function [c, a, cov] = svdfit(S, y, wcrit)
  [U, W, V] = svd(S);

  [m n] = size(S);
  M = zeros(n, m);
  p = min([m n]);
  for i = 1:p
    if (W(i,i) > wcrit)
      M(i,i) = 1/W(i,i);
    else
      W(i,i) = 0;
    endif
  endfor

  a = V * M * U' * y;  
  c = sumsq(S*a - y);

  if (nargout > 2)
    cov = V * M(1:n,1:n) * V';
  endif
endfunction

function [c, a, cov, k] = prevodnost ( wcrit, shrani )
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

  [c, a, cov, k] = svdbrute(C*S, C*y, wcrit);
  if (shrani)
    B = [y S*a];
    global alpha
    save (["prev_" int2str(10*alpha) ".dat"], "B");
    [a sqrt(diag(cov))]
  endif
endfunction

function a = kristal(name)
  P = dlmread(["wrboost_" name ".dat"]);
  [x0, x1, X] = normalize(P(:,1));
  [y0, y1, Y] = normalize(P(:,2));
  m = length(X);
  n = 16;

  S = zeros(m, n);
  for i = 0:n-1
      S(:,i+1) = X.^i;
  endfor

  Sigma = ones(m,1) * 2 / (y1 - y0);
  C = diag(Sigma .^ -1);

  [c, a, cov, k] = svdbrute(C*S, C*Y);
  name, c, [a sqrt(diag(cov))]

  B = [X Y S*a];
  save (["kris_" name ".dat"], "B");
endfunction

function [c, a, cov, k] = svdbrute(S, y)
  bestchi = Inf;
  [m,N] = size(S);
  u = 1:N;
  bestk = Inf;
  a = zeros(N,1);
  cov = zeros(N);
  for n = 1:min(N, 10)
    for used = nchoosek(1:N, n)'
      [chi, ta] = svdfit(S(:,used), y, 1e-5);
      if (cena(chi, n) < cena(bestchi, bestk))
	bestchi = chi/(m-n);
	bestk = n;
	a = zeros(N,1);
	cov = zeros(N);
	[chi, a(used), cov(used,used)] = svdfit(S(:,used), y, 1e-5);
	k = n;
      endif
    endfor
  endfor
  c = bestchi;
  k = bestk;
endfunction

function preva()
  global alpha
  A = [];
  for p = linspace(0,2,100)
    alpha = p;
    [c, a, cov, k] = prevodnost(1e-5, 0);
    A = [A; alpha c k];
  endfor
  A
  save prev_a.dat A
endfunction

function kristali()
  K = [kristal("eight") kristal("omega") kristal("theta")]
  save kristal.dat K
endfunction

# preva()

alpha = 0.1;
# kristali()
alpha = 0;
kristal("skupaj");

alpha = 1.5;
prevodnost(1e-5, 1)
