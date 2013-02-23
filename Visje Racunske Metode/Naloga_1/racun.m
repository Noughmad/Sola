global L
global N
global tau
global h

# Precej v redu vrednosti so L=20, N=400, potem je lahko |z| nekje do 5

L = 20;
N = 400;
h = 2*L/N;
tau = h^2

function p = hermitovp(n)
  if n == 0
    p = [1];
  elseif n == 1
    p = [2, 0];
  else
    m = hermitovp(n-1);
    p = [2*m, 0] - [0, 0, polyder(m)];
  endif
endfunction

function x = space()
  global L;
  global N;
  x = linspace(-L, L, N)';
endfunction

function psi = stanje(n)
  psi = koherentno_stanje(n, 0);
endfunction

function psi = stanje(n)
  global N;
  global L;
  x = space();
  psi = polyval(hermitovp(n), x) .* exp(-0.5 * x.^2) / pi^(1/4) / sqrt(2^n * factorial(n));
endfunction

function psi = koherentno_stanje(z)
  global N;
  psi = zeros(N, 1);
  C = exp(-abs(z^2)/2);
  for n = 0:30
    psi = psi + C * z^n / sqrt(factorial(n)) * stanje(n);
  endfor
endfunction

function H = hamiltonian(lambda)
  global N;
  global h;
  global L;

  H = zeros(N,N);
  H = H + diag(-2*ones(N,1));
  H = H + diag(ones(N-1,1),-1);
  H = H + diag(ones(N-1,1),1);
  H = -H / h^2 / 2;
  
  x = space();
  H = H + diag(0.5 * x.^2 + lambda * x.^4);
endfunction

function U,V = implicitna(H)
  global tau;
  U = eye(size(H)) + i*tau/2*H;
  V = eye(size(H)) - i*tau/2*H;
endfunction

function U = eksplicitna(H)
  global tau;

  U = zeros(size(H));
  F = eye(size(H));
  for k = 1:10
    U = U + F;
    F = F * (-i*tau)/k * H;
  endfor
endfunction

function n = norma(psi)
  global h;
  n = psi' * psi * h;
endfunction

function E = energija(H, psi)
  global h;
  E = (psi' * H * psi * h) / (norma(psi));
endfunction

function primerjava(lambda, z)
  H = hamiltonian(lambda);
  psi_e = koherentno_stanje(z);
  psi_i = psi_e;
  U = eksplicitna(H);
  [W, V] = implicitna(H);
  x = space();
  
  for step = 1:1000
    psi_e = U * psi_e;
    psi_i = W \ (V * psi_i);
    plot(x, abs(psi_e), x, abs(psi_i));
    usleep(10);
  endfor
endfunction