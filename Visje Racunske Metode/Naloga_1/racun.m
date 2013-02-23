global L
global N
global tau
global h

L = 50;
N = 1000;
h = 2*L/N;
tau = h^2

function x = space()
  global L;
  global N;
  x = linspace(-L, L, N)';
endfunction

function psi = stanje(n)
  psi = koherentno_stanje(n, 0);
endfunction

function psi = koherentno_stanje(n, z)
  global N;
  global L;
  x = space() - z;
  H = zeros(N, 1);
  if n == 0
    H = ones(size(x));
  elseif n == 1
    H = 2*x;
  elseif n == 2
    H = 4*x.^2 - 2;
  elseif n == 3
    H = 8*x.^3 - 12*x;
  endif
  psi = H .* exp(-0.5 * x.^2) / pi^(1/4) / sqrt(2^n * factorial(n));
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
  for k = 1:100
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

function rac_eksp(lambda, n, z)
  H = hamiltonian(lambda);
  psi = koherentno_stanje(n, z);
  U = eksplicitna(H);
  x = space();
  
  
  for step = 1:1000
    psi = U * psi;
    plot(x, real(psi), x, imag(psi));
    n = norma(psi)
    E = energija(H, psi)
    usleep(30);
  endfor
endfunction

function rac_imp(lambda, n, z)
  H = hamiltonian(lambda);
  psi = koherentno_stanje(n,z);
  [U,V] = implicitna(H);
  x = space();
  
  for step = 1:1000
    psi = U \ (V * psi);
    plot(x, real(psi), x, imag(psi));
    n = norma(psi);
    E = energija(H, psi);
    usleep(30);
  endfor
endfunction

# rac_eksp(0, 0);
# rac_imp(0, 0);