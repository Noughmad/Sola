global L
global N
global tau
global h

# Precej v redu vrednosti so L=20, N=400, potem je lahko |z| nekje do 5

function set_sizes(l, n)
  global L
  global N
  global tau
  global h
  
  L = l;
  N = n;
  h = 2*L/N;
  tau = h^2
endfunction

# Za 1D: dobre vrednosti je L=20, N = 400, |z| do 10
# Za 2D: Treba manjse, naprimer L=6 in N=100, |z| do 2
set_sizes(20, 100)

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
  x = space();
  psi = polyval(hermitovp(n), x) .* exp(-0.5 * x.^2) / pi^(1/4) / sqrt(2^n * factorial(n));
endfunction

function psi = koherentno_stanje(z)
  n = 0;
  x = space() - real(z);
  psi = polyval(hermitovp(n), x) .* exp(-0.5 * x.^2 + i*imag(z)*x) / pi^(1/4) / sqrt(2^n * factorial(n));
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
  H = sparse(H);
endfunction

function H = ham2d(lambda)
  global N;
  global h;
  global L;
  
  H = -4*speye(N*N, N*N);

  H = spdiags(ones(N*N,4),[-N, -1, 1, N], H);
  H = -H / h^2 / 2;
  
  x = space();
  
  T = 0.5 * ones(N,1) * x'.^2 + 0.5 * x.^2 * ones(1,N) + lambda * x.^2 * x'.^2;
  H = H + spdiags(reshape(T, N*N, 1), 0, speye(N*N));
endfunction

function U,V = implicitna(H)
  global tau;
  U = sparse(eye(size(H)) + i*tau/2*H);
  V = sparse(eye(size(H)) - i*tau/2*H);
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
  n = real(psi' * psi * h);
endfunction

function p = pricakovana(A, psi)
  global h;
  s = size(A);
  if s(1) == 1 || s(2) == 1
    A = diag(A);
  endif
  p = real(psi' * A * psi * h) / norma(psi);
endfunction

function plot2d(psi)
  global N;

  x = space();
  [xx, yy] = meshgrid(x, x);
  contourf(xx, yy, reshape(psi, N, N));
endfunction

function s = sirina(A, psi)
  s = sqrt(pricakovana(A.^2, psi) - pricakovana(A, psi)^2);
endfunction

function primerjava(lambda, z)
  H = hamiltonian(lambda);
  psi_e = koherentno_stanje(z);
  psi_i = psi_e;
  U = eksplicitna(H);
  [W, V] = implicitna(H);
  x = space();
  
  Pr_i = [];
  Pr_e = [];
  
  for step = 1:10000
    psi_e = U * psi_e;
    psi_i = W \ (V * psi_i);
    # plot(x, abs(psi_i));
    Pr_i = [Pr_i; step, norma(psi_i), pricakovana(x, psi_i), pricakovana(H, psi_i), sirina(x, psi_i)];
    Pr_e = [Pr_e; step, norma(psi_e), pricakovana(x, psi_e), pricakovana(H, psi_e), sirina(x, psi_e)];
    # usleep(30);
  endfor
  
  csvwrite(["g_implicitna_" num2str(lambda) ".csv"], Pr_i);
  csvwrite(["g_eksplicitna_" num2str(lambda) ".csv"], Pr_e);
endfunction

function racun2d(lambda, a, b);
  global N;
  H = ham2d(lambda);
  [W, V] = implicitna(H);
  psi = reshape(koherentno_stanje(a) * koherentno_stanje(b)', N*N, 1);
  
  for step = 1:300
    for i = 1:3
      psi = W \ (V * psi);
    endfor
    plot2d(abs(psi));
    axis([-6, 6, -6, 6], "square")
    print("-dpng", "-S500,480", sprintf("g_animation_%g_%.3d.png", lambda, step));
  endfor
endfunction

function krozenje(lambda)
# Lepa animacija krozenja, ce je lambda 0
# Za razpade naj bo lambda med 0.1 (blizu krozenja, vidno razpadanje) in 1 (hitro razpade)
  set_sizes(6, 100);
  racun2d(lambda, 2, 2i);
endfunction

function stabilnost()
  set_sizes(20, 330)
  primerjava(0, 2)
  primerjava(0.1, 2)
  primerjava(1, 2)
endfunction

# stabilnost()
krozenje(0)
krozenje(0.1)
krozenje(1)