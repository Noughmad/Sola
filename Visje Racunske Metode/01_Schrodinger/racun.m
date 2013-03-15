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
# Za 2D: Treba manjse, na primer L=6 in N=100, |z| do 2
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

function psi = hermitov(n, x)
  m = 0:floor(n/2);
  psi = zeros(size(x));
  for i=1:size(x,1)
    v = (-1).^m .* bincoeff(n, m) .* factorial(n-m) ./ factorial(n-2*m) .* (2*x(i)).^(n-2*m);
    psi(i) = sum(v);
  endfor
endfunction


function psi = stanje(n)
  x = space();
  psi = hermitov(n, x) .* exp(-0.5 * x.^2) / pi^(1/4) / sqrt(2^n * factorial(n));
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
  
  csvwrite(["data/g_implicitna_" num2str(lambda) ".csv"], Pr_i);
  csvwrite(["data/g_eksplicitna_" num2str(lambda) ".csv"], Pr_e);
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
    print("-dpng", "-S500,480", sprintf("data/g_animation_2D_%data/g_%.3d.png", lambda, step));
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

function racun(lambda, z, S, M)
  global N;
  H = hamiltonian(lambda);
  [W, V] = implicitna(H);
  psi = koherentno_stanje(z);
  x = space();
  
  for step = 1:300
    for i = 1:S
      psi = W \ (V * psi);
    endfor
    plot(x, abs(psi));
    axis([-10, 10, 0, M]);
    print("-dpng", "-S500,480", sprintf("data/g_animation_1D_%data/g_%.3d.png", lambda, step));
  endfor
endfunction

function nihanje()
# animacija nihanja v 1D
  set_sizes(10, 400)
  racun(0, 3, 20, 1);
  racun(0.1, 3, 10, 2);
  racun(1, 3, 5, 2);
endfunction
  

# stabilnost()
# krozenje(0)
# krozenje(0.1)
# krozenje(1)

global I
global V

function p = numerov_end(E)
  global N
  global h
  global V
  global I
  
  p1 = I(1);
  p2 = I(2);
  
  K = 2 * (E - V);

  hh = h*h;

  for i=2:N-1
    t = ( 2*(1-5*hh/12*K(i))*p2 - (1+hh/12*K(i-1))*p1 ) / (1+hh/12*K(i+1));
    p1 = p2;
    p2 = t;
  endfor
  p = t;
endfunction

function En = numerov(n, lambda)
  global N;
  global L;
  global h;
  global I;
  global V;

  E = n + 1/2;
  x = linspace(0, L, N);
  V = lambda * x .* x .* x .* x;
  h = L/N;
  
  if rem(n, 2) == 0
    I = [1, 0];
  else
    I = [0, 1];
  endif
  
  En = fsolve(@numerov_end, E)
endfunction

function p = mat_element(psi1, A, psi2)
  global h;
  s = size(A);
  if s(1) == 1 || s(2) == 1
    A = diag(A);
  endif
  p = real(psi1' * A * psi2 * h) / sqrt(norma(psi1) * norma(psi2));
endfunction

function p = mat_element_c(psi1, A, psi2)
  global h;
  s = size(A);
  if s(1) == 1 || s(2) == 1
    A = diag(A);
  endif
  p = (psi1' * A * psi2 * h) / sqrt(norma(psi1) * norma(psi2));
endfunction

function Hb = ham_matrika(N, lambda)
  F = lambda / 4;
  D = zeros(N,5);
  n = 0:N-1;
  D(:,1) = F * sqrt(n.*(n-1).*(n-2).*(n-3));
  D(:,2) = F * (2*(2*n-1) .* sqrt(n .* (n-1)));
  D(:,3) = n + 0.5 + F * (1 + 4*n + 4*n.^2 + (n+1).*(n+2) + n .* (n-1));
  D(:,4) = F * 2*(2*n+3) .* sqrt((n+1) .* (n+2));
  D(:,5) = F * sqrt((n+1).*(n+2).*(n+3).*(n+4));
  Hb = spdiags(D, [4, 2, 0, -2, -4], zeros(N));
endfunction

function Hb = mat_v_bazi(baza, lambda)
  n = size(baza, 2);
  Hb = zeros(n);
  H = hamiltonian(0);
  x = space();
  Vd = lambda * x .^ 4;
  
  for i=1:n
    for j=1:n
      Hb(i,j) = mat_element(baza(:,i), Vd, baza(:,j));
    endfor
    Hb(i,i) = Hb(i,i) + i - 0.5;
  endfor
endfunction

function Hb = lan_matrika(n, lambda)
  baza = lanczos_baza(n, lambda);
  Hb = spalloc(n, n, 3*n-2);
  
  H = hamiltonian(lambda);

  for i=1:n
    if i < n
      Hb(i,i+1) = mat_element(baza(:,i), H, baza(:,i+1));
    endif
    Hb(i,i) = mat_element(baza(:,i), H, baza(:,i));
    if i > 1
      Hb(i,i-1) = mat_element(baza(:,i), H, baza(:,i-1));
    endif
  endfor

endfunction

function Em = matricni(Hb, st)
  n = size(Hb,1);
  Em = eigs(Hb, st, "sm");
  if Em(1) > Em(st)
      Em = flipud(Em);
  endif
endfunction

function st = st_skonvergiranih(E)
  n = size(E,1)
  st = zeros(n, 1);
  for i = 1:n
    st(i) = sum(abs(E(i,:) - E(n,:)) < 1e-3);
  endfor
endfunction

function konvergenca_x(lambda)
  E = [];
  set_sizes(20, 400);
  Bb = 25:5:175;
  for b=Bb
    E = [E; matricni(mat_v_bazi(ho_baza(b), lambda), 20)'];
  endfor
  
  E = [Bb', st_skonvergiranih(E), E];
  csvwrite(["data/g_konvergenca_x_" num2str(lambda) ".csv"], E);
endfunction

function konvergenca_ho(lambda)
  E = [];
  Bb = [25:5:225, 1000];
  for b=Bb
    E = [E; matricni(ham_matrika(b, lambda), 20)'];
  endfor
  
  E = [Bb', st_skonvergiranih(E), E];
  csvwrite(["data/g_konvergenca_ho_" num2str(lambda) ".csv"], E);
endfunction

function konvergenca_L(lambda)
  set_sizes(20, 2000);
  E = [];
  Bb = 25:10:275;
  for b=Bb
    E = [E; matricni(lan_matrika(b, lambda), 20)'];
  endfor
  
  E = [Bb', st_skonvergiranih(E), E];
  csvwrite(["data/g_konvergenca_L_" num2str(lambda) ".csv"], E);
endfunction

function vse_konvergence()
  for lambda=[0.1, 1]
    konvergenca_ho(lambda);
    # konvergenca_x(lambda);
    konvergenca_L(lambda);
  endfor
endfunction

function B = ho_baza(n)
  B = [];
  for i=1:n
    B = [B, stanje(i-1)];
  endfor
endfunction

function L = lanczos_baza(n, lambda)
  global N
  
  H = hamiltonian(lambda);
  L = [];
  psi = stanje(0);
  psi = psi / sqrt(norma(psi));
  L = psi;
  
  psi = H*psi;
  psi = psi / sqrt(norma(psi));
  L = [L, psi];
  
  for i = 3:n
    pj = L(:,i-1);
    pm = L(:,i-2);
    psi = H*pj - pj * pricakovana(H, pj) - pm * mat_element(pm, H, pj);
    psi = psi / sqrt(norma(psi));
    L = [L, psi];
  endfor
endfunction

function En = lastne_energije(n, N)
  En = (0:n-1)';
  for lambda = [0 1e-4 3e-4 1e-3 3e-3 1e-2 3e-2 1e-1 3e-1 1]
    H = ham_matrika(N, lambda);
    En = [En, matricni(H, n)];
  endfor

  csvwrite(["data/g_energije.csv"], En);
endfunction

function psi = stanje_ob_casu(baza, energije, stanje, t)
  [N, M] = size(baza);
  psi = zeros(N,1);
  for i = 1:M
    psi = psi + baza(:,i) .* stanje(i) .* exp(j*energije(i)*t);
  endfor
endfunction

function razvoj(z, lambda)
  global tau;
  global N;
  M = 40;
  set_sizes(10, 400);
  
  [V, D] = eig(ham_matrika(M, lambda));
  [S, I] = sort(diag(D));
  energije = S;
  V = V(:,I);
  
  
  psi = koherentno_stanje(z);
  
  bazaHo = ho_baza(M);
  enHo = (1:M) - 0.5;
  baza = zeros(N, M);
  S = zeros(M,1);
  for i = 1:M
    P = stanje_ob_casu(bazaHo, enHo, V(:,i), 0);
    P = P / sqrt(norma(P));
    baza(:,i) = P;
    S(i) = mat_element_c(psi, eye(N), P);
  endfor
  
    
  x = space();
  H = hamiltonian(lambda);
  [W, V] = implicitna(H);
  
  # psi = stanje_ob_casu(baza, energije, S, 0);

  T = 10;
  for step=1:1000
    for i = 1:T
      psi = W \ (V * psi);
    endfor
    plot(x, abs(stanje_ob_casu(baza, energije, S, step * T * tau)), x, abs(psi));
    legend("Spektralna", "Direktna");
    axis([-10, 10, -1, 1]);
    print("-dpng", "-S500,480", sprintf("anim/g_razvoj_%g_%.3d.png", lambda, step));
  endfor
endfunction
