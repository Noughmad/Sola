global UU

function M = block_diagonal(block, m)
  [a,b] = size(block);
  M = spalloc(m*a,m*b,m*nnz(block));
  for i=1:m
    M((i-1)*a+1:i*a,(i-1)*b+1:i*b) = block;
  endfor
endfunction

function UN = napihni(U, n, j)
  A = {};
  if (j == n)
    UN = spalloc(2^n, 2^n, 2^(n-2)*nnz(U));
    offset = 2^(n-1);
    
    A = U(1:2,1:2);
    B = U(1:2,3:4);
    C = U(3:4,1:2);
    D = U(3:4,3:4);
    
    T = spalloc(1,1,1);
    T(1,1) = 1;
    F = kron(speye(2^(n-2)), T);
    
    UN = [kron(F,A), kron(F,B); kron(F,C), kron(F,D)];
  else
    block = kron(U, speye(2^(j-1)));
    UN = kron(speye(2^(n-j-1)), block);
  endif
  
  assert(nnz(UN) == 2^(n-2)*nnz(U));
  assert(size(UN) == [2^n, 2^n]);
endfunction

function S = base_state(I, n)
  S = dec2bin(I-1, n);
endfunction

function I = state_index(S)
  I = bin2dec(S);
endfunction

function U = dvodelcni(z)
  U = spalloc(4, 4, 6);
  U(1,1) = exp(2*z);
  U(4,4) = exp(2*z);
  U(2,2) = U(3,3) = cosh(2*z);
  U(2,3) = U(3,2) = sinh(2*z);
  U = exp(-z) * U;
endfunction

function Uc = dvodelcni_cel(j, n, z)
  Uc = napihni(dvodelcni(z), n, j);
endfunction

function make_matrices(n, shema)
  global UU

  UU = {};
  
  for i = 1:length(shema)  
    for j = 1:n
      UU{i}{j} = dvodelcni_cel(j, n, shema(i));
    endfor
  endfor
endfunction

function psiprime = lihi_korak(psi, n, i)
  global UU

  psiprime = psi;
  for j = 1:2:n-1
    U = UU{i}{j};
    psiprime = U * psiprime;
  endfor
endfunction

function psiprime = sodi_korak(psi, n, i)
  global UU

  psiprime = psi;
  for j = 2:2:n-1
    U = UU{i}{j};
    psiprime = U * psiprime;
  endfor
endfunction

function psiprime = korak(psi, n, s)
  psiprime = psi;
  lih = 1;
  for i = 1:s
    if lih
      psiprime = lihi_korak(psiprime, n, i);
    else
      psiprime = sodi_korak(psiprime, n, i);
    endif
    lih = 1-lih;
  endfor
endfunction

function s = shema_s4()
  r = 2^(1/3);
  x0 = -r / (2 - r);
  x1 = 1 / (2 - r);
  
  s = [x1/2 x1 (x1+x0)/2 x0 (x1+x0)/2 x1 x1/2];
endfunction

function H = hamiltonian(n, j)
  h = spalloc(4, 4, 6);
  h(1,1) = h(4,4) = 1;
  h(2,2) = h(3,3) = -1;
  h(2,3) = h(3,2) = 2;
  
  H = napihni(h, n, j);
endfunction

function F, E = povprecna_energija(n, Npsi, Nkorakov, beta)
  Z = zeros(1,Npsi);
  E = zeros(1,Npsi);
  
  shema = shema_s4();
  z = beta / 2 / Nkorakov;
  make_matrices(n, shema .* z)
  s = length(shema);
  N = 2^n;
  
  H = spalloc(N, N, 4*N);
  for j = 1:n-1
    H += hamiltonian(n, j);
  endfor
  H = sparse(H);
  
  for p = 1:Npsi
    psi = rand(N,1);
    psi = psi / sqrt(psi' * psi);

    for k = 1:Nkorakov
      psi = korak(psi, n, s);
    endfor
  
    Z(p) = psi' * psi;
    E(p) = psi' * H * psi / Z(p);
  endfor
  
  F = -log(Z) / beta;
endfunction

function R = prosta_energija_vse()
  n = 16;
  Npsi = 20;
  N = 2^n;
  R = [];
  for beta = [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30]
    Nkorakov = max(10, 1000 * beta);
    [F, E] = povprecna_energija(n, Npsi, Nkorakov, beta);
    avgF = sum(F) / Npsi;
    stdevF = sqrt(sumsq(F - avgF) / Npsi);
    avgE = sum(E) / Npsi;
    stdevE = sqrt(sumsq(E - avgE) / Npsi);
    
    R = [R; beta, avgF, stdevF, avgE, stdevE];
  endfor
  dlmwrite("g_FE.dat", R);
endfunction

function vse_profile()
  profile off
  profile on
  prosta_energija_vse();
  T = profile('info');
  profile off
  
  profshow(T)
endfunction

function T = casovna_odvisnost()
  T = [];
  for n = [6 8 10 12 14 16]
    t = time();
    povprecna_energija(n, 10, 10, 1);
    T = [T; n, time() - t];
  endfor
  
  dlmwrite("g_time.dat", T);
endfunction
