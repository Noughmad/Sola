global UU
global UL
global US

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
  global UL
  global US

  UU = {};
  
  for i = 1:length(shema)  
    UL{i} = speye(2^n);
    US{i} = speye(2^n);
    for j = 1:n
      UU{i}{j} = dvodelcni_cel(j, n, shema(i));
      if mod(j,2)
        UL{i} = UL{i} * UU{i}{j};
      else
        US{i} = US{i} * UU{i}{j};
      endif
    endfor
  endfor
endfunction

function psiprime = lihi_korak(psi, n, i)
  global UU
  global UL
  
  # psiprime = UL{i} * psi;
  # return

  psiprime = psi;
  for j = 1:2:n-1
    U = UU{i}{j};
    psiprime = U * psiprime;
  endfor
endfunction

function psiprime = sodi_korak(psi, n, i)
  global UU
  global US

  # psiprime = US{i} * psi;
  # return

  psiprime = psi;
  for j = 2:2:n
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

function s = shema_s5()
  p1 = 1/4 + i*sqrt(3)/12;
  p2 = 2*p1;
  p3 = 1/2;
  
  s = [p1, p2, p3, p2', p2', p2', p3', p2, p1] / 2;
endfunction

function H = hamiltonian(n, j)
  h = spalloc(4, 4, 6);
  h(1,1) = h(4,4) = 1;
  h(2,2) = h(3,3) = -1;
  h(2,3) = h(3,2) = 2;
  
  H = napihni(h, n, j);
endfunction

function F, E = povprecna_energija(n, Npsi, Nkorakov, beta)
  shema = shema_s5();
  z = beta / 2 / Nkorakov;
  make_matrices(n, shema .* z)
  s = length(shema);
  N = 2^n;
  
  H = spalloc(N, N, 2*N);
  for j = 1:n
    H += hamiltonian(n, j);
  endfor
  H = sparse(H);
  
  psi = rand(N,Npsi);
  psi = psi ./ norm(psi,2,"columns");

  for k = 1:Nkorakov
    psi = korak(psi, n, s);
  endfor
  
  Z = diag(psi' * psi)';
  E = real(diag(psi' * H * psi))';
  
  Zp = sum(Z) / Npsi;
  F = -log(Z) / beta;
  E = E ./ Z;
endfunction

function R = prosta_energija_vse()
  n = 10;
  Npsi = 50;
  N = 2^n;
  R = [];
  for beta = [1e-4, 3e-4, 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100]
    Nkorakov = max(10, 1000*beta);
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

function T = casovna_odvisnost(name)
  T = [];
  for n = [6 8 10 12 14 16]
    t = time();
    povprecna_energija(n, 10, 10, 1);
    T = [T; n, time() - t];
  endfor
  
  dlmwrite(["g_time_" name ".dat"], T);
endfunction

function M = magnetizacija_mat(n)
  S = sparse([1, 0; 0 -1]);
  M = kron(speye(2^(n-1)), S);
endfunction

function S,D = mag_korelacija(n, Npsi, Nkorakov, Nizmerkov, t)
  
  shema = shema_s5();
  z = -i * t / Nkorakov;
  make_matrices(n, shema .* z)
  s = length(shema);
  N = 2^n;
  
  M = magnetizacija_mat(n);
  C = [];
  
  psi = rand(N,Npsi);
  psi = psi ./ norm(psi,2,"columns");
  chi = M * psi;
  
  C = [C; real(diag(psi' * M * chi))'];
  
  for i = 1:Nizmerkov
    for k = 1:Nkorakov
      psi = korak(psi, n, s);
      chi = korak(chi, n, s);
    endfor
    
    C = [C; real(diag(psi' * M * chi))'];
    # dlmwrite("g_korelacija_raw.dat", C);
  endfor
  
  T = (0:Nizmerkov) * t;
  
  S = sum(C,2) / Npsi;
  D = sqrt(sumsq(C-S,2) / Npsi);
  dlmwrite("g_korelacija.dat", [T', S, D]);
endfunction

function J = tok_mat_d(n, j)
  sx = [0, 1; 1, 0];
  sy = [0, -i; i, 0];
  Jd = kron(sx, sy) - kron(sy, sx);
  J = napihni(Jd, n, j);
endfunction

function J = tok_mat(n)
  J = spalloc(2^n, 2^n, 2*2^n);
  for j = 1:n
    J += tok_mat_d(n, j);
  endfor
endfunction

function S,D = tok_korelacija(n, Npsi, Nkorakov, Nizmerkov, t)
  
  shema = shema_s4();
  z = -i * t / Nkorakov;
  make_matrices(n, shema .* z)
  s = length(shema);
  N = 2^n;
  
  M = tok_mat(n);
  C = [];
  
  psi = rand(N,Npsi);
  psi = psi ./ norm(psi,2,"columns");
  chi = tok_mat_d(n, 1) * psi;
  
  C = [C; real(diag(psi' * M * chi))'];
  
  for i = 1:Nizmerkov
    for k = 1:Nkorakov
      psi = korak(psi, n, s);
      chi = korak(chi, n, s);
    endfor
    
    C = [C; real(diag(psi' * M * chi))'];
    # dlmwrite("g_tok_raw.dat", C);
  endfor
  
  T = (0:Nizmerkov) * t;
  
  S = sum(C,2) / Npsi;
  D = sqrt(sumsq(C-S,2) / Npsi);
  dlmwrite("g_tok_long.dat", [T', S, D]);
endfunction

function C = komutator(A,B)
  C = A*B - B*A;
endfunction

function F = four_trans(R)
  C = R(:,2);
  L = length(C);
  t = R(2,1) * L;
  
  F = [(0:L-1)' / t, abs(fft(C))];
endfunction