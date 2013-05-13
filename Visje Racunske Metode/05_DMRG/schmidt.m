function E = entanglement_entropy(psi, d, nA)
    N = length(psi);
    n = round(log(N) / log(d));
    
    if (d^n != N)
        error("Size of psi must be d^n");
    endif
    
    nB = n - nA;
    
    # TODO: To dela samo, ce je prvih nA spinov v A, ostali pa v B
    # TODO: Ugotovi kako je treba reshape()-at, da so lahko poljubni spini v A in ostali v B
    Psi = reshape(psi, d^nA, d^nB);
    
    [U, S, V] = svd(Psi);
    lmu = diag(S) .^ 2;
    
    E = -sum(lmu .* log(lmu));
endfunction

function Psi = reshape_noncompact(psi, d, A)
    N = length(psi);
    n = round(log(N) / log(d));

    nA = length(A);
    nB = n - nA;
    Psi = zeros(d^nA, d^nB);
    
    for i = 1:N
      bits = dec2base(i-1, d, n)
      row = base2dec(bits(A), d)+1
      bits(A) = [q];
      column = base2dec(bits, d)+1
      Psi(row, column) = psi(i);
    endfor
endfunction

function E = entanglement_entropy_noncompact(psi, d, A)
    Psi = reshape_noncompact(psi, d, A);
    
    [U, S, V] = svd(Psi);
    lmu = diag(S) .^ 2;
    
    E = -sum(lmu .* log(lmu));
endfunction

function test_entent_random()
    d = 2;
    n = 12;
    N = d^n;
    
    nA_array = 1:n-1;
    E_array = [];
    Sig_array = [];
    
    for nA = 1:n-1
      v = [];
      for i = 1:20
        psi = rand(N, 1);
        psi /= norm(psi);
        E = entanglement_entropy(psi, d, nA);
        v = [v; E];
      endfor
      E_array = [E_array; mean(v)];
      Sig_array = [Sig_array; std(v)];
    endfor
    
    R = [nA_array' E_array Sig_array];
    dlmwrite("g_entent_random.dat", R);
endfunction

function test_entent_ground()
    d = 2;
    n = 12;
    N = d^n;
    
    nA_array = 1:n-1;
    E_array = [];
    I_array = [];
    
    for nA = nA_array
        psi = ground_state(n, d);
        E = entanglement_entropy(psi, d, nA);
        E_array = [E_array; E];
    endfor
    R = [nA_array', E_array];
    dlmwrite("g_entent_ground.dat", R);
endfunction

function UN = napihni(U, n, d, j)
  A = {};
  if (j == n)
    UN = spalloc(d^n, d^n, d^(n-2)*nnz(U));
    offset = d^(n-1);
    
    A = U(1:d,1:d);
    B = U(1:d,d+1:d+d);
    C = U(d+1:d+d,1:d);
    D = U(d+1:d+d,d+1:d+d);
    
    T = spalloc(1,1,1);
    T(1,1) = 1;
    F = kron(speye(d^(n-2)), T);
    
    UN = [kron(F,A), kron(F,B); kron(F,C), kron(F,D)];
  else
    block = kron(U, speye(d^(j-1)));
    UN = kron(speye(d^(n-j-1)), block);
  endif
  
  assert(nnz(UN) == d^(n-2)*nnz(U));
  assert(size(UN) == [d^n, d^n]);
endfunction

function H = hamiltonian(n, d, pbc)
    if (d == 2)
        h = spalloc(4, 4, 6);
        h(1,1) = h(4,4) = 1;
        h(2,2) = h(3,3) = -1;
        h(2,3) = h(3,2) = 2;
    else
        # TODO: Matrike za visje spine
        error("Only d=2 supported now")
    endif
    
    N = d^n;
    H = spalloc(N, N, 2*N);
    for j = 1:n-1
        H += napihni(h, n, d, j);
    endfor
    if (pbc)
        H += napihni(h, n, d, n);
    endif
    H = sparse(H);
endfunction

function psi = ground_state(n, d, pbc = 1)
    H = hamiltonian(n, d, pbc);
    [V, D] = eigs(H, 6, "sa");
    psi = V(:,6);
endfunction

function S = neumann_entropy(rho)
    S = -trace(rho * logm(rho));
endfunction

function I = qmi(psi, d, nA)
    N = length(psi);
    n = round(log(N) / log(d));
    nB = n - nA;
    Psi = reshape(psi, d^nA, d^nB);
    
    [U, S, V] = svd(Psi);
    
    rhoA = U*S*S'*U';
    rhoB = V*S'*S*V';
    rho = psi * psi';
    
    I = neumann_entropy(rhoA) + neumann_entropy(rhoB) - neumann_entropy(rho);
endfunction

function velikost_sistema()
  d = 2;
  
  nn = [8 10 12 14];
  R = [];
  
  for n = nn
    N = d^n;
    nA = n/2;
    
    psi = ground_state(n, d);
    E = entanglement_entropy(psi, d, 1);
    En = entanglement_entropy(psi, d, n/2);
    R = [R; n, E, En];
  endfor

  dlmwrite("g_entent_n.dat", R);
endfunction

function nekompaktna()
  d = 2;
  n = 12;
  
  R = [];
  psi = ground_state(n, d);
  for p = [2 3 4 6]
    A = 1:p:N;
    E = entanglement_entropy_noncompact(psi, d, A);
    R = [R; p, E];
  endfor
  
  dlmwrite("g_entent_nonc.dat", R);
endfunction