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

function test_entent_random()
    d = 2;
    n = 12;
    N = d^n;
    
    nA_array = [];
    E_array = [];
    
    for i = 1:100
        nA = randi(n-1);
        psi = rand(N, 1);
        psi /= norm(psi);
        E = entanglement_entropy(psi, d, nA);
        nA_array = [nA_array; nA];
        E_array = [E_array; E];
    endfor
    
    [nA_array, I] = sort(nA_array);
    E_array = E_array(I);
    
    plot(nA_array, E_array);
endfunction

function test_entent_ground()
    d = 2;
    n = 12;
    N = d^n;
    
    nA_array = [];
    E_array = [];
    
    for nA = 1:n-1
        psi = ground_state(n, d);
        E = entanglement_entropy(psi, d, nA);
        nA_array = [nA_array; nA];
        E_array = [E_array; E];
    endfor
    
    plot(nA_array, E_array);
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