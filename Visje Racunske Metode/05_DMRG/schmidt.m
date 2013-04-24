function E = entanglement_entropy(psi, d, nA)
    N = length(psi);
    n = round(log(N) / log(d));
    
    if (d^n != N)
        error("Size of psi must be d^n");
    endif
    
    nB = n - nA;
    
    # TODO: To dela samo, ce je prvih nA spinov v A, ostali pa v B
    # TODO: Ugotovi kako je treba reshape()-at, da so lahko poljubni spini v A in ostali v B
    Psi = reshape(d^nA, d^nB);
    
    [U, S, V] = svd(Psi);
    d = diag(S) .^ 2;
    
    E = -sum(d .* log(d));
endfunction

