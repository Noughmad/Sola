function A = mpa(psi, d)
    N, o = size(psi);
    if (N != 1 && o != 1)
        error("psi must be a vector");
    endif
    if (N == 1 && o != 1)
        psi = psi';
        N, o = size(psi);
    endif

    n = log(N) / log(d);
    if (floor(x) != x)
        error("psi must be 2^d - dimensional");
    endif

    A = {};
    psi = reshape(psi, d, d^(n-1));
    U, S, V = svd(psi);
    for i = 1:d
        A{1, i} = U();
    endfor
    psi = neki(S, V);

endfunction