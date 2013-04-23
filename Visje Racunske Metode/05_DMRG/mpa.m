warning ("off", "Octave:broadcast");

function i = psi_index(S, d)
    i = 1;
    j = 0;
    for s = S'
        i = i + (s-1)*d^j;
        j = j + 1;
    endfor
endfunction

function A = mpa(psi, d)
    [N, o] = size(psi);
    if (N != 1 && o != 1)
        error("psi must be a vector");
    endif
    if (N == 1 && o != 1)
        psi = psi';
        N, o = size(psi);
    endif

    n = round(log(N) / log(d))
    if (d^n != N)
        N,n
        error("psi must be d^n - dimensional");
    endif
    

    A = {};
    M = {};
    
    psi = reshape(psi, d, d^(n-1));
    [U, S, V] = svd(psi);
    M{1} = d;
    for s = 1:d
        A{1, s} = U(s,:);
    endfor
    psi = diag(S) .* V'(1:M{1},:);
    
    for j = 2:n-1
        psi = reshape(psi, M{j-1}*d, d^(n-j));
        [U, S, V] = svd(psi);
        m = min(size(S));
        M{j} = min(size(S));
        
        for s = 1:d
            A{j, s} = U(M{j-1}*(s-1)+1:M{j-1}*s,1:m);
        endfor

        psi = diag(S) .* V'(1:M{j},:);
    endfor
    
    psi = reshape(psi, M{n-1}, d);
    for s = 1:d
        A{n, s} = psi(:,s);
    endfor
endfunction

function p = mpa_psi_element(A, S)
    [n, d] = size(A);
    p = 1;
    for j = 1:n
        p = p * A{j, S(j)};
    endfor
endfunction

function Psi, Mpa = test_mpa(d, n)
    psi = rand(d^n, 1);
    A = mpa(psi, d);
    
    S = randi(d, n, 1);
    i = psi_index(S, d);
    
    Psi = psi(i);
    Mpa = mpa_psi_element(A, S);
endfunction

function test_many()
    for n = [2 6 10 14]
        for i = 1:10
            [p, m] = test_mpa(2, n);
            NumError = log(abs(p/m - 1));
            assert(NumError < -16);
        endfor
    endfor
    
    for d = 2:6
        for i = 1:10
            [p, m] = test_mpa(d, 5);
            NumError = log(abs(p/m - 1));
            assert(NumError < -16);
        endfor
    endfor
endfunction