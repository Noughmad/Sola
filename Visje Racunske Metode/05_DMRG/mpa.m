warning ("off", "Octave:broadcast");

global MinSchmidt
global MaxM

MinSchmidt = 1e-3
MaxM = 40

function i = psi_index(S, d)
    i = 1;
    j = 0;
    for s = S'
        i = i + (s-1)*d^j;
        j = j + 1;
    endfor
endfunction

function [U, S, V, e] = svd_limited(A)
  global MinSchmidt
  global MaxM

  [U, S, V] = svds(A, MaxM);
  s = length(diag(S));
  l = sum(diag(S) > MinSchmidt);
  if s > l
    trunc = diag(S)(l+1:s);
    e = sumsq(trunc);
  else
    e = 0;
  endif
  U = U(:,1:l);
  S = S(1:l,1:l);
  V = V(:,1:l);
endfunction

function [MPA, truncation_error] = mpa(psi, d)
    truncation_error = 0

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
    L = {};
    
    psi = reshape(psi, d, d^(n-1));
    [U, S, V, e] = svd_limited(psi);
    truncation_error += e;
    M{1} = length(S);
    L{1} = S;
    for s = 1:d
        B{1, s} = U(s,:);
    endfor
    psi = diag(S) .* V';
    
    for j = 2:n-1
        psi = reshape(psi, M{j-1}*d, d^(n-j));
        [U, S, V, e] = svd_limited(psi);
        truncation_error += e;
        M{j} = length(S);
        
        L{j} = S;
        for s = 1:d
            B{j, s} = L{j-1} \ U(M{j-1}*(s-1)+1:M{j-1}*s,:);
        endfor

        psi = diag(S) .* V';
    endfor
    
    psi = reshape(psi, M{n-1}, d);
    for s = 1:d
        B{n, s} = L{n-1} \ psi(:,s);
    endfor
    
    MPA = struct();
    MPA.d = d;
    MPA.n = n;
    MPA.B = B;
    MPA.L = L;
endfunction

function p = mpa_psi_element(B, L, S)
    [n, d] = size(B);
    p = B{1, S(1)};
    for j = 2:n
        p = p * L{j-1} * B{j, S(j)};
    endfor
endfunction

function Psi, Mpa = test_mpa(d, n)
    psi = rand(d^n, 1);
    psi /= norm(psi);
    [MPA, e] = mpa(psi, d);
    
    S = randi(d, n, 1);
    i = psi_index(S, d);
    
    Psi = psi(i);
    Mpa = mpa_psi_element(MPA.B, MPA.L, S);
endfunction

function test_mpa_many()
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