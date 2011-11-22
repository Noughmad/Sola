function [m, M, data] = normalize(D)
    m = min(D);
    M = max(D);
    if (m == M)
        M = m + 1;
    endif
    data = (D - m)/(M-m);
endfunction

## Preberemo podatke in jih reskaliramo, da bodo vsi med 0 in 1
Data = dlmread("prevodnost.dat");
[t0, t1, T] = normalize(Data(:,1));
[p0, p1, P] = normalize(Data(:,2));
[y0, y1, y] = normalize(Data(:,3));

## Ustvarimo matriko S da bodo nasi parametri resitve sistema S*a = y
N = 3;
n = N^2;
S = zeros(size(T) .* [1 n]);
for i = 0:N-1
    for j = 0:N-1
        k = j*N+i+1;
        S(:,k) = T.^i .* P.^j;
    endfor
endfor

a = S\y

[U, W, V] = svd(S);
w = diag(W)