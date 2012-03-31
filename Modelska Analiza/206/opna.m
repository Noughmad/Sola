
function n = st(N,i,j)
    n = (i-1)*N+j;
endfunction

function [i,j] = indeks(N,n)
    i = floor(n/N)+1;
    j = mod(n,N);
endfunction

function A = matrika(n)
    # Stirice na diagonali
    A = 4*speye(n*n);
    
    # Prispevki vseh stirih sosedov zaradi Nabla^2
    ob1 = mod(1:n*n-1,n) != 0;
    A = A - diag(ob1,1) - diag(ob1,-1);
    A = A - diag(ones(n*n-n,1),n) - diag(ones(n*n-n,1),-n);
    
    # Robni pogoj: Na vseh stirih straneh je dolocena vrednost => po 1 pristejemo
    for i = 1:n
        in = st(n,i,1);
        A(in,in) = A(in,in) + 1;
        in = st(n,i,n);
        A(in,in) = A(in,in) + 1;
        in = st(n,1,i);
        A(in,in) = A(in,in) + 1;
        in = st(n,n,i);
        A(in,in) = A(in,in) + 1;
    endfor
endfunction

function B = masa(n, c)
    m = n/4;
    T = ones(m);
    B = [c*ones(n/2, n/4), ones(n/2, n/4), c*ones(n/2, n/2);  ones(n/4,n); c*ones(n/4, n/4), ones(n/4, n/4), c*ones(n/4, n/2) ];
    B = sparse(diag(reshape(B,n*n,1)));
endfunction

function izracun(n,c,st)
    A = matrika(n);
    B = masa(n,c);
    [Vektorji, Vrednosti] = eig(A,B);
    for i = 1:st
        V = reshape(Vektorji(:,i),n,n);
        save(["g_opna_" int2str(n) "_" num2str(c) "_" int2str(i) ".dat"], "V");
        Vrednosti(i,i)
    endfor
endfunction

izracun(32, 2, 10);
izracun(32, 20, 10);