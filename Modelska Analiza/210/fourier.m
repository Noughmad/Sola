function F = dst2(P)
    F = dst(dst(P)')';
endfunction

function F = idst2(P)
    F = idst(idst(P)')';
endfunction

function S = sum2(P)
    S = sum(sum(P,1),2);
endfunction

function R = dva(P)
    F = dst2(P);
    [m,n] = size(P);
    h = 1.0/m;
    k = 1.0/n;
    for i=1:m
        for j=1:n
            F(i,j) = h*k*F(i,j) / 2 / (cos(pi*i/m) + cos(pi*j/n) - 2);
        endfor
    endfor
    R = idst2(F);
endfunction

function R = ena(P)
    [m,n] = size(P);
    F = dst(P)';
    X = [];
    hk = 1.0 / m / n;
    for j=1:m
        D = (4-2*cos(j*pi/m)) * ones(n,1);
        A = sparse(-diag(D) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
        b = hk * F(:,j);
        X = [X, (A\b)];
    endfor
    R = idst(X');
endfunction

function R = valj(T, m, n)
    P = zeros(m,n);
    P(:,n) = T*m*n * ones(1,m);
    P(:,1) = T*m*n * ones(1,m);
    [m,n] = size(P);
    P
    F = dst(P)'
    X = [];
    for j=1:m
        D = (4-2*cos(j*pi/m)) * ones(n,1);
        r = linspace(-1,1,n);
        A = sparse(-diag(D) + diag(1 - 0.5 ./ n ./ r(1:n-1), 1) + diag(1 + 0.5 ./ n ./ r(2:n), -1));
        FullA = full(A)
        b = 1.0/m/n* F(:,j);
        X = [X, (A\b)];
    endfor
    R = idst(X');
endfunction
    
function plot_3d(U)
    [m,n] = size(U);
    x = 1:m;
    y = 1:n;
    [xx,yy] = meshgrid(x,y);
    mesh(xx, yy, U);
endfunction

function opna(G, size, name)
    U = dva(G)
    save(["g_opna_" name "_" int2str(size) "_2d.dat"], "U")
    V = ena(G)
    save(["g_opna_" name "_" int2str(size) "_1d.dat"], "V")
endfunction

function prazne_opne(Ns)
    for n=Ns
        opna(ones(n,n), n, "prazna");
    endfor
endfunction

function obtezene_opne(Ns, c)
    for n=Ns
        m = n/4;
        T = ones(m);
        G = [c*ones(n/2, n/4), ones(n/2, n/4), c*ones(n/2, n/2);  ones(n/4,n); c*ones(n/4, n/4), ones(n/4, n/4), c*ones(n/4, n/2) ];
        opna(G, n, ["obtezena_" num2str(c)])
    endfor
endfunction

function p = pretok(R)
    p = sum(sum(R,1,'extra'),2,'extra');
    [m,n] = size(R);
    p = abs(p) / m / n;
endfunction

function R = sor(G)
    [m,n] = size(G);
    # omega = 2.0/(1.0 + pi/n)
    popravek = 1;
    rho = -1.0/(m+1)/(n+1);
    meja = 1e-8;
    N = n*m;
    Res = zeros(N,1);
    D1 = [];
    for i=1:n
        D1 = [D1; ones(n-1,1); 0];
    endfor
    D1 = D1(1:N-1);
    
    D2 = [ones(N-n,1)];
    Dd = -4*ones(N,1);
    M = spdiags([[D2; zeros(n,1)], [D1;0], Dd, [0;D1], [zeros(n,1); D2]], [-n -1 0 1 n], N, N);
    while popravek > meja
        P = M*Res + rho;
        popravek = sumsq(P);
        Res = Res + P*0.25;
    endwhile
    R = reshape(Res,m,n);
endfunction

function cajt()
	R = [];
	for n=[16 32 64 128 256 512 1024 2048 4096]
		T = [n];
		G = ones(n,n);
        
		a = time();
		E = ena(G);
                d = time()-a;
		T = [T pretok(E) d];
	
                a = time();
                D = dva(G);
                d = time()-a;
                T = [T pretok(D) d];
                
                if n < 1000
                    a = time();
                    S = sor(G);
                    d = time()-a;
                    T = [T pretok(S) d];
                else
                    T = [T 0 0];
                endif
                
                T
                R = [R; T];
	endfor
	save g_cas.dat R
endfunction

# obtezene_opne([16 32 64 128 256 512 1024], 3);

# cajt()