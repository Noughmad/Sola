function F = dst2(P)
    F = dst(dst(P)')';
endfunction

function F = idst2(P)
    F = idst(idst(P)')';
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
    F = dst(P)
    
    F = F'
    X = [];
    for j=1:m
        D = (4-2*cos(j*pi/m)) * ones(n,1);
        A = sparse(-diag(D) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
        b = 1.0/m/n* F(:,j)
        X = [X, (A\b)];
    endfor
    X = X'
    # plot_3d(X)
    R = idst(X)
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

function cajt()
	R = [];
	for n=[16 32 64 128]
		T = [n];
		G = ones(n,n);
		a = time();
		ena(G);
		T = [T (time()-a)];
		a = time();
		dva(G);
		T = [T (time()-a)];
		R = [R; T];
	endfor
	save g_cas.dat R
endfunction

# obtezene_opne([16 32 64 128 256 512 1024], 3);

cajt()