
function n = st(N,i,j)
    n = (i-1)*N+j;
endfunction

function [i,j] = indeks(N,n)
    i = floor(n/N)+1;
    j = mod(n,N);
endfunction

function r = radij(N,i)
	r = (i-0.5)/N;

endfunction

function A = matrik_valj(nr,nfi)
	c = (pi*nr/nfi).^-2;
	D = ones(2/c);
	n = nr*nfi;
	A = zeros(n);
	for i=1:nr
		for j=1:nfi
			s = st(nfi,i,j);
			## Diagonalni clen
			A(s,s) = 2*c/radij(nr,i)^2 + 2;

			# Robni pogoji, vsi so prve vrste
			if i == 1
				A(s,s) = A(s,s) + 1;
			else
				A(s,st(nfi,i-1,j)) = -(1-0.5/i);
			endif

			if i == nr
				A(s,s) = A(s,s) + 1;
			else
				A(s,st(nfi,i+1,j)) = -(1+0.5/i);
			endif

			if j == 1
				A(s,s) = A(s,s) + c/radij(nr,i)^2;
			else
				A(s,st(nfi,i,j-1)) = -c/radij(nr,i)^2;
			endif

			if j == nfi
				A(s,s) = A(s,s) + c/radij(nr,i)^2;
			else
				A(s,st(nfi,i,j+1)) = -c/radij(nr,i)^2;
			endif
		endfor
	endfor
	A = sparse(A);
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

function izracun(n,c,stv)
    A = matrika(n);
    B = masa(n,c);
	SizeA = size(A)
	SizeB = size(B)
    [Vektorji, Vrednosti] = eigs(A,B,stv,'sm');
    for i = 1:st
        V = reshape(Vektorji(:,i),n,n);
        save(["g_opna_" int2str(n) "_" num2str(c) "_" int2str(stv-i+1) ".dat"], "V");
        Vrednosti(i,i)
    endfor
	M = reshape(full(diag(B)),n,n);
	save(["g_masa_" int2str(n) "_" num2str(c) ".dat"], "M");
endfunction

function izracun_valj(n,st)
	nr = n;
	nfi = 2*n;
    A = matrik_valj(nr,nfi);
    [Vektorji, Vrednosti] = eigs(A,st,'sm');
	diag(Vrednosti)
    for i = 1:st
        V = reshape(Vektorji(:,i),nfi,nr);
        save(["g_valj_" int2str(n) "_" int2str(st-i+1) ".dat"], "V");
    endfor
endfunction
