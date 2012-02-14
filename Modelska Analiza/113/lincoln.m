global TaPrav

function [X, Y, V] = pgmread(s)
	fid = fopen(s, "r");
	fscanf(fid, "%c", 3);
	[V, count] = fscanf(fid, "%d");
	L = max(size(V));
	X = V(1);
	Y = V(2);
	X, Y, X*Y, count
	assert(X*Y == count-3)
	Max = V(3);
	V = V(4:L);
	V = V ./ Max;

	# Branje datoteke gre po vrsticah, zato jo najprej zlozimo v matriko enake dimenzije kot slika
	V = reshape(V, X, Y);

	# Sama slika pa je odbrana po stolpcih, zato jo transponiramo in zlozimo v en velik stolpec
	V = reshape(V', X*Y, 1);
	fclose(fid)
endfunction

function pgmwrite(s, V, X, Y)
	# TODO
	fid = fopen(s, "w");
	fprintf(fid, "P2\n%d %d\n%d\n", X, Y, 255);
	V = 255 .* V;
	V = reshape(V, Y, X);
	V = reshape(V', X*Y, 1);
	X, Y
	for i = 1:X*Y
		V(i) = max(min(V(i), 255), 0);
		if (mod(i, 19) != 0)
			fprintf(fid, "%d ", V(i));
		else
			fprintf(fid, "%d\n", V(i));
		endif
	endfor
	fclose(fid)
endfunction

function [X, Y, V] = linread(i)
	[X, Y, V] = pgmread(["lincoln_L30_N" int2str(i) "0.pgm"]);
endfunction

function Iz = lincoln(i, popravi)
	global TaPrav
	[X, Y, L] = linread(i);
	SizeOfL = size(L)

	# Odbrano je po stolpcih => zlozimo matriko v en velik stolpec
	N = X*Y;
	SL = reshape(L, N, 1);
	FL = fft(SL);

	tau = 30;

	R = fft(exp(-linspace(N, 0, N)/tau)/tau)';
	# R(R < eps) = 1;
	
	if (popravi)
		Iz = FL ./ R .* abs(TaPrav./FL).^2;
		pgmwrite(["g_lincoln_popravljen_" int2str(i) ".pgm"], abs(ifft(Iz)), X, Y);
	else
		if (i == 0)
			TaPrav = FL;
		endif
		Iz = FL ./ R;
		pgmwrite(["g_lincoln_" int2str(i) ".pgm"], abs(ifft(Iz)), X, Y);
	endif
endfunction

for i=0:4
	lincoln(i,0);
	lincoln(i,1);
endfor
