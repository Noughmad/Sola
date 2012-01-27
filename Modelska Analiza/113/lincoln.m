function [X, Y, V] = pgmread(s)
	V = dlmread(s)';
	L = max(size(V));
	X = V(2);
	Y = V(3);
	V = V(5:L) ./ V(4);
	V = reshape(reshape(V, X, Y)', X*Y, 1);
endfunction

function pgmwrite(s, V, X, Y)
	# TODO
endfunction

function [X, Y, V] = linread(i)
	[X, Y, V] = pgmread(["Lin" int2str(i) ".dat"]);
endfunction

function lincoln(i)

	[X, Y, L] = linread(i);
	SizeOfL = size(L)

	# Odbrano je po stolpcih => zlozimo matriko v en velik stolpec
	N = X*Y;
	SL = reshape(L, N, 1);
	FL = fft(SL);

	semilogy(linspace(0,1,N), abs(FL).^2);
	xlabel ("Frekvenca $f$")
	ylabel ("$|C(f)|^2$")
	print -depslatex g_fft_lincoln "-S640,480" 

	R = fft(exp(-linspace(0, N, N)/30))';
	R(R < eps) = 1;

	Iz = abs(ifft(FL ./ R));

	pgmwrite(["g_lin_" int2str(i) ".dat"], Iz, X, Y);
endfunction

lincoln(0)