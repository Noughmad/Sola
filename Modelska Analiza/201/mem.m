pkg load tsa

function [S,P] = mem(Podatki, K, F)
    N = max(size(Podatki));
    [AR,RC,PE] = durlev(acovf(Podatki, K));

    P = ar2poly(AR);

    x = 1:N;
    z = exp(-2*pi*i*x/N);
    S = 1 ./ abs(polyval(P,z)).^2;
    
    S = S .* mean(F ./ S, "a");
	S = S(1:N/2);
    P = ar2poly(AR);

	assert(abs(roots(P)) < 1);

endfunction

function [D,n] = podatki(P, trend)
	n = max(size(P));
	m = min(size(P));
	D = P;
	if (m > 1)
		D = D(m,:);
	endif
	if trend
		D = D(D > 0);
		n = max(size(D));
		co = polyfit(1:n, D, 1);
		D = D - polyval(co, 1:n);
		save g_co2_popravljen.dat D
	endif
endfunction

function naredi(Ime)
	[D, n] = podatki(dlmread([Ime ".dat"])', strcmp(Ime,"co2") || strcmp(Ime,"borza"));

    F = abs(fft(D)).^2;
	
	x = 1:n/2;
	
    semilogy(x, F(x), x, mem(D,30,F), x, mem(D,20,F), x, mem(D,10,F));
    legend("Podatki", "$m=30$", "$m=20$", "$m=10$");
    print("-depslatex", "-S640,480", ["g_" Ime "_psd.tex"])
endfunction

function L = loci(Frekvence, m)
	n = 512;
	x = 1:n;
	y = sum(sin(Frekvence' * x), 1);
	for f = Frekvence
		y = y + sin(f * x);
	endfor
	F = abs(fft(y)).^2;
	[S,P] = mem(y,m,F);
	E = [];
	for i = 2:255
		if ( (S(i-1) > S(i)) == (S(i+1) > S(i)) && S(i) > 10)
			E = [E S(i)];
		endif
	endfor
	L = (size(E,2) > 2);
endfunction

function plotroots(Ime)
	[D, n] = podatki(dlmread([Ime ".dat"])', strcmp(Ime,"co2") || strcmp(Ime,"borza"));
	K = 20;

	[S,P] = mem(D, K, abs(fft(D).^2));
	R = roots(P);

    x = 1:n;
    z = exp(-2*pi*i*x/n);

	plot(real(R), imag(R), "x", real(z), imag(z), "-");
	legend("Poli", "Enotska kro\"znica")
    print("-depslatex", "-S640,480", ["g_" Ime "_roots.tex"])
endfunction

function napoved(Ime)
	[D, n] = podatki(dlmread([Ime ".dat"])', strcmp(Ime,"co2") || strcmp(Ime,"borza"));
	h = floor(n/2)
	hD = D(1:h);
	m = min(size(D));
	F = abs(fft(hD).^2);
	K = [30 20 10];
	
	Napoved = [];
	for k = K
		[S,P] = mem(hD, k, F);
		P = poly2ar(P);
		N = D;
		for i = 1:h
			N(h+i) = N(h+i-1:-1:h+i-k) * P';
		endfor
		Napoved = [Napoved; N(h+1:n)];
	endfor

	
	M = n;
	x = h+1:M;
	I = 1:M-h;

	plot(1:M, D(1:M), x, Napoved(1,I), x, Napoved(2,I), x, Napoved(3,I));
    legend("Podatki", "$m=30$", "$m=20$", "$m=10$");
    print("-depslatex", "-S640,480", ["g_" Ime "_napoved.tex"])

	x = h-9:h+40;
	px = h+1:h+40;
	I = 1:40;
	plot(x, D(x), px, Napoved(1,I), px, Napoved(2,I), px, Napoved(3,I));
    legend("Podatki", "$m=30$", "$m=20$", "$m=10$");
    print("-depslatex", "-S640,480", ["g_" Ime "_napoved_zoom.tex"])
endfunction

function prva()
	naredi("val2");
	naredi("val3");
	naredi("co2");
	plotroots("co2");
endfunction

function napovedi()
	napoved("val2");
	napoved("val3");
	napoved("co2");
	napoved("luna");
	napoved("borza");
endfunction

function L = locljivosti()
	L = zeros(1,50);
	f1 = 0.2;
	for m = 5:50
		for f2 = linspace(f1, f1 + 0.6, 200)
			if (loci([f1, f2], m))
				L(m) = f2 - f1;
				break;
			endif
		endfor
	endfor

	x = 5:50;
	plot(x, L(x));
	xlabel("\"Stevilo polov");
	ylabel("Najmanj\"sa $\\Delta\\nu$")
    print("-depslatex", "-S640,480", ["g_loc.tex"]);
endfunction

prva()
napovedi()
naredi("borza")