pkg load tsa

function [S,P] = mem(Podatki, K, F)
    N = max(size(Podatki));
    [AR,RC,PE] = durlev(acovf(Podatki, K));

    P = ar2poly(AR);

    x = 1:N;
    z = exp(-2*pi*i*x/N);
    S = 1 ./ abs(polyval(P,z)).^2;
    
    S = S .* mean(F ./ S, "h");
	S = S(1:N/2);
    P = ar2poly(AR);

endfunction

function [D,n] = podatki(P)
	n = max(size(P));
	m = min(size(P));
	D = P;
	if (m > 1)
		D = D(m,:);
		co = polyfit(1:n, D, 1);
		D = D - polyval(co, 1:n);
	endif
endfunction

function naredi(Ime)
	[D, n] = podatki(dlmread([Ime ".dat"])');

    F = abs(fft(D)).^2;
	
	x = 1:n/2;
	
    semilogy(x, F(x), x, mem(D,30,F), x, mem(D,20,F), x, mem(D,10,F));
    legend("Podatki", "$m=30$", "$m=20$", "$m=10$");
    print("-depslatex", "-S640,480", ["g_" Ime "_psd.tex"])
endfunction

function plotroots(Ime)
	[D, n] = podatki(dlmread([Ime ".dat"])');
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
	D = dlmread([Ime ".dat"])';
	n = max(size(D));
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
endfunction

naredi("val2");
naredi("val3");
naredi("co2");

plotroots("co2");

napoved("val2");
napoved("val3");
napoved("co2");
