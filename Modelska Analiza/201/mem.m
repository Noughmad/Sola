pkg load tsa

function S = mem(Podatki, K, F)
    N = max(size(Podatki));
    [AR,RC,PE] = durlev(acovf(Podatki, K));
    P = ar2poly(AR);
    R = roots(P);

    x = 1:N;
    z = exp(-2*pi*i*x/N);
    S = 1 ./ abs(polyval(P,z)).^2;
    
    S = S .* mean(F ./ S, "h");
	S = S(1:N/2);
endfunction

function naredi(Ime)
    D = dlmread([Ime ".dat"])';

	n = max(size(D));
	m = min(size(D));

	if (m > 1)
		D = D(m,:);
		co = polyfit(1:n, D, 1);
		D = D - polyval(co, 1:n);
	endif
    F = abs(fft(D)).^2;
	
	x = 1:n/2;
	
    semilogy(x, F(x), x, mem(D,30,F), x, mem(D,20,F), x, mem(D,10,F));
    legend("Podatki", "$m=30$", "$m=20$", "$m=10$");
    print("-depslatex", "-S640,480", ["g_" Ime "_psd.tex"])
endfunction

function plotroots(Ime)
	D = dlmread([Ime ".dat"])';
	K = 20;

	n = max(size(D));
	m = min(size(D));

	if (m > 1)
		D = D(m,:);
		co = polyfit(1:n, D, 1);
		D = D - polyval(co, 1:n);
	endif

    [AR,RC,PE] = durlev(acovf(D, K));
    P = ar2poly(AR);
    R = roots(P);

    x = 1:n;
    z = exp(-2*pi*i*x/n);

	plot(real(R), imag(R), "x", real(z), imag(z), "-");
	legend("Poli", "Enotska kro\"znica")
    print("-depslatex", "-S640,480", ["g_" Ime "_roots.tex"])
endfunction

naredi("val2");
naredi("val3");
naredi("co2");

plotroots("co2");