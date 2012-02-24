pkg load tsa

function [P,R] = mem(Podatki, K)
    N = max(size(Podatki));
    [AR,RC,PE] = durlev(acovf(Podatki, K));
    P = ar2poly(AR);
    R = roots(P);
endfunction

function naredi(Ime, K)
    D = dlmread([Ime ".dat"])';
    [P,R] = mem(D, K);
    
    n = max(size(D))
    x = 1:n;
    z = exp(-2*pi*i*x/n);
    S = 1 ./ abs(polyval(P,z)).^2;
    
    F = abs(fft(D)).^2;
    
    S = S .* mean(F ./ S, "g");
    
    plot(x, S, x, F);
    legend("Priblizek", "Podatki");
    print("-depslatex", "S480,320", ["g_" Ime "_psd.tex"])
endfunction

naredi("val2", 20);
naredi("val3", 20);
