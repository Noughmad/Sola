global n
global RF
global A
global utez

R = dlmread("signal0.dat");
RF = fft(R);
A = abs(RF);
n = size(R,1)

beta = 1;

utez = 1 - abs(linspace(-1,1,n));

function x = slabost(data)
    global utez
    x = utez * abs(data);
    size(x)
endfunction

function x = frekvenca(data)
    x = sumsq(abs(data));
endfunction

function x = sprem(data)
    global n
    d = ifft(data);
    x = 0;
    for i = 1:(n-1)
        x = x + abs(d(i) - d(i+1))^2;
    endfor
endfunction

function x = amplituda(data)
    global n
    [v, omega] = max(data);
    d = abs(ifft(data));
    x = 0;
    for i = 1:(n-1)
        x = x + omega^2 * d(i)^2 + (d(i+1) - d(i))^2;
    endfor
endfunction

function S = ftsignal(beta)
    global n
    global RF
    global A

    x = linspace(-1,1,n)';
    G = fft(exp(-beta*abs(x)));
    SF = zeros(size(G));
    
    sigma = 1;
    SF(A>=sigma) = RF(A>=sigma) ./ G(A>=sigma);
    SF(A<sigma) = 0;
    
    S = SF;
endfunction

function S = sl_beta(beta)    
    S = slabost(ftsignal(beta));
endfunction

function S = am_beta(beta)    
    S = amplituda(ftsignal(beta));
endfunction

function S = sprem_beta(beta)    
    S = sprem(ftsignal(beta));
endfunction

function S = frek_beta(beta)    
    S = frekvenca(ftsignal(beta));
endfunction

Beta1 = fsolve("sl_beta", beta)
Beta2 = fsolve("sprem_beta", beta)
Beta3 = fsolve("am_beta", beta)

x = linspace(-1,1,n);
plot(x, R/10, x, ifft(ftsignal(Beta1)), x, ifft(ftsignal(Beta2)), x, ifft(ftsignal(Beta3)))
xlabel ("$t$");
ylabel ("Signal");
legend ("Izhodni $R(t)$", "Sprememba", "Frekvenca", "Amplituda")
orient landscape 
print -depslatex g_decon_signal "-S640,480"


B = linspace(0,5);
Z = [];
for b = B
    Z = [Z real(ifft(ftsignal(exp(b))))];
endfor
[xx, bb] = meshgrid(x, B);
mesh(xx,bb,Z')

xlabel ("$t$")
ylabel ("$\\log \\beta$")
zlabel ("Signal")
print -depslatex g_decon_3d "-S800,600"
