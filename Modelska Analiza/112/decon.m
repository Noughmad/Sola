global n
global RF
global A

R = dlmread("signal0.dat");
RF = fft(R);
A = abs(RF);
n = size(R,1)

beta = 1e3;

function x = slabost(data)
    x = 0;
    for i = 1:size(data,1)-1
        x = x + abs(data(i) - data(i+1))^2;
    endfor
endfunction

function S = signal(beta)
    global n
    global RF
    global A

    x = linspace(-1,1,n)';
    G = fft(exp(-beta*abs(x)));
    SF = zeros(size(G));
    
    sigma = 1;
    SF(A>=sigma) = RF(A>=sigma) ./ G(A>=sigma);
    SF(A<sigma) = 0;
    
    S = ifft(SF);
endfunction

function S = sl_beta(beta)    
    S = slabost(signal(beta));
endfunction

Beta = fsolve("sl_beta", beta)
plot(linspace(0,1,n), signal(Beta), linspace(0,1,n), R/10)
xlabel ("$t$");
ylabel ("Signal");
legend ("Vhodni $S(t)$", "Izhodni $R(t)$")
orient landscape 
print -depslatex g_decon_signal "-S640,480"

