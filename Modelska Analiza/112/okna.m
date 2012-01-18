val2 = dlmread("val2.dat");
val3 = dlmread("val3.dat");
N = max(size(val2));

orient landscape 
x = linspace(0,1,N);
semilogy(x, abs(fft(val2)), x, abs(fft(val3)))
legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
xlabel("$\\omega$")
ylabel("$S(\\omega)$")

print -depslatex g_val_fft "-S640,480"

hamm2 = hamming(N) .* val2;
hamm3 = hamming(N) .* val3;

semilogy(x, abs(fft(hamm2)), x, abs(fft(hamm3)))
legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
xlabel("$\\omega$")
ylabel("$S(\\omega)$")
print -depslatex g_val_hamm "-S640,480" 

hamm2 = sin(2*pi*x)' .* val2;
hamm3 = sin(2*pi*x)' .* val3;

semilogy(x, abs(fft(hamm2)), x, abs(fft(hamm3)))
legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
xlabel("$\\omega$")
ylabel("$S(\\omega)$")
print -depslatex g_val_cos "-S640,480" 

hamm2 = blackman(N) .* val2;
hamm3 = blackman(N) .* val3;

semilogy(x, abs(fft(hamm2)), x, abs(fft(hamm3)))
legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
xlabel("$\\omega$")
ylabel("$S(\\omega)$")
print -depslatex g_val_black "-S640,480" 

for M =  [64 128 256 512]
    x = linspace(0,1,M);
    semilogy(x, abs(fft(val2(1:M))), x, abs(fft(val3(1:M))))
    legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
    axis("nolabel");
    print("-depslatex", ["g_val_" int2str(M)], "-S440,360")

    H = hamming(M);
    semilogy(x, abs(fft(H .* val2(1:M))), x, abs(fft(H .* val3(1:M))))
    legend("\\texttt{val2.dat}", "\\texttt{val3.dat}")
    axis("nolabel");
    print("-depslatex", ["g_val_hamm_" int2str(M)], "-S440,360")
endfor