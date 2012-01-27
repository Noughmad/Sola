c = [];
for i = 0:3
	c = [c, dlmread(["signal" int2str(i) ".dat"])];
endfor

m = size(c, 1);

C = fft(c);
D = abs(C).^2;
x = linspace(0, 1, m)';

semilogy(x, D(:,1), x, D(:,2), x, D(:,3), x, D(:,4))
legend("signal0.dat", "signal1.dat", "signal2.dat", "signal3.dat")
xlabel("Frekvenca $f$")
ylabel("$|C(f)|^2$")
print -depslatex g_fft_c "-S640,480" 

function r = rob(x,y)
	r = zeros(x,y);
endfunction

N = [];
N(:,1) = zeros(m, 1);
N(:,2) = [0.01 * rob(50, 1); .99*D(51:462,1); 0.01*rob(50,1)];
N(:,3) = [rob(30, 1); .99*D(31:482,2); rob(30,1)];
N(:,4) = [100 * rob(20, 1); D(21:492,3); 100*rob(20,1)];

Phi = ones(size(D)) - N ./ D

r = exp(-abs(linspace(-255,256, m))/16)';
R = fft(r);
R(abs(R)<=eps) = 1;

s = abs(ifft(C .* Phi ./ [R R R R]));

clf
plot(x, s(:,1), x, s(:,2))
legend("signal0.dat", "signal1.dat", "signal2.dat", "signal3.dat")
xlabel("Cas $t$")
ylabel("$s(t)$")
print -depslatex g_signal "-S640,480" 

