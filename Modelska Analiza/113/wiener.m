global OrigSignal;
global R;
global CurrentC;
global OrigS;

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

function [Sig, Phi] = dekonsignal(C, Phi, popravi)
	global R
	global CurrentC
	D = abs(C).^2;
	if (popravi)
		CurrentC = C;
		Phi = sqp(Phi, @napaka_suma, [], []);
	endif
	Sig = abs(ifft(C .* Phi ./ R));
endfunction

function n = napaka_suma(Phi)
	global R;
	global OrigSignal;
	global CurrentC;
	s = abs(ifft(CurrentC .* Phi ./ R));
	n = sumsq(s - OrigSignal);
endfunction

N = [];
N(:,1) = zeros(m, 1);
N(:,2) = [0.01 * rob(50, 1); D(51:462,1); 0.01*rob(50,1)];
N(:,3) = [rob(30, 1); D(31:482,2); rob(30,1)];
N(:,4) = [100 * rob(20, 1); D(21:492,3); 100*rob(20,1)];

r = exp(-abs(linspace(-255,256, m))/16)';
R = fft(r);
R(abs(R)<=eps) = 1;


[s, p] = dekonsignal(C(:,1), ones(m,1), 0);
OrigSignal = s;

for i = 2:4
	[ss, pp] = dekonsignal(C(:,i), abs( C(:,1) ./ C(:,i) ).^2, 0);
	s = [s, ss];
	p = [p, pp];
endfor

Splosno = dekonsignal(C(:,1), ones(m,1), 0);
CR = conj(R);
for i = 2:4
	Splosno = [Splosno, dekonsignal(C(:,i), ones(m,1) - N(:,i) ./ D(:,i), 0)];
endfor

clf
plot(x, s(:,1), x, s(:,2))
legend("signal0.dat", "signal1.dat", "signal2.dat", "signal3.dat")
xlabel("Cas $t$")
ylabel("$s(t)$")
print -depslatex g_signal "-S640,480" 

save g_wiener_input.dat D
save g_wiener_phi.dat p
save g_wiener_rezultat.dat s
save g_wiener_splosno.dat Splosno

