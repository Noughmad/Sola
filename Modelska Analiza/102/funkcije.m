global A
global amax
global amin
global vmax
global N
global v0

function x = final_x(xx)
  x = [0;0;0;xx];
endfunction

function y = phi_t(xx)
  x = final_x(xx);
  global N
  y = 0;
  for i = 1:N
    for j = i+1:N
      ct = cos((x(2*i-1)-x(2*j-1))*2);
      sf = sin(x(2*i))*sin(x(2*j));
      cf = cos(x(2*i))*cos(x(2*j));
      y = y + 1/(1 - ct*sf - cf);
    endfor
  endfor
endfunction

N = 5;

A = diag(ones(N+1,1)) - diag(ones(N,1),-1);
A = A(2:N+1,1:N+1);
amax = 0.2;
amin = 1;
vmax = 2;
v0 = 2

function y = semafor(x)
  y = sumsq(pospesek(x));
endfunction

function a = pospesek(x)
  global A
  global v0
  a = A * [v0; x];
endfunction

function r = max_posp(x)
  global amax
  global amin
  p = pospesek(x);
  o = ones(size(p));
  r = [amax.*o-p; p-amin.*o]; 
endfunction

function r = pot(x)
  global N
  r = sum(x) - N;
endfunction

# Resitev za semafor
Y0 = linspace(0,1,N)';
Y = sqp(Y0, @semafor, @pot, [], 0, vmax);

# Resitev za naboje
x0 = linspace(0.2,pi-0.2,2*N-3)'; # Zagotovimo, da nobena dva nista na istem mestu
final_x(x0)./pi*180
X = sqp(x0, @phi_t, [], [], 0, pi);
R = final_x(X)./pi*180;
R = [2.*R(1:2:2*N) 90-R(2:2:2*N)]
save(['naboji_' int2str(N) '.dat'], 'R')