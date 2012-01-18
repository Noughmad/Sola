function c = conv_fft(data, n)
  ds = size(data,2);
  data(n*ds) = 0;
  c = ifft(fft(data).^n);
endfunction

function c = conv_normal(data, n)
  c = data;
  for i = 1:(n-1)
    c = [conv(c, data), 0];
  endfor
endfunction

R = [];
N = 1

for i = 1:18
  N = 2*N;
  a = 2/N;
  data = a * ( 1 - linspace(1, 0, N));

  [t,u,s] = cputime;
  T = [t u s];
  C = conv_fft(data, 3);
  [t,u,s] = cputime;
  T1 = [t u s];
  Tfft = T1 - T;

  [t,u,s] = cputime;
  T = [t u s];
  D = conv_normal(data, 3);
  [t,u,s] = cputime;
  T1 = [t u s];
  Tnorm = T1 - T;
  
  R = [R; N Tfft Tnorm sumsq(D-C)];
  disp(["Solved for " int2str(N)]);
endfor

save g_konvolucija_hitrost.dat R
