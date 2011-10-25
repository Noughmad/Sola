
# Resitev za naboje


function x = final_x(xx)
  x = [0;0;0;xx];
endfunction

function y = phi_t(xx)
  x = final_x(xx);
  N = max(size(x))/2;
  y = 0;
  for i = 1:N
    for j = i+1:N
      ct = cos((x(2*i-1)-x(2*j-1))*2);
      sf = sin(x(2*i))*sin(x(2*j));
      cf = cos(x(2*i))*cos(x(2*j));
      if ( ct*sf + cf == 1 )
	y = inf;
	return
      else
	y = y + 1/(1 - ct*sf - cf);
      endif
    endfor
  endfor
endfunction

function razporeditev(N)
  x0 = rand(2*N-3,1) * pi;
  X = sqp(x0, @phi_t, [], [], 0, pi);
  R = final_x(X)./pi*180;
  R = [2.*R(1:2:2*N) 90-R(2:2:2*N)];
  saveTo = ['naboji_' int2str(N) '.dat']
  save(saveTo, 'R')
endfunction

for st = 2:12
  razporeditev(st)
endfor