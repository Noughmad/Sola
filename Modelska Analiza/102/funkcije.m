global A
global amax
global amin
global vmax
global N
global v0
global B
global beta

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

N = 20;

amax = 0.2;
amin = -1;
vmax = 10;
v0 = 2

function y = semafor(x)
  y = sumsq(pospesek(x));
endfunction

function y = semafor_p(x)
  y = sumsq(pospesek_p(x));
endfunction

function y = semafor_a(x)
  global beta
  y = semafor(x) + exp( - beta * pot(x));
end

function y = semafor_a_p(x)
  global beta
  y = semafor_p(x) + exp( - beta * pot_p(x));
end

function a = pospesek(x)
  global A
  global v0
  global N
  a = A * [v0; x] .* N;
endfunction

function a = pospesek_p(x)
  global B
  global v0
  global N
  a = B * [v0; x; v0] .* (N+1);
endfunction

function r = max_posp(x)
  global amax
  global amin
  p = pospesek(x);
  o = ones(size(p));
  r = [amax.*o-p; p - amin.*o]; 
endfunction

function r = max_posp_p(x)
  global amax
  global amin
  p = pospesek_p(x);
  o = ones(size(p));
  r = [amax.*o-p; p - amin.*o]; 
endfunction

function r = pot(x)
  global N
  global v0
  r = v0/2 + sum(x) - x(N)/2 - N;
endfunction

function r = pot_p(x)
  global N
  global v0
  r = v0/2 + sum(x) + v0/2 - N - 1;
end

function r = posp_in_pot(x)
  r = [-pot(x); max_posp(x)];
end

function r = posp_in_pot_p(x)
  r = [-pot_p(x); max_posp_p(x)];
end

function s = semaforji()
  global N
  global vmax
  global A
  global v0
  global B
  global amax
  global amin
  global beta


  N = 20
  beta = 10
    
  A = diag(ones(N+1,1)) - diag(ones(N,1),-1);
  A = A(2:N+1,1:N+1);

  B = diag(ones(N+2,1)) - diag(ones(N+1,1),-1);
  B = B(2:N+2,1:N+2);

  amax = 3
  amin = -23
  
  for v0 = [0 0.5 1.5 2 3 4]
    for vmax = [9.9 1.3 1.2 1.1]

      Y0 = ones(N,1);
      y = sqp(Y0, @semafor, @pot, [], 0, vmax);
      Y = [ linspace(0,1,N+1)' ./ (pot(y)+1) [v0; y] ];
      saveTo = ['semafor_' int2str(10*v0) '_' int2str(10*vmax) '.dat']
      save(saveTo, 'Y')

      y = sqp(Y0, @semafor_p, @pot_p, [], 0, vmax);

      Y = [ linspace(0,1,N+2)' ./ (pot_p(y)+1) [v0; y; v0] ];
      saveTo = ['semafor_p_' int2str(10*v0) '_' int2str(10*vmax) '.dat']
      save(saveTo, 'Y')

      try
  
	Y0 = ones(N,1);
	y = sqp(Y0, @semafor_a, [], @posp_in_pot, 0, vmax);
	Y = [ linspace(0,1,N+1)' ./ (pot(y)+1) [v0; y] ];
	saveTo = ['semafor_a_' int2str(10*v0) '_' int2str(10*vmax) '.dat']
	save(saveTo, 'Y')
	disp(['Pot za ' saveTo ' je enaka ' pot(Y)])

      catch
	disp(['Error with ' int2str(10*v0) ' and ' int2str(10*vmax) '!'])
      end

      try
	Y0 = ones(N,1);
	y = sqp(Y0, @semafor_a_p, [], @posp_in_pot_p, 0, vmax);
	Y = [ linspace(0,1,N+2)' ./ (pot_p(y)+1) [v0; y; v0] ];
	saveTo = ['semafor_a_p_' int2str(10*v0) '_' int2str(10*vmax) '.dat']
	save(saveTo, 'Y')
      catch
	disp(['Error (periodic) with ' int2str(10*v0) ' and ' int2str(10*vmax) '!'])
      end

    endfor
  endfor
endfunction

# Resitev za semafor
semaforji()

# Resitev za naboje
%{
x0 = linspace(0.2,pi-0.2,2*N-3)'; # Zagotovimo, da nobena dva nista na istem mestu
final_x(x0)./pi*180
X = sqp(x0, @phi_t, [], [], 0, pi);
R = final_x(X)./pi*180;
R = [2.*R(1:2:2*N) 90-R(2:2:2*N)]
save(['naboji_' int2str(N) '.dat'], 'R')
%}
