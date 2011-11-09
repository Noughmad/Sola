function semafor_lin(N, vmax, amax, amin, y0, utez)
  if utez == 1
    c = linspace(0,1,N)';
  else
    c = [ zeros(N-1,1); 1];
  endif

  O = ones(N,1);

  Aa = diag(ones(N,1)) - diag(ones(N-1,1),-1);
  bvu = amax * O;
  bvu(1) = bvu(1) + y0;
  bvl = amin * O;
  bvl(1) = bvl(1) + y0;
  
  As = ones(1,N);
  bs = N;

  A = [ Aa; Aa; As ];
  b = [ bvu; bvl; bs ];

  ctype = [];
  for i = 1:N
    ctype = [ctype "U"];
  end

  for i = 1:N
    ctype = [ctype "L"];
  end

  ctype = [ctype "S"];

  vartype = [];
  for i = 1:N
    vartype = [ vartype "C" ];
  endfor

  [X, F, status, extra] = glpk(c, A, b, zeros(100,1), vmax * ones(100,1), ctype, vartype, -1);

  saveTo = [ "lin_" int2str(10*y0) "_" int2str(utez) ".dat" ]
  Y = [linspace(0,1,N+1)' [y0; X ] ];
  save(saveTo, "Y");
end

for y0 = [0 5 9 15 20]
  semafor_lin(100, 2, 0.03, -0.05, y0/10, 1)
  semafor_lin(100, 2, 0.03, -0.05, y0/10, 0)
endfor